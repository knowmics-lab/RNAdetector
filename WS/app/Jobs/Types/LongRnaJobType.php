<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsBamToFastqTrait;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Utils;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class LongRnaJobType extends AbstractJob
{
    use ConvertsBamToFastqTrait, RunTrimGaloreTrait;

    private const FASTQ             = 'fastq';
    private const BAM               = 'BAM';
    private const VALID_INPUT_TYPES = [self::FASTQ, self::BAM];
    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'paired'               => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'       => 'Required, input file for the analysis. FASTQ or BAM',
            'secondInputFile'      => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'            => 'Required, type of the input file (fastq, bam)',
            'convertBam'           => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'           => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'customFASTATranscriptome'    => 'An optional Transcriptome to employ for custom annotation (Not needed for human hg19 mRNAs and lncRNAs)',
            'customTranscriptomeName'     => 'An optional name for the custom transcriptome',
            'threads'              => 'Number of threads for this analysis (Default 1)',
        ];
    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [
            'outputFile' => 'Formatted read counts files (If multiple files a zip archive is returned)',
        ];
    }

    /**
     * Returns an array containing rules for input validation.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return array
     */
    public static function validationSpec(Request $request): array
    {
        return [
            'paired'               => ['filled', 'boolean'],
            'firstInputFile'       => ['required', 'string'],
            'secondInputFile'      => [
                Rule::requiredIf(
                    static function () use ($request) {
                        return $request->get('parameters.inputType') === self::FASTQ && ((bool)$request->get(
                                'parameters.paired',
                                false
                            )) === true;
                    }
                ),
                'string',
            ],
            'inputType'                 => ['required', Rule::in(self::VALID_INPUT_TYPES)],
            'convertBam'                => ['filled', 'boolean'],
            'trimGalore'                => ['filled', 'array'],
            'trimGalore.enable'         => ['filled', 'boolean'],
            'trimGalore.quality'        => ['filled', 'integer'],
            'trimGalore.length'         => ['filled', 'integer'],
            'customFASTATranscriptome'  => ['filled', 'string'],
            'customTranscriptomeName'   => ['filled', 'alpha_num'],
            'threads'                   => ['filled', 'integer'],
        ];
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $customFASTATranscriptome= $this->model->getParameter('customFASTATranscriptome');
        if (!in_array($inputType, self::VALID_INPUT_TYPES, true)) {
            return false;
        }
        $disk = Storage::disk('public');
        $dir = $this->model->getJobDirectory() . '/';
        if (!$disk->exists($dir . $firstInputFile)) {
            return false;
        }
        if ($paired && $inputType === self::FASTQ && (empty($secondInputFile) || !$disk->exists(
                    $dir . $secondInputFile
                ))) {
            return false;
        }
        if (!empty($customFASTATranscriptome) && !$disk->exists($dir . $customFASTATranscriptome)) {
            return false;
        }
        return true;
    }

    /**
     * Run Salmon analysis
     *
     */

    private function runSalmon(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        int $threads = 1,
        ?string $customFASTATranscriptome =null,
        ?string $customTranscriptomeName = null
    ): array {
        if ($customFASTATranscriptome === null) {
            $transcriptomeDir = env('HUMAN_SALMON_INDEXED_TRANSCRIPTOME');
            $index = false;
        } else {
            $transcriptomeDir = env('CUSTOM_TRANSCRIPTOME_PATH') . '/' . $customTranscriptomeName;
            if (!file_exists($transcriptomeDir) || !is_dir($transcriptomeDir)) {
                $index = true;
            }
        }
        if ($index) {
            // Call Salmon index script
            if (!file_exists($transcriptomeDir) && !is_dir($transcriptomeDir)) {
                throw new ProcessingJobException('Unable to create indexed transcriptome');
            }
        }
        $salmonOutputRelative = $this->model->getJobTempFile('salmon_output', '.txt');
        $salmonOutput = $this->model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = \Storage::disk('public')->url($salmonOutput);
        // Call Salmon counting scripts
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        return [$salmonOutputRelative, $salmonOutputUrl];
    }

    /**
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    public function handle(): void
    {
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $convertBam = (bool)$this->model->getParameter('convertBam', false);
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->model->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->model->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 14);
        $customFASTATranscriptome = $this->model->getParameter('customFASTATranscriptome');
        $customTranscriptomeName = $this->model->getParameter(
            'customTranscriptomeName',
            ($customFASTATranscriptome !== null) ? pathinfo($customFASTATranscriptome, PATHINFO_FILENAME) : null
        );
        $threads = (int)$this->model->getParameter('threads', 1);
        if ($inputType === self::BAM && $convertBam) {
            $inputType = self::FASTQ;
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq($this->model, $paired, $firstInputFile);
        }
        if ($inputType === self::FASTQ) {
            if ($trimGaloreEnable) {
                [$firstTrimmedFastq, $secondTrimmedFastq] = self::runTrimGalore(
                    $this->model,
                    $paired,
                    $firstInputFile,
                    $secondInputFile,
                    $trimGaloreQuality,
                    $trimGaloreLength
                );
            } else {
                [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
            }
        }
        [$salmonOutput,$salmonOutputUrl] = $this->runSalmon(
            $paired,
            $firstInputFile,
            $secondInputFile,
            $threads,
            $customFASTATranscriptome,
            $customTranscriptomeName
        );
        $this->model->setOutput(
            [
                'outputFile' => [
                    'path' => $salmonOutput,
                    'url'  => $salmonOutputUrl,
                ],
            ]
        );
        $this->model->save();
    }


    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Runs mRNAs and\or lncRNAs (long RNAs reads) analysis from sequencing data';
    }
}
