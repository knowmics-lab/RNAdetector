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

class CircRnaJobType extends AbstractJob
{

    use ConvertsBamToFastqTrait, RunTrimGaloreTrait;

    private const FASTQ             = 'fastq';
    private const BAM               = 'BAM';
    private const SAM               = 'SAM';
    private const VALID_INPUT_TYPES = [self::FASTQ, self::BAM, self::SAM];

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'paired'               => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'       => 'Required, input file for the analysis',
            'secondInputFile'      => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'            => 'Required, type of the input file (fastq, bam, sam)',
            'convertBam'           => 'If inputType is bam converts input in another format: fastq or sam.',
            'trimGalore'           => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'customGTFFile'        => 'An optional GTF file for custom annotation of reads (Not needed for human hg19)',
            'customGenomeName'     => 'An optional name for the custom genome',
            'customFASTAGenome'    => 'An optional Genome to employ for custom annotation (Not needed for human hg19)',
            'threads'              => 'Number of threads for this analysis (Default 1)',
            'ciriSpanningDistance' => 'The maximum spanning distance used in CIRI (Default 500000)',
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
            'inputType'            => ['required', Rule::in(self::VALID_INPUT_TYPES)],
            'convertBam'           => ['filled', 'boolean'],
            'trimGalore'           => ['filled', 'array'],
            'trimGalore.enable'    => ['filled', 'boolean'],
            'trimGalore.quality'   => ['filled', 'integer'],
            'trimGalore.length'    => ['filled', 'integer'],
            'customGTFFile'        => ['filled', 'string'],
            'customGenomeName'     => ['filled', 'alpha_num'],
            'customFASTAGenome'    => ['filled', 'string'],
            'threads'              => ['filled', 'integer'],
            'ciriSpanningDistance' => ['filled', 'integer'],
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
        $customGTFFile = $this->model->getParameter('customGTFFile');
        $customFASTAGenome = $this->model->getParameter('customFASTAGenome');
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
        if (!empty($customGTFFile) && !$disk->exists($dir . $customGTFFile)) {
            return false;
        }
        if (!empty($customFASTAGenome) && !$disk->exists($dir . $customFASTAGenome)) {
            return false;
        }

        return true;
    }

    /**
     * Runs BWA analysis
     *
     * @param bool        $paired
     * @param string      $firstInputFile
     * @param string|null $secondInputFile
     * @param int         $threads
     * @param string|null $customGTFFile
     * @param string|null $customFASTAGenome
     * @param string|null $customGenomeName
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runBWA(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        int $threads = 1,
        ?string $customGTFFile = null,
        ?string $customFASTAGenome = null,
        ?string $customGenomeName = null
    ): string {
        if ($customGTFFile === null) {
            $customGTFFile = env('HUMAN_GTF_PATH');
        }
        if ($customFASTAGenome === null) {
            $genomeDir = env('HUMAN_BWA_GENOME');
            $index = false;
        } else {
            $genomeDir = env('CUSTOM_GENOME_PATH') . '/' . $customGenomeName;
            if (!file_exists($genomeDir) || !is_dir($genomeDir)) {
                $index = true;
            }
        }
        if ($index) {
            // Call BWA index script
            if (!file_exists($genomeDir) && !is_dir($genomeDir)) {
                throw new ProcessingJobException('Unable to create indexed genome');
            }
        }
        $samOutput = $this->model->getJobTempFileAbsolute('bwa_output', '.sam');
        // Call BWA
        if (!file_exists($samOutput)) {
            throw new ProcessingJobException('Unable to create BWA output file');
        }

        return $samOutput;
    }

    /**
     * Runs CIRI analysis
     *
     * @param bool        $paired
     * @param string      $ciriInputFile
     * @param string|null $customGTFFile
     * @param string|null $customFASTAGenome
     * @param int         $spanningDistance
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runCIRI(
        bool $paired,
        string $ciriInputFile,
        ?string $customGTFFile = null,
        ?string $customFASTAGenome = null,
        int $spanningDistance = 500000
    ): array {
        if ($customGTFFile === null) {
            $customGTFFile = env('HUMAN_GTF_PATH');
        }
        if ($customFASTAGenome === null) {
            $customFASTAGenome = env('HUMAN_FASTA_PATH');
        }
        $ciriOutputRelative = $this->model->getJobTempFile('ciri_output', '.txt');
        $ciriOutput = $this->model->absoluteJobPath($ciriOutputRelative);
        $ciriOutputUrl = \Storage::disk('public')->url($ciriOutput);
        // Call CIRI
        if (!file_exists($ciriOutput)) {
            throw new ProcessingJobException('Unable to create CIRI output file');
        }

        return [$ciriOutputRelative, $ciriOutputUrl];
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
        $customGTFFile = $this->model->getParameter('customGTFFile');
        $customFASTAGenome = $this->model->getParameter('customFASTAGenome');
        $customGenomeName = $this->model->getParameter(
            'customGenomeName',
            ($customFASTAGenome !== null) ? pathinfo($customFASTAGenome, PATHINFO_FILENAME) : null
        );
        $threads = (int)$this->model->getParameter('threads', 1);
        $ciriSpanningDistance = (int)$this->model->getParameter('ciriSpanningDistance', 500000);
        $ciriInputFile = null;
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
            $ciriInputFile = $this->runBWA(
                $paired,
                $firstTrimmedFastq,
                $secondTrimmedFastq,
                $threads,
                $customGTFFile,
                $customFASTAGenome,
                $customGenomeName
            );
        } elseif ($inputType === self::BAM) {
            $ciriInputFile = $this->model->getJobTempFileAbsolute('bam2sam', '.sam');
            // Call conversion script
            if (!file_exists($ciriInputFile)) {
                throw new ProcessingJobException('Unable to create converted BAM file');
            }
        } else {
            $ciriInputFile = $firstInputFile;
        }
        [$ciriOutput, $ciriOutputUrl] = $this->runCIRI(
            $paired,
            $ciriInputFile,
            $customGTFFile,
            $customFASTAGenome,
            $ciriSpanningDistance
        );
        $this->model->setOutput(
            [
                'outputFile' => [
                    'path' => $ciriOutput,
                    'url'  => $ciriOutputUrl,
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
        return 'Runs circRNA analysis from sequencing data';
    }
}
