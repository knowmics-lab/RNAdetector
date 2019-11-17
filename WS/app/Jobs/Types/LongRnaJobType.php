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
            'paired'                   => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'           => 'Required, input file for the analysis. FASTQ or BAM',
            'secondInputFile'          => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'                => 'Required, type of the input file (fastq, bam)',
            'convertBam'               => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'               => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'customFASTATranscriptome' => 'An optional Transcriptome to employ for custom annotation (Not needed for human hg19 mRNAs and lncRNAs)',
            'customTranscriptomeName'  => 'An optional name for the custom transcriptome',
            'threads'                  => 'Number of threads for this analysis (Default 1)',
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
            'paired'                   => ['filled', 'boolean'],
            'firstInputFile'           => ['required', 'string'],
            'secondInputFile'          => [
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
            'inputType'                => ['required', Rule::in(self::VALID_INPUT_TYPES)],
            'convertBam'               => ['filled', 'boolean'],
            'trimGalore'               => ['filled', 'array'],
            'trimGalore.enable'        => ['filled', 'boolean'],
            'trimGalore.quality'       => ['filled', 'integer'],
            'trimGalore.length'        => ['filled', 'integer'],
            'customFASTATranscriptome' => ['filled', 'string'],
            'customTranscriptomeName'  => ['filled', 'alpha_num'],
            'threads'                  => ['filled', 'integer'],
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
        $customFASTATranscriptome = $this->model->getParameter('customFASTATranscriptome');
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
     * @param bool        $paired
     * @param string      $firstInputFile
     * @param string|null $secondInputFile
     * @param string      $inputType
     * @param int         $threads
     * @param string|null $customFASTATranscriptome
     * @param string|null $customTranscriptomeName
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runSalmon(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        string $inputType,
        int $threads = 1,
        ?string $customFASTATranscriptome = null,
        ?string $customTranscriptomeName = null
    ): array {
        $index = false;
        if ($customFASTATranscriptome === null) {
            $transcriptomeDir = env('HUMAN_SALMON_INDEXED_TRANSCRIPTOME');
            $this->log('Using Human mRNAs\lncRNAs indexed transcriptome.');
        } else {
            $transcriptomeDir = env('CUSTOM_TRANSCRIPTOME_PATH') . '/' . $customTranscriptomeName;
            if (!file_exists($transcriptomeDir) || !is_dir($transcriptomeDir)) {
                $index = true;
            }
        }
        if ($index) {
            $this->log('Indexing custom transcriptome');
            // Call Salmon index script
            $command = [
                'bash',
                env('BASH_SCRIPT_PATH') . '/salmon_index_2.sh',
                '-r',
                $customFASTATranscriptome,
                '-i',
                $transcriptomeDir
            ];
            $output = AbstractJob::runCommand(
                $command,
                $this->model->getAbsoluteJobDirectory(),
                null,
                null,
                [
                    3 => 'FASTA file with transcripts does not exist.',
                    4 => 'Indexed trascriptome folder does not exist.',
                ]
            );
            $this->log('Custom transcriptome indexed');
            if (!file_exists($transcriptomeDir) && !is_dir($transcriptomeDir)) {
                throw new ProcessingJobException('Unable to create indexed transcriptome');
            }
            $this->log($output);
            return $transcriptomeDir;
        }
        $salmonOutputRelative = $this->model->getJobTempFile('salmon_output', '.txt');
        $salmonOutput = $this->model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = \Storage::disk('public')->url($salmonOutput);
        switch ($inputType) {
            case self::FASTQ:
                // TODO: Call salmon counting on fastq
                $command = [
                    'bash',
                    env('BASH_SCRIPT_PATH') . '/salmon_counting.sh',
                    '-i',
                    $transcriptomeDir,
                    '-f',
                    $firstInputFile,
                    '-t',
                    $threads,
                    '-o',
                    $salmonOutput
                ];
                if ($paired) {
                    $command[] = '-s';
                    $command[] = $secondInputFile;
                }
                $output = AbstractJob::runCommand(
                    $command,
                    $this->model->getAbsoluteJobDirectory(),
                    null,
                    null,
                    [
                        3 => 'Input file does not exist.',
                        4 => 'Second input file does not exist.',
                        5 => 'Output file must be specified.',
                        6 => 'Output directory is not writable.',
                        7 => 'Indexed trascriptome does not exist.',
                        8 => 'Unable to find output file.',
                    ]
                );
                break;
            case self::BAM:
                // TODO: Call salmon counting on bam
                $command = [
                    'bash',
                    env('BASH_SCRIPT_PATH') . '/salmon_counting_bam.sh',
                    '-r',
                    $customFASTATranscriptome,
                    '-i',
                    $firstInputFile,
                    '-t',
                    $threads,
                    '-o',
                    $salmonOutput
                ];
                $output = AbstractJob::runCommand(
                    $command,
                    $this->model->getAbsoluteJobDirectory(),
                    null,
                    null,
                    [
                        3 => 'Input file does not exist.',
                        4 => 'FASTA transcriptome file does not exist.',
                        5 => 'Output directory must be specified.',
                        6 => 'Output directory is not writable.',
                        7 => 'Unable to find output file.',
                    ]
                );
                break;
            default:
                throw new ProcessingJobException('Unsupported input type');
        }
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $this->log($output);
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
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq($this->model, $paired, $firstInputFile);
            $inputType = self::FASTQ;
        }
        [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
        if ($inputType === self::FASTQ && $trimGaloreEnable) {
            [$firstTrimmedFastq, $secondTrimmedFastq] = self::runTrimGalore(
                $this->model,
                $paired,
                $firstInputFile,
                $secondInputFile,
                $trimGaloreQuality,
                $trimGaloreLength
            );
        }
        [$salmonOutput, $salmonOutputUrl] = $this->runSalmon(
            $paired,
            $firstTrimmedFastq,
            $secondTrimmedFastq,
            $inputType,
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
