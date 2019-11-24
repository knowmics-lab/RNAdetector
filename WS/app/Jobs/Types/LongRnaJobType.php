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
use App\Models\Reference;
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
            'paired'          => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'  => 'Required, input file for the analysis. FASTQ or BAM',
            'secondInputFile' => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'       => 'Required, type of the input file (fastq, bam)',
            'convertBam'      => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'      => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'transcriptome'   => 'An optional transcriptome to employ for annotation (Not needed for human hg19 mRNAs and lncRNAs)',
            'threads'         => 'Number of threads for this analysis (Default 1)',
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
            'paired'             => ['filled', 'boolean'],
            'firstInputFile'     => ['required', 'string'],
            'secondInputFile'    => [
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
            'inputType'          => ['required', Rule::in(self::VALID_INPUT_TYPES)],
            'convertBam'         => ['filled', 'boolean'],
            'trimGalore'         => ['filled', 'array'],
            'trimGalore.enable'  => ['filled', 'boolean'],
            'trimGalore.quality' => ['filled', 'integer'],
            'trimGalore.length'  => ['filled', 'integer'],
            'transcriptome'      => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
            'threads'            => ['filled', 'integer'],
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

        return true;
    }

    /**
     * @param bool                  $paired
     * @param string                $firstInputFile
     * @param string|null           $secondInputFile
     * @param string                $inputType
     * @param \App\Models\Reference $transcriptome
     * @param int                   $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runSalmon(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        string $inputType,
        Reference $transcriptome,
        int $threads = 1
    ): array {
        if (!$transcriptome->isAvailableFor('salmon')) {
            throw new ProcessingJobException('The specified reference sequence is not indexed for salmon analysis.');
        }
        $this->log('Computing counts using Salmon.');
        $salmonOutputRelative = $this->model->getJobTempFile('salmon_output', '_sa.txt');
        $salmonOutput = $this->model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = Storage::disk('public')->url($salmonOutputRelative);
        switch ($inputType) {
            case self::FASTQ:
                $command = [
                    'bash',
                    self::scriptPath('salmon_counting.sh'),
                    '-i',
                    $transcriptome->basename(),
                    '-f',
                    $firstInputFile,
                    '-t',
                    $threads,
                    '-o',
                    $salmonOutput,
                ];
                if ($paired) {
                    $command[] = '-s';
                    $command[] = $secondInputFile;
                }
                $output = self::runCommand(
                    $command,
                    $this->model->getAbsoluteJobDirectory(),
                    null,
                    null,
                    [
                        3 => 'Input file does not exist.',
                        4 => 'Second input file does not exist.',
                        5 => 'Output file must be specified.',
                        6 => 'Output directory is not writable.',
                        7 => 'Indexed transcriptome does not exist.',
                        8 => 'Unable to find output file.',
                    ]
                );
                break;
            case self::BAM:
                $output = self::runCommand(
                    [
                        'bash',
                        self::scriptPath('salmon_counting_bam.sh'),
                        '-r',
                        $transcriptome->basename(),
                        '-i',
                        $firstInputFile,
                        '-t',
                        $threads,
                        '-o',
                        $salmonOutput,
                    ],
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
        $this->log($output);
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $this->log('Count computation completed.');

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
        $this->log('Starting analysis.');
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $convertBam = (bool)$this->model->getParameter('convertBam', false);
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->model->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->model->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 14);
        $transcriptomeName = $this->model->getParameter('transcriptome', env('HUMAN_TRANSCRIPTOME_NAME'));
        $transcriptome = Reference::whereName($transcriptomeName)->firstOrFail();
        $threads = (int)$this->model->getParameter('threads', 1);
        if ($inputType === self::BAM && $convertBam) {
            $this->log('Converting BAM to FASTQ.');
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq($this->model, $paired, $firstInputFile);
            $inputType = self::FASTQ;
            $this->log('BAM converted to FASTQ.');

        }
        [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
        if ($inputType === self::FASTQ && $trimGaloreEnable) {
            $this->log('Trimming reads using TrimGalore.');
            [$firstTrimmedFastq, $secondTrimmedFastq] = self::runTrimGalore(
                $this->model,
                $paired,
                $firstInputFile,
                $secondInputFile,
                $trimGaloreQuality,
                $trimGaloreLength
            );
            $this->log('Trimming completed.');
        }
        [$salmonOutput, $salmonOutputUrl] = $this->runSalmon(
            $paired,
            $firstTrimmedFastq,
            $secondTrimmedFastq,
            $inputType,
            $transcriptome,
            $threads
        );
        $this->model->setOutput(
            [
                'outputFile' => [
                    'path' => $salmonOutput,
                    'url'  => $salmonOutputUrl,
                ],
            ]
        );
        $this->log('Analysis completed.');
        $this->model->save();
    }


    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Runs mRNAs and/or lncRNAs (long RNAs reads) analysis from sequencing data';
    }
}
