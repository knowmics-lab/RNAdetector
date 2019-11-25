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
use App\Models\Annotation;
use App\Models\Reference;
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
            'genome'               => 'An optional name for a reference genome (Default human hg19)',
            'annotation'           => 'An optional name for a genome annotation (Default human hg19)',
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
            'genome'               => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
            'annotation'           => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
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
     * Runs BWA analysis
     *
     * @param bool                   $paired
     * @param string                 $firstInputFile
     * @param string|null            $secondInputFile
     * @param \App\Models\Reference  $genome
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runBWA(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        Reference $genome,
        Annotation $annotation,
        int $threads = 1
    ): string {
        if (!$genome->isAvailableFor('bwa')) {
            throw new ProcessingJobException('The specified genome is not indexed for BWA analysis.');
        }
        $samOutput = $this->model->getJobTempFileAbsolute('bwa_output', '.sam');
        $command = [
            'bash',
            self::scriptPath('bwa.bash'),
            '-a',
            $annotation->path,
            '-g',
            $genome->basename(),
            '-t',
            $threads,
            '-o',
            $samOutput,
            '-f',
            $firstInputFile,
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
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Second input file does not exist.',
                6 => 'Output file must be specified.',
                7 => 'Output directory is not writable.',
                8 => 'Unable to find bwa output file.',
            ]
        );
        if (!file_exists($samOutput)) {
            throw new ProcessingJobException('Unable to create BWA output file');
        }
        $this->log($output);

        return $samOutput;
    }

    /**
     * Runs CIRI analysis
     *
     * @param bool                   $paired
     * @param string                 $ciriInputFile
     * @param \App\Models\Reference  $genome
     * @param \App\Models\Annotation $annotation
     * @param int                    $spanningDistance
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runCIRI(
        bool $paired,
        string $ciriInputFile,
        Reference $genome,
        Annotation $annotation,
        int $spanningDistance = 500000
    ): array {
        $ciriOutputRelative = $this->model->getJobTempFile('ciri_output', '_ci.txt');
        $ciriOutput = $this->model->absoluteJobPath($ciriOutputRelative);
        $ciriOutputUrl = \Storage::disk('public')->url($ciriOutputRelative);
        $output = AbstractJob::runCommand(
            [
                'bash',
                self::scriptPath('ciri.bash'),
                '-a',
                $annotation->path,
                '-f',
                $genome->path,
                '-s',
                ($paired) ? 'paired' : 'single',
                '-m',
                $spanningDistance,
                '-i',
                $ciriInputFile,
                '-o',
                $ciriOutput,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Sequencing strategy not valid.',
                6 => 'Output file must be specified.',
                7 => 'Output directory is not writable.',
                8 => 'FASTA file does not exist.',
                9 => 'Unable to find CIRI output file.',
            ]
        );
        if (!file_exists($ciriOutput)) {
            throw new ProcessingJobException('Unable to create CIRI output file');
        }
        $this->log($output);

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
        $this->log('Starting CircRNA analysis.');
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $convertBam = (bool)$this->model->getParameter('convertBam', false);
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->model->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->model->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 14);
        $genomeName = $this->model->getParameter('genome', env('HUMAN_GENOME_NAME'));
        $annotationName = $this->model->getParameter('annotation', env('HUMAN_CIRI_ANNOTATION_NAME'));
        $threads = (int)$this->model->getParameter('threads', 1);
        $ciriSpanningDistance = (int)$this->model->getParameter('ciriSpanningDistance', 500000);
        $genome = Reference::whereName($genomeName)->firstOrFail();
        $annotation = Annotation::whereName($annotationName)->firstOrFail();
        $ciriInputFile = null;
        if ($inputType === self::BAM && $convertBam) {
            $inputType = self::FASTQ;
            $this->log('Converting BAM to FASTQ.');
            [$firstInputFile, $secondInputFile, $bashOutput] = self::convertBamToFastq(
                $this->model,
                $paired,
                $firstInputFile
            );
            $this->log($bashOutput);
            $this->log('BAM converted to FASTQ.');
        }
        if ($inputType === self::FASTQ) {
            [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
            if ($trimGaloreEnable) {
                $this->log('Trimming reads using TrimGalore');
                [$firstTrimmedFastq, $secondTrimmedFastq, $bashOutput] = self::runTrimGalore(
                    $this->model,
                    $paired,
                    $firstInputFile,
                    $secondInputFile,
                    $trimGaloreQuality,
                    $trimGaloreLength,
                    true
                );
                $this->log($bashOutput);
                $this->log('Trimming completed');
            }
            $this->log('Aligning reads with BWA');
            $ciriInputFile = $this->runBWA(
                $paired,
                $firstTrimmedFastq,
                $secondTrimmedFastq,
                $genome,
                $annotation,
                $threads
            );
            $this->log('Alignment completed.');
        } elseif ($inputType === self::BAM) {
            $ciriInputFile = $this->model->getJobTempFileAbsolute('bam2sam', '.sam');
            $this->log('Converting BAM to SAM.');
            $output = self::runCommand(
                [
                    'bash',
                    self::scriptPath('bam2sam.sh'),
                    '-b',
                    $firstInputFile,
                    '-o',
                    $ciriInputFile,
                ],
                $this->model->getAbsoluteJobDirectory(),
                null,
                null,
                [
                    3 => 'Input file does not exist.',
                    4 => 'Output file must be specified.',
                    5 => 'Output directory is not writable.',
                ]
            );
            $this->log($output);
            $this->log('BAM converted to SAM.');
            if (!file_exists($ciriInputFile)) {
                throw new ProcessingJobException('Unable to create converted BAM file');
            }
        } else {
            $ciriInputFile = $firstInputFile;
        }
        $this->log('Computing counts of CircRNA using CIRI.');
        [$ciriOutput, $ciriOutputUrl] = $this->runCIRI(
            $paired,
            $ciriInputFile,
            $genome,
            $annotation,
            $ciriSpanningDistance
        );
        $this->log('CircRNA Analysis completed!');
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
