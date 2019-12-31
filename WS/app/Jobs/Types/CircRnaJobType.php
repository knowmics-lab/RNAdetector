<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsBamToFastqTrait;
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Models\Annotation;
use App\Models\Reference;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class CircRnaJobType extends AbstractJob
{

    use HasCommonParameters, ConvertsBamToFastqTrait, RunTrimGaloreTrait;

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return array_merge(
            self::commonParametersSpec(),
            [
                'trimGalore'           => [
                    'enable'   => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                    'quality'  => 'Minimal PHREAD quality for trimming (Default 20)',
                    'length'   => 'Minimal reads length (Default 40)',
                    'hardTrim' => 'A boolean value to indicate if reads should be trimmed to the same size (Default true)',
                ],
                'genome'               => 'An optional name for a reference genome (Default human hg19)',
                'annotation'           => 'An optional name for a GTF genome annotation (Default human hg19)',
                'threads'              => 'Number of threads for this analysis (Default 1)',
                'useFastqPair'         => 'A boolean value to indicate whether to use fastq_pair or bbmap repair (Default false=bbmap repair)',
                'ciriSpanningDistance' => 'The maximum spanning distance used in CIRI (Default 200000)',
                'ciriUseVersion1'      => 'A boolean value to indicate whether to use CIRI 1 or CIRI 2 (Default false)',
            ]
        );
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
        return array_merge(
            self::commonParametersValidation($request),
            [
                'trimGalore.length'    => ['filled', 'integer', 'min:40'],
                'trimGalore.hardTrim'  => ['filled', 'boolean'],
                'genome'               => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'annotation'           => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
                'threads'              => ['filled', 'integer'],
                'ciriSpanningDistance' => ['filled', 'integer'],
            ]
        );
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        return $this->validateCommonParameters($this->model, self::VALID_INPUT_TYPES, self::FASTQ);
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
     * @param bool                   $useFastqPair
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
        int $threads = 1,
        bool $useFastqPair = false
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
            if ($useFastqPair) {
                $command[] = '-p';
            }
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
     * @param string                 $ciriInputFile
     * @param \App\Models\Reference  $genome
     * @param \App\Models\Annotation $annotation
     * @param int                    $spanningDistance
     * @param int                    $threads
     * @param bool                   $paired
     * @param bool                   $useCiri1
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runCIRI(
        string $ciriInputFile,
        Reference $genome,
        Annotation $annotation,
        int $spanningDistance = 200000,
        int $threads = 1,
        bool $paired = false,
        bool $useCiri1 = false
    ): array {
        $ciriOutputRelative = $this->model->getJobTempFile('ciri_output', '_ci.txt');
        $ciriOutput = $this->model->absoluteJobPath($ciriOutputRelative);
        $ciriOutputUrl = \Storage::disk('public')->url($ciriOutputRelative);
        $command = [
            'bash',
            self::scriptPath('ciri.bash'),
            '-a',
            $annotation->path,
            '-f',
            $genome->path,
            '-t',
            $threads,
            '-m',
            $spanningDistance,
            '-i',
            $ciriInputFile,
            '-o',
            $ciriOutput,
        ];
        if ($paired) {
            $command[] = '-p';
        }
        if ($useCiri1) {
            $command[] = '-1';
        }
        $output = AbstractJob::runCommand(
            $command,
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'CIRI returned non-zero exit code.',
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
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 40);
        $trimGaloreHardTrim = (bool)$this->model->getParameter('trimGalore.hardTrim', true);
        $genomeName = $this->model->getParameter('genome', env('HUMAN_GENOME_NAME'));
        $annotationName = $this->model->getParameter('annotation', env('HUMAN_CIRI_ANNOTATION_NAME'));
        $threads = (int)$this->model->getParameter('threads', 1);
        $ciriSpanningDistance = (int)$this->model->getParameter('ciriSpanningDistance', 200000);
        $useFastqPair = (bool)$this->model->getParameter('useFastqPair', false);
        $useCiri1 = (bool)$this->model->getParameter('ciriUseVersion1', false);
        $genome = Reference::whereName($genomeName)->firstOrFail();
        $annotation = Annotation::whereName($annotationName)->firstOrFail();
        if ($annotation->type !== 'gtf') {
            throw new ProcessingJobException('You must select only GTF annotations');
        }
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
                    $trimGaloreHardTrim,
                    $threads
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
                $threads,
                $useFastqPair
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
            $ciriInputFile,
            $genome,
            $annotation,
            $ciriSpanningDistance,
            $threads,
            $paired,
            $useCiri1
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
