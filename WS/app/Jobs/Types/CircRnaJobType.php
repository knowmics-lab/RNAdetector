<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsBamToFastqTrait;
use App\Jobs\Types\Traits\ConvertsSamToBamTrait;
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Models\Annotation;
use App\Models\Reference;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class CircRnaJobType extends AbstractJob
{

    use HasCommonParameters, ConvertsBamToFastqTrait, RunTrimGaloreTrait, ConvertsSamToBamTrait;

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
                'ciriQuant'            => 'An optional boolean to indicate whether to use ciri-quant instead of ciri (Default false)',
                'genome'               => 'An optional name for a reference genome (Default human hg19)',
                'annotation'           => 'An optional name for a GTF genome annotation (Default human hg19)',
                'bedAnnotation'        => 'An optional annotation in BED format for circRNAs junctions. Required for ciriQuant analysis.',
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
            'outputFile' => 'Formatted read counts files',
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

        $parameters = (array)$request->get('parameters', []);
        $requiredIfQuant = static function () use ($parameters) {
            return (bool)data_get($parameters, 'ciriQuant', false);
        };

        return array_merge(
            self::commonParametersValidation($request),
            [
                'trimGalore.length'    => ['filled', 'integer', 'min:40'],
                'trimGalore.hardTrim'  => ['filled', 'boolean'],
                'genome'               => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'annotation'           => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
                'threads'              => ['filled', 'integer'],
                'ciriSpanningDistance' => ['filled', 'integer'],
                'ciriQuant'            => ['filled', 'boolean'],
                'bedAnnotation'        => [Rule::requiredIf($requiredIfQuant), Rule::exists('annotations', 'name')],
                'paired'               => [
                    'filled',
                    'boolean',
                    static function ($attribute, $value, $fail) use ($parameters) {
                        $quant = (bool)data_get($parameters, 'ciriQuant', false);
                        if ($quant && !$value) {
                            $fail('Only paired-end sequencing is supported for CIRIquant analysis');
                        }
                    },
                ],
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
        if (!$this->validateCommonParameters($this->model, self::VALID_INPUT_TYPES, self::FASTQ)) {
            return false;
        }
        $quant = (bool)$this->getParameter('ciriQuant', false);
        if ($quant) {
            $paired = (bool)$this->getParameter('paired', false);
            if (!$paired) {
                return false;
            }
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
        $samOutput = $this->model->getJobFileAbsolute('bwa_output_', '.sam');
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
            function ($type, $buffer) {
                $this->log(trim($buffer));
            },
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

        // $this->log($output);

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
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $ciriOutputRelative = $this->model->getJobFile('ciri_output_', '_ci.txt');
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
            function ($type, $buffer) {
                $this->log(trim($buffer));
            },
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

        // $this->log($output);

        return [$ciriOutputRelative, $ciriOutputUrl];
    }

    /**
     * Make CIRIquant configuration file and returns its path
     *
     * @param \App\Models\Reference  $reference
     * @param \App\Models\Annotation $annotation
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function makeQuantConfig(Reference $reference, Annotation $annotation): string
    {
        if (!$reference->isAvailableFor('bwa')) {
            throw new ProcessingJobException('The selected reference has not been indexed for BWA.');
        }
        if (!$reference->isAvailableFor('hisat')) {
            throw new ProcessingJobException('The selected reference has not been indexed for HISAT.');
        }
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $configFile = $this->getJobFile('quant_config_', '.yml');
        $name = basename($configFile, '.yml');
        $template = file_get_contents(resource_path('templates/quant_config.yml'));
        $template = str_replace(
            [
                '{NAME}',
                '{GENOME_FASTA}',
                '{GENOME_GTF}',
                '{GENOME_INDEX}',
            ],
            [
                $name,
                $reference->path,
                $annotation->path,
                $reference->basename(),
            ],
            $template
        );
        @file_put_contents($configFile, $template);
        @chmod($configFile, 0777);

        return $configFile;
    }

    /**
     * Runs CIRI analysis
     *
     * @param string                 $firstInputFile
     * @param string                 $secondInputFile
     * @param string                 $configFile
     * @param \App\Models\Annotation $bedAnnotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runCIRIQuant(
        string $firstInputFile,
        string $secondInputFile,
        string $configFile,
        Annotation $bedAnnotation,
        int $threads = 1
    ): array {
        if (!$bedAnnotation->isBed()) {
            throw new ProcessingJobException('The BED annotation file is not valid.');
        }
        $quantOutputRelative = $this->model->getJobFile('quant_output_', '.gtf');
        $quantOutput = $this->model->absoluteJobPath($quantOutputRelative);
        $quantOutputUrl = \Storage::disk('public')->url($quantOutputRelative);
        $command = [
            'bash',
            self::scriptPath('ciri_quant.sh'),
            '-f',
            $firstInputFile,
            '-s',
            $secondInputFile,
            '-c',
            $configFile,
            '-b',
            $bedAnnotation->path,
            '-t',
            $threads,
            '-o',
            $quantOutput,
        ];
        AbstractJob::runCommand(
            $command,
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log(trim($buffer));
            },
            [
                3  => 'First input file does not exist.',
                4  => 'Second input file does not exist.',
                5  => 'Configuration file does not exist.',
                6  => 'BED annotation file does not exist.',
                7  => 'Output directory is not specified.',
                8  => 'Output directory is not writeable.',
                9  => 'Unknown error during CIRIquant execution.',
                10 => 'Unable to find CIRIquant output file.',
            ]
        );
        if (!file_exists($quantOutput)) {
            throw new ProcessingJobException('Unable to create CIRIquant output file');
        }

        return [$quantOutputRelative, $quantOutputUrl];
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
        $paired = (bool)$this->getParameter('paired', false);
        $inputType = $this->getParameter('inputType');
        $convertBam = (bool)$this->getParameter('convertBam', false);
        $firstInputFile = $this->getParameter('firstInputFile');
        $secondInputFile = $this->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->getParameter('trimGalore.length', 40);
        $trimGaloreHardTrim = (bool)$this->getParameter('trimGalore.hardTrim', true);
        $genomeName = $this->getParameter('genome', env('HUMAN_GENOME_NAME'));
        $annotationName = $this->getParameter('annotation', env('HUMAN_CIRC_ANNOTATION_NAME'));
        $threads = (int)$this->getParameter('threads', 1);
        $ciriSpanningDistance = (int)$this->getParameter('ciriSpanningDistance', 200000);
        $useFastqPair = (bool)$this->getParameter('useFastqPair', false);
        $useCiri1 = (bool)$this->getParameter('ciriUseVersion1', false);
        $ciriQuant = (bool)$this->getParameter('ciriQuant', false);
        $bedAnnotationName = $this->getParameter('bedAnnotation');
        $bedAnnotation = null;
        if ($ciriQuant) {
            $bedAnnotation = Annotation::whereName($bedAnnotationName)->firstOrFail();
        }
        $genome = Reference::whereName($genomeName)->firstOrFail();
        $annotation = Annotation::whereName($annotationName)->firstOrFail();
        if ($annotation->type !== 'gtf') {
            throw new ProcessingJobException('You must select only GTF annotations');
        }
        $ciriInputFile = null;
        if ($inputType === self::SAM && $ciriQuant) {
            $inputType = self::BAM;
            $firstInputFile = self::convertSamToBam($this->model, $firstInputFile);
        }
        if ($inputType === self::BAM && ($convertBam || $ciriQuant)) {
            $inputType = self::FASTQ;
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq(
                $this->model,
                $paired,
                $firstInputFile
            );
        }
        if ($ciriQuant) {
            if ($inputType !== self::BAM) {
                throw new ProcessingJobException('Only FASTQ files are supported for CIRIquant analysis.');
            }
            if ($bedAnnotation === null) {
                throw new ProcessingJobException('No valid BED file specified.');
            }
            $this->log('Building CIRIquant config file');
            $configFile = $this->makeQuantConfig($genome, $annotation);
            $this->log('Starting CIRIquant analysis');
            [$circOutput, $circOutputUrl] = $this->runCIRIQuant($firstInputFile, $secondInputFile, $configFile, $bedAnnotation, $threads);
        } else {
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
                $ciriInputFile = $this->model->getJobFileAbsolute('bam2sam_', '.sam');
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
                    function ($type, $buffer) {
                        $this->log(trim($buffer));
                    },
                    [
                        3 => 'Input file does not exist.',
                        4 => 'Output file must be specified.',
                        5 => 'Output directory is not writable.',
                    ]
                );
                // $this->log($output);
                $this->log('BAM converted to SAM.');
                if (!file_exists($ciriInputFile)) {
                    throw new ProcessingJobException('Unable to create converted BAM file');
                }
            } else {
                $ciriInputFile = $firstInputFile;
            }
            $this->log('Computing counts of CircRNA using CIRI.');
            [$circOutput, $circOutputUrl] = $this->runCIRI(
                $ciriInputFile,
                $genome,
                $annotation,
                $ciriSpanningDistance,
                $threads,
                $paired,
                $useCiri1
            );
        }
        $this->log('CircRNA Analysis completed!');
        $this->model->setOutput(
            [
                'outputFile' => [
                    'path' => $circOutput,
                    'url'  => $circOutputUrl,
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
