<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsSamToBamTrait;
use App\Jobs\Types\Traits\HandlesCompressedFastqTrait;
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\IndexesBAMTrait;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Jobs\Types\Traits\UseAlignmentTrait;
use App\Jobs\Types\Traits\UseCountingTrait;
use App\Jobs\Types\Traits\UseGenome;
use App\Jobs\Types\Traits\UseGenomeAnnotation;
use App\Jobs\Types\Traits\UsesJBrowseTrait;
use App\Jobs\Types\Traits\UseTranscriptome;
use App\Models\Annotation;
use App\Models\Job;
use App\Models\Reference;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class LongRnaJobType extends AbstractJob
{
    use HasCommonParameters, ConvertsSamToBamTrait, RunTrimGaloreTrait, UseAlignmentTrait, UseCountingTrait, HandlesCompressedFastqTrait;
    use UseTranscriptome, UseGenome, UseGenomeAnnotation, IndexesBAMTrait, UsesJBrowseTrait;

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
                'algorithm'         => 'Alignment/quantification algorithm: salmon, hisat2, star (Default: star)',
                'countingAlgorithm' => 'Counting algorithm in case of alignment: htseq, feature-counts, or salmon (Default feature-counts)',
                'genome'            => 'An optional genome for star or hisat2 (Default: human hg19)',
                'annotation'        => 'An optional annotation file for counting (Default: human hg19 mRNAs and lncRNAs)',
                'transcriptome'     => 'An optional transcriptome for quantification with salmon (Default: human hg19 mRNAs and lncRNAs)',
                'threads'           => 'Number of threads for this analysis (Default 1)',
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
            'outputFile'                => 'Raw output file',
            'harmonizedFile'            => 'Harmonized output file',
            'harmonizedTranscriptsFile' => 'Harmonized output file for transcripts expression',
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
        $requiredAlignment = static function () use ($parameters) {
            $algorithm = data_get($parameters, 'algorithm', self::SALMON);

            return $algorithm === self::HISAT2 || $algorithm === self::STAR;
        };

        return array_merge(
            self::commonParametersValidation($request),
            [
                'algorithm'         => ['filled', Rule::in(self::VALID_ALIGN_QUANT_METHODS)],
                'countingAlgorithm' => [Rule::requiredIf($requiredAlignment), Rule::in(self::VALID_COUNTING_METHODS)],
                'genome'            => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'annotation'        => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
                'transcriptome'     => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'threads'           => ['filled', 'integer'],
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
        $algorithm = $this->model->getParameter('algorithm', self::STAR);
        if (!in_array($algorithm, self::VALID_ALIGN_QUANT_METHODS, true)) {
            return false;
        }
        $countingAlgorithm = $this->model->getParameter('countingAlgorithm', self::FEATURECOUNTS_COUNTS);
        if (!in_array($countingAlgorithm, self::VALID_COUNTING_METHODS, true)) {
            return false;
        }

        return true;
    }

    /**
     * Checks if any valid genome and annotation has been provided for the analysis
     *
     * @param \App\Models\Reference|null  $genome
     * @param \App\Models\Annotation|null $annotation
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkGenomeAndAnnotation(?Reference $genome, ?Annotation $annotation): void
    {
        if ($genome === null) {
            throw new ProcessingJobException('An invalid genome has been provided');
        }
        if ($annotation === null) {
            throw new ProcessingJobException('An invalid genome annotation has been provided');
        }
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
        $threads = (int)$this->getParameter('threads', 1);
        $algorithm = $this->getParameter('algorithm', self::STAR);
        $countingAlgorithm = $this->getParameter('countingAlgorithm', self::FEATURECOUNTS_COUNTS);
        if ($inputType === self::SAM) {
            $firstInputFile = self::convertSamToBam($this->model, $firstInputFile);
            $inputType = self::BAM;
        }
        if ($inputType === self::BAM && $convertBam) {
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq($this->model, $paired, $firstInputFile);
            $inputType = self::FASTQ;
        }
        $outputFile = '';
        $outputUrl = '';
        $harmonizedGeneFile = '';
        $harmonizedGeneUrl = '';
        $harmonizedTxFile = null;
        $harmonizedTxUrl = null;
        $countingInputFile = '';
        $bamOutput = null;
        $count = true;
        if ($inputType === self::FASTQ) {
            $this->log('Checking if input is compressed...');
            $firstInputFile = self::checksForCompression($this->model, $firstInputFile);
            $secondInputFile = self::checksForCompression($this->model, $secondInputFile);
            [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
            if ($trimGaloreEnable) {
                [$firstTrimmedFastq, $secondTrimmedFastq] = self::runTrimGalore(
                    $this->model,
                    $paired,
                    $firstInputFile,
                    $secondInputFile,
                    $trimGaloreQuality,
                    $trimGaloreLength,
                    false,
                    $threads
                );
            }
            switch ($algorithm) {
                case self::HISAT2:
                    $countingInputFile = $this->runHisat(
                        $this->model,
                        $paired,
                        $firstTrimmedFastq,
                        $secondTrimmedFastq,
                        $this->getGenome(),
                        $threads
                    );
                    break;
                case self::STAR:
                    $countingInputFile = $this->runSTAR(
                        $this->model,
                        $paired,
                        $firstTrimmedFastq,
                        $secondTrimmedFastq,
                        $this->getGenome(),
                        $this->getGenomeAnnotation(),
                        $threads
                    );
                    break;
                case self::SALMON:
                    [
                        $outputFile,
                        $outputUrl,
                        $harmonizedGeneFile,
                        $harmonizedGeneUrl,
                        $harmonizedTxFile,
                        $harmonizedTxUrl,
                    ] = $this->runSalmon(
                        $this->model,
                        $paired,
                        $firstTrimmedFastq,
                        $secondTrimmedFastq,
                        $inputType,
                        $this->getTranscriptome(),
                        $threads
                    );
                    $count = false;
                    break;
            }
            if ($countingInputFile) {
                $bamOutput = $this->getFilePaths($countingInputFile);
            }
        } else {
            $countingInputFile = $firstInputFile;
            $bamOutput = $this->indexBAM($this->model, $countingInputFile, true, $threads);
        }
        if ($count && $countingInputFile) {
            switch ($countingAlgorithm) {
                case self::HTSEQ_COUNTS:
                    [$outputFile, $outputUrl, $harmonizedGeneFile, $harmonizedGeneUrl] = $this->runHTSEQ(
                        $this->model,
                        $countingInputFile,
                        $this->getGenomeAnnotation(),
                        $threads
                    );
                    break;
                case self::FEATURECOUNTS_COUNTS:
                    [$outputFile, $outputUrl, $harmonizedGeneFile, $harmonizedGeneUrl] = $this->runFeatureCount(
                        $this->model,
                        $countingInputFile,
                        $this->getGenomeAnnotation(),
                        $threads
                    );
                    break;
                case self::SALMON:
                    [
                        $outputFile,
                        $outputUrl,
                        $harmonizedGeneFile,
                        $harmonizedGeneUrl,
                        $harmonizedTxFile,
                        $harmonizedTxUrl,
                    ] = $this->runSalmonCount(
                        $this->model,
                        $paired,
                        $countingInputFile,
                        $this->getTranscriptome(),
                        $threads
                    );
                    break;
                default:
                    throw new ProcessingJobException('Invalid counting algorithm');
            }
        }
        $jbrowseConfig = $this->makeJBrowseConfig($this->model, $bamOutput, $this->getGenome(), $this->getGenomeAnnotation());
        $output = [
            'type'              => self::OUT_TYPE_ANALYSIS_HARMONIZED,
            'outputFile'        => [
                'path' => $outputFile,
                'url'  => $outputUrl,
            ],
            'outputBamFile'     => $this->getFilePathsForOutput($bamOutput),
            'outputJBrowseFile' => $this->getFilePathsForOutput($jbrowseConfig),
            'harmonizedFile'    => [
                'path' => $harmonizedGeneFile,
                'url'  => $harmonizedGeneUrl,
            ],
        ];
        if ($harmonizedTxFile !== null) {
            $output['type'] = self::OUT_TYPE_ANALYSIS_HARMONIZED_TRANSCRIPTS;
            $output['harmonizedTranscriptsFile'] = [
                'path' => $harmonizedTxFile,
                'url'  => $harmonizedTxUrl,
            ];
        }
        $this->model->setOutput($output);
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

    /**
     * @inheritDoc
     */
    public static function sampleGroupFunctions(): ?array
    {
        return [
            static function (Job $job) {
                return $job->sample_code;
            },
            static function (Job $job) {
                $output = $job->getOutput('harmonizedFile');

                return $job->absoluteJobPath($output['path']);
            },
            static function (Job $job) {
                $output = $job->getOutput('outputFile');

                return $job->absoluteJobPath($output['path']);
            },
            static function (Job $job) {
                $output = $job->getOutput('harmonizedTranscriptsFile');
                if ($output === null) {
                    return null;
                }

                return $job->absoluteJobPath($output['path']);
            },
            static function (Job $job) {
                $algorithm = $job->getParameter('algorithm', self::SALMON);
                $countingAlgorithm = $job->getParameter('countingAlgorithm', self::FEATURECOUNTS_COUNTS);
                $transcriptome = $job->getParameter('transcriptome', env('HUMAN_TRANSCRIPTOME_NAME'));
                $annotation = $job->getParameter('annotation', env('HUMAN_RNA_ANNOTATION_NAME'));

                return ($algorithm === self::SALMON || $countingAlgorithm === self::SALMON) ?
                    Reference::whereName($transcriptome)->firstOrFail()->map_path :
                    Annotation::whereName($annotation)->firstOrFail()->map_path;
            },
        ];
    }

    /**
     * @inheritDoc
     */
    public static function displayName(): string
    {
        return 'RNA-seq Analysis';
    }
}
