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
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Jobs\Types\Traits\UseAlignmentTrait;
use App\Jobs\Types\Traits\UseCountingTrait;
use App\Models\Annotation;
use App\Models\Reference;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class SmallRnaJobType extends AbstractJob
{
    use HasCommonParameters, ConvertsSamToBamTrait, RunTrimGaloreTrait, UseAlignmentTrait, UseCountingTrait;

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
                'algorithm'         => 'Alignment/quantification algorithm: salmon, tophat, hisat2 (Default: hisat2)',
                'countingAlgorithm' => 'The counting algorithm htseq, feature-counts, or salmon (Default: feature-counts)',
                'genome'            => 'An optional name for a reference genome (Default human hg19)',
                'transcriptome'     => 'An optional name for a transcriptome if counting algorithm is salmon (Default human hg19)',
                'annotation'        => 'An optional name for a GTF genome annotation (Default human hg19)',
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
        $parameters = (array)$request->get('parameters', []);

        return array_merge(
            self::commonParametersValidation($request),
            [
                'algorithm'         => ['filled', Rule::in(self::VALID_ALIGN_QUANT_METHODS)],
                'countingAlgorithm' => ['filled', Rule::in(self::VALID_COUNTING_METHODS)],
                'genome'            => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'transcriptome'     => [
                    Rule::requiredIf(
                        static function () use ($parameters) {
                            $countingSalmon = ($parameters['countingAlgorithm'] ?? self::HTSEQ_COUNTS) === self::SALMON;
                            $quantifSalmon = ($parameters['algorithm'] ?? self::HISAT2) === self::SALMON;

                            return $quantifSalmon || $countingSalmon;
                        }
                    ),
                    'alpha_dash',
                    Rule::exists('references', 'name'),
                ],
                'annotation'        => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
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
        $algorithm = $this->model->getParameter('algorithm', self::HISAT2);
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
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    public function handle(): void
    {
        $this->log('Starting small ncRNAs analysis.');
        $paired = (bool)$this->getParameter('paired', false);
        $inputType = $this->getParameter('inputType');
        $convertBam = (bool)$this->getParameter('convertBam', false);
        $firstInputFile = $this->getParameter('firstInputFile');
        $secondInputFile = $this->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->getParameter('trimGalore.length', 14);
        $genomeName = $this->getParameter('genome', env('HUMAN_GENOME_NAME'));
        $annotationName = $this->getParameter('annotation', env('HUMAN_SNCRNA_ANNOTATION_NAME'));
        $threads = (int)$this->getParameter('threads', 1);
        $algorithm = $this->getParameter('algorithm', self::HISAT2);
        $countingAlgorithm = $this->getParameter('countingAlgorithm', self::FEATURECOUNTS_COUNTS);
        $transcriptomeName = $this->getParameter('transcriptome', env('HUMAN_TRANSCRIPTOME_SNCRNA_NAME'));
        $transcriptome = Reference::whereName($transcriptomeName)->firstOrFail();
        $genome = Reference::whereName($genomeName)->firstOrFail();
        $annotation = Annotation::whereName($annotationName)->firstOrFail();
        if ($annotation->type !== 'gtf') {
            throw new ProcessingJobException('You must select only GTF annotations');
        }
        if ($inputType === self::SAM) {
            $firstInputFile = self::convertSamToBam($this->model, $firstInputFile);
            $inputType = self::BAM;
        }
        if ($inputType === self::BAM && $convertBam) {
            [$firstInputFile, $secondInputFile] = self::convertBamToFastq(
                $this->model,
                $paired,
                $firstInputFile
            );
            $inputType = self::FASTQ;
        }
        $outputFile = '';
        $outputUrl = '';
        $countingInputFile = '';
        $count = true;
        if ($inputType === self::FASTQ) {
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
                case self::TOPHAT:
                    $countingInputFile = $this->runTophat(
                        $this->model,
                        $paired,
                        $firstTrimmedFastq,
                        $secondTrimmedFastq,
                        $genome,
                        $annotation,
                        $threads
                    );
                    break;
                case self::HISAT2:
                    $countingInputFile = $this->runHisat($this->model, $paired, $firstTrimmedFastq, $secondTrimmedFastq, $genome, $threads);
                    break;
                case self::SALMON:
                    [$outputFile, $outputUrl] = $this->runSalmon(
                        $this->model,
                        $paired,
                        $firstTrimmedFastq,
                        $secondTrimmedFastq,
                        $inputType,
                        $transcriptome,
                        $threads
                    );
                    $count = false;
                    break;
            }
        } else {
            $countingInputFile = $firstInputFile;
        }
        if ($count) {
            switch ($countingAlgorithm) {
                case self::HTSEQ_COUNTS:
                    [$outputFile, $outputUrl] = $this->runHTSEQ($this->model, $countingInputFile, $annotation, $threads);
                    break;
                case self::FEATURECOUNTS_COUNTS:
                    [$outputFile, $outputUrl] = $this->runFeatureCount($this->model, $countingInputFile, $annotation, $threads);
                    break;
                case self::SALMON:
                    [$outputFile, $outputUrl] = $this->runSalmonCount($this->model, $paired, $countingInputFile, $transcriptome, $threads);
                    break;
                default:
                    throw new ProcessingJobException('Invalid counting algorithm');
            }
        }
        $this->model->setOutput(
            [
                'outputFile' => ['path' => $outputFile, 'url' => $outputUrl],
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
        return 'Runs small ncRNAs analysis from sequencing data';
    }
}
