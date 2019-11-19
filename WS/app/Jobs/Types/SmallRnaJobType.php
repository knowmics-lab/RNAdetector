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
use App\Utils;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class SmallRnaJobType extends AbstractJob
{
    use ConvertsBamToFastqTrait, RunTrimGaloreTrait;

    private const FASTQ                = 'fastq';
    private const BAM                  = 'BAM';
    private const HTSEQ_COUNTS         = 'htseq';
    private const FEATURECOUNTS_COUNTS = 'feature-counts';
    private const VALID_INPUT_TYPES    = [self::FASTQ, self::BAM];
    private const VALID_COUNTS_METHOS  = [self::HTSEQ_COUNTS, self::FEATURECOUNTS_COUNTS];

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'paired'            => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'    => 'Required, input file for the analysis',
            'secondInputFile'   => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'         => 'Required, type of the input file (fastq, bam)',
            'convertBam'        => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'        => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'countingAlgorithm' => 'The counting algorithm htseq or feature-counts (Default htseq)',
            'genome'            => 'An optional name for a reference genome (Default human hg19)',
            'annotation'        => 'An optional name for a genome annotation (Default human hg19)',
            'threads'           => 'Number of threads for this analysis (Default 1)',
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
            'countingAlgorithm'  => ['filled', Rule::in(self::VALID_COUNTS_METHOS)],
            'genome'             => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
            'annotation'         => ['filled', 'alpha_dash', Rule::exists('annotations', 'name')],
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
        $countingAlgorithm = $this->model->getParameter('countingAlgorithm', self::HTSEQ_COUNTS);
        if (!in_array($inputType, self::VALID_INPUT_TYPES, true)) {
            return false;
        }
        if (!in_array($countingAlgorithm, self::VALID_COUNTS_METHOS, true)) {
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
     * Runs TopHat
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
    private function runTophat(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        Reference $genome,
        Annotation $annotation,
        int $threads = 1
    ): string {
        $bamOutput = $this->model->getJobTempFileAbsolute('bowtie_output', '.bam');
        // Call TopHat
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create TopHat output file');
        }

        return $bamOutput;
    }

    /**
     * Runs HTseq-count
     *
     * @param string                 $htseqInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runHTSEQ(string $htseqInputFile, Annotation $annotation, int $threads = 1): array
    {
        $htseqOutputRelative = $this->model->getJobTempFile('htseq_output', '.txt');
        $htseqOutput = $this->model->absoluteJobPath($htseqOutputRelative);
        $htseqOutputUrl = \Storage::disk('public')->url($htseqOutput);
        // Runs HTseq-count
        if (!file_exists($htseqOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }

        return [$htseqOutputRelative, $htseqOutputUrl];
    }

    /**
     * Runs FeatureCount
     *
     * @param string                 $featurecountInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runFeatureCount(string $featurecountInputFile, Annotation $annotation, int $threads = 1): array
    {
        $featurecountOutputRelative = $this->model->getJobTempFile('featurecount_output', '.txt');
        $featurecountOutput = $this->model->absoluteJobPath($featurecountOutputRelative);
        $featurecountOutputUrl = \Storage::disk('public')->url($featurecountOutput);
        // Runs FeatureCount
        if (!file_exists($featurecountOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }

        return [$featurecountOutputRelative, $featurecountOutputUrl];
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

    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'A greeting to the user';
    }
}
