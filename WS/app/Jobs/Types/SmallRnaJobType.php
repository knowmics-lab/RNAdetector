<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsBamToFastqTrait;
use App\Jobs\Types\Traits\ConvertsSamToBamTrait;
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Models\Annotation;
use App\Models\Reference;
use App\Utils;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class SmallRnaJobType extends AbstractJob
{
    use HasCommonParameters, ConvertsBamToFastqTrait, ConvertsSamToBamTrait, RunTrimGaloreTrait;

    private const HTSEQ_COUNTS         = 'htseq';
    private const FEATURECOUNTS_COUNTS = 'feature-counts';
    private const SALMON               = 'salmon';
    private const VALID_COUNTS_METHODS = [self::HTSEQ_COUNTS, self::FEATURECOUNTS_COUNTS, self::SALMON];

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
                'countingAlgorithm' => 'The counting algorithm htseq, feature-counts, or salmon (Default htseq)',
                'genome'            => 'An optional name for a reference genome (Default human hg19)',
                'transcriptome'     => 'An optional name for a transcriptome if counting algorithm is salmon (Default human hg19)',
                'annotation'        => 'An optional name for a genome annotation (Default human hg19)',
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
                'countingAlgorithm' => ['filled', Rule::in(self::VALID_COUNTS_METHODS)],
                'genome'            => ['filled', 'alpha_dash', Rule::exists('references', 'name')],
                'transcriptome'     => [
                    Rule::requiredIf(
                        static function () use ($parameters) {
                            return ($parameters['countingAlgorithm'] ?? self::HTSEQ_COUNTS) === self::SALMON;
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
        $countingAlgorithm = $this->model->getParameter('countingAlgorithm', self::HTSEQ_COUNTS);
        if (!in_array($countingAlgorithm, self::VALID_COUNTS_METHODS, true)) {
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
        $bamOutput = $this->model->getJobTempFileAbsolute('tophat_output', '.bam');
        $command = [
            'bash',
            self::scriptPath('tophat.bash'),
            '-a',
            $annotation->path,
            '-g',
            $genome->basename(),
            '-t',
            $threads,
            '-f',
            $firstInputFile,
            '-o',
            $bamOutput,
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
                8 => 'Unable to find output bam file.',
            ]
        );
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create TopHat output file');
        }
        $this->log($output);

        return $bamOutput;
    }

    /**
     * Runs HTseq-count
     *
     * @param string                 $countingInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runHTSEQ(string $countingInputFile, Annotation $annotation, int $threads = 1): array
    {
        $htseqOutputRelative = $this->model->getJobTempFile('htseq_output', '_ht.txt');
        $htseqOutput = $this->model->absoluteJobPath($htseqOutputRelative);
        $htseqOutputUrl = \Storage::disk('public')->url($htseqOutputRelative);
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('htseqcount.bash'),
                '-a',
                $annotation->path,
                '-b',
                $countingInputFile,
                '-t',
                $threads,
                '-o',
                $htseqOutput,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
            ]
        );
        if (!file_exists($htseqOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }
        $this->log($output);

        return [$htseqOutputRelative, $htseqOutputUrl];
    }

    /**
     * Runs FeatureCount
     *
     * @param string                 $countingInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runFeatureCount(string $countingInputFile, Annotation $annotation, int $threads = 1): array
    {
        $featurecountOutputRelative = $this->model->getJobTempFile('featurecount_output', '_fc.txt');
        $featurecountOutput = $this->model->absoluteJobPath($featurecountOutputRelative);
        $featurecountOutputUrl = \Storage::disk('public')->url($featurecountOutputRelative);
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('htseqcount.bash'),
                '-a',
                $annotation->path,
                '-b',
                $countingInputFile,
                '-t',
                $threads,
                '-o',
                $featurecountOutput,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
            ]
        );

        if (!file_exists($featurecountOutput)) {
            throw new ProcessingJobException('Unable to create FeatureCount output file');
        }
        $this->log($output);

        return [$featurecountOutputRelative, $featurecountOutputUrl];
    }

    /**
     * Run salmon for SmallRNA counting
     *
     * @param bool                  $paired
     * @param string                $topHatInputFile
     * @param \App\Models\Reference $transcriptome
     * @param int                   $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runSalmonCount(bool $paired, string $topHatInputFile, Reference $transcriptome, int $threads = 1): array
    {
        if (!$transcriptome->isAvailableFor('salmon')) {
            throw new ProcessingJobException('The specified reference sequence is not indexed for salmon analysis.');
        }
        $this->log('Converting TopHat output to FastQ');
        [$firstTempFastQ, $secondTempFastQ, $output] = self::convertBamToFastq($this->model, $paired, $topHatInputFile);
        $this->log($output);
        $this->log('TopHat output converted to FastQ');
        $this->log('Running salmon');
        $salmonOutputRelative = $this->model->getJobTempFile('salmon_output', '_sa.txt');
        $salmonOutput = $this->model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = \Storage::disk('public')->url($salmonOutputRelative);
        $command = [
            'bash',
            self::scriptPath('salmon_counting.sh'),
            '-i',
            $transcriptome->basename(),
            '-t',
            $threads,
            '-o',
            $salmonOutput,
            '-f',
            $firstTempFastQ,
        ];
        if ($paired) {
            $command[] = '-s';
            $command[] = $secondTempFastQ;
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
                7 => 'Indexed trascriptome does not exist.',
                8 => 'Unable to find output file.',
                9 => 'An error occurred during salmon quant execution.',
            ]
        );
        $this->log($output);
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $this->log('Count computation completed.');
        if (file_exists($firstTempFastQ)) {
            @unlink($firstTempFastQ);
        }
        if (file_exists($secondTempFastQ)) {
            @unlink($secondTempFastQ);
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
        $this->log('Starting small ncRNAs analysis.');
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $convertBam = (bool)$this->model->getParameter('convertBam', false);
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->model->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->model->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 14);
        $genomeName = $this->model->getParameter('genome', env('HUMAN_GENOME_NAME'));
        $annotationName = $this->model->getParameter('annotation', env('HUMAN_SNCRNA_ANNOTATION_NAME'));
        $threads = (int)$this->model->getParameter('threads', 1);
        $countingAlgorithm = $this->model->getParameter('countingAlgorithm', self::HTSEQ_COUNTS);
        $genome = Reference::whereName($genomeName)->firstOrFail();
        $annotation = Annotation::whereName($annotationName)->firstOrFail();
        if ($inputType === self::SAM) {
            $inputType = self::BAM;
            $this->log('Converting SAM to BAM.');
            [$firstInputFile, $bashOutput] = self::convertSamToBam($this->model, $firstInputFile);
            $this->log($bashOutput);
            $this->log('SAM converted to BAM.');
        }
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
        [$firstTrimmedFastq, $secondTrimmedFastq] = [$firstInputFile, $secondInputFile];
        if ($inputType === self::FASTQ && $trimGaloreEnable) {
            $this->log('Trimming reads using TrimGalore.');
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
            $this->log('Trimming completed.');
        }
        $this->log('Aligning reads with TopHat');
        $countingInputFile = $this->runTophat(
            $paired,
            $firstTrimmedFastq,
            $secondTrimmedFastq,
            $genome,
            $annotation,
            $threads
        );
        $this->log('Alignment completed.');
        switch ($countingAlgorithm) {
            case self::HTSEQ_COUNTS:
                $this->log('Starting reads counting with HTseq-count');
                [$outputFile, $outputUrl] = $this->runHTSEQ($countingInputFile, $annotation, $threads);
                $this->log('Reads counting completed');
                break;
            case self::FEATURECOUNTS_COUNTS:
                $this->log('Starting reads counting with FeatureCount');
                [$outputFile, $outputUrl] = $this->runFeatureCount($countingInputFile, $annotation, $threads);
                $this->log('Reads counting completed');
                break;
            case self::SALMON:
                $transcriptomeName = $this->model->getParameter('transcriptome', env('HUMAN_TRANSCRIPTOME_SNCRNA_NAME'));
                $transcriptome = Reference::whereName($transcriptomeName)->firstOrFail();
                $this->log('Starting reads counting with Salmon');
                [$outputFile, $outputUrl] = $this->runSalmonCount($paired, $countingInputFile, $transcriptome, $threads);
                $this->log('Reads counting completed');
                break;
            default:
                throw new ProcessingJobException("Invalid counting algorithm");
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
