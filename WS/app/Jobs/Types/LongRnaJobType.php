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
            'paired'               => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'       => 'Required, input file for the analysis',
            'secondInputFile'      => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputBAM'             => 'Optional input BAM file if FASTQ files are not provided',
            'inputType'            => 'Required, type of the input file (fastq, bam)',
            'convertBam'           => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'           => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
            'customFASTATranscriptome'    => 'An optional Transcriptome to employ for custom annotation (Not needed for human hg19 mRNAs and lncRNAs)',
            'customTranscriptomeName'     => 'An optional name for the custom transcriptome',
            'threads'              => 'Number of threads for this analysis (Default 1)',
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
            'firstInputFile'       => ['required', 'string'], # quà c'è da aggiungere la condizione che è richiesto solo se non viene fornito in alternativa al FASTQ un file BAM
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
            'inputBAM'                  => [], #quà c'è da aggiungere la condizione che il file BAM è richiesto qualora l'utente non fornisca alcun file FASTQ
            'inputType'                 => ['required', Rule::in(self::VALID_INPUT_TYPES)],
            'convertBam'                => ['filled', 'boolean'],
            'trimGalore'                => ['filled', 'array'],
            'trimGalore.enable'         => ['filled', 'boolean'],
            'trimGalore.quality'        => ['filled', 'integer'],
            'trimGalore.length'         => ['filled', 'integer'],
            'customFASTATranscriptome'  => ['filled', 'string'],
            'customTranscriptomeName'   => ['filled', 'alpha_num'],
            'threads'                   => ['filled', 'integer'],
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
        $inputBAM = $this->model->getParameter('inputBAM');
        $customFASTATranscriptome= $this->model->getParameter('customFASTATranscriptome');
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
     * Run Salmon analysis
     *
     */

    private function runSalmon(
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        int $threads = 1,
        ?string $customTranscriptomeName = null
    ): string {
        if ($customTranscriptomeName === null) {
            $transcriptomeDir = env('HUMAN_SALMON_INDEXED_TRANSCRIPTOME');
            $index = false;
        } else {
            $transcriptomeDir = env('CUSTOM_TRANSCRIPTOME_PATH') . '/' . $customTranscriptomeName;
            if (!file_exists($transcriptomeDir) || !is_dir($transcriptomeDir)) {
                $index = true;
            }
        }
        if ($index) {
            // Call Salmon index script
            if (!file_exists($transcriptomeDir) && !is_dir($transcriptomeDir)) {
                throw new ProcessingJobException('Unable to create indexed transcriptome');
            }
        }
        $salmonOutput = $this->model->getJobTempFileAbsolute('salmon_output', '.txt');
        // Call Salmon counting scripts
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        return $salmonOutput;
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
        try {
            $name = $this->model->getParameter('name', \Auth::user()->name);
            $this->model->setOutput('greetings', 'Hello ' . $name . '!!');
        } catch (\Exception $e) {
            throw new ProcessingJobException('An error occurred during job processing.', 0, $e);
        }
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
