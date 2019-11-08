<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class CircRnaJobType extends AbstractJob
{

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
            'customGTFFile'        => 'An optional GTF file for custom annotation of reads (Not needed for human hg19)',
            'customFASTAGenome'    => 'An optional Genome to employ for custom annotation (Not needed for human hg19)',
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
            'customGTFFile'        => ['filled', 'string'],
            'customFASTAGenome'    => ['filled', 'string'],
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
        $customGTFFile = $this->model->getParameter('customGTFFile');
        $customFASTAGenome = $this->model->getParameter('customFASTAGenome');
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
        if (!empty($customGTFFile) && !$disk->exists($dir . $customGTFFile)) {
            return false;
        }
        if (!empty($customFASTAGenome) && !$disk->exists($dir . $customFASTAGenome)) {
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
        $paired = (bool)$this->model->getParameter('paired', false);
        $inputType = $this->model->getParameter('inputType');
        $convertBam = (bool)$this->model->getParameter('convertBam', false);
        $firstInputFile = $this->model->getParameter('firstInputFile');
        $secondInputFile = $this->model->getParameter('secondInputFile');
        $trimGaloreEnable = (bool)$this->model->getParameter('trimGalore.enable', $inputType === self::FASTQ);
        $trimGaloreQuality = (int)$this->model->getParameter('trimGalore.quality', 20);
        $trimGaloreLength = (int)$this->model->getParameter('trimGalore.length', 14);
        $customGTFFile = $this->model->getParameter('customGTFFile');
        $customFASTAGenome = $this->model->getParameter('customFASTAGenome');
        $threads = (int)$this->model->getParameter('threads', 1);
        $ciriSpanningDistance = (int)$this->model->getParameter('ciriSpanningDistance', 500000);


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
