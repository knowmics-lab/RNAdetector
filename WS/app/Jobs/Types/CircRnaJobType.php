<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use Illuminate\Validation\Rule;

class CircRnaJobType extends AbstractJob
{

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'paired' => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile' => 'Required, input file for the analysis',
            'secondInputFile' => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType' => 'Type of the input file (fastq, bam, sam)',
            'convertBam' => 'If inputType is bam converts input in another format: fastq or sam.',
            'trimGalore' => [
                'enable' => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length' => 'Minimal reads length (Default 14)',
            ],
            'customGTFFile' => 'An optional GTF file for custom annotation of reads (Not needed for human hg19)',
            'customFASTAGenome' => 'An optional Genome to employ for custom annotation (Not needed for human hg19)',
            'threads' => 'Number of threads for this analysis (Default 1)',
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
     * @return array
     */
    public static function validationSpec(): array
    {
        return [
            'paired' => ['filled', 'boolean'],
            'firstInputFile' => ['required', 'string'],
            'secondInputFile' => [Rule::requiredIf(function () {

            })],
            'inputType' => 'Type of the input file (fastq, bam, sam)',
            'convertBam' => 'If inputType is bam converts input in another format: fastq or sam.',
            'trimGalore' => [
                'enable' => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length' => 'Minimal reads length (Default 14)',
            ],
            'customGTFFile' => 'An optional GTF file for custom annotation of reads (Not needed for human hg19)',
            'customFASTAGenome' => 'An optional Genome to employ for custom annotation (Not needed for human hg19)',
            'threads' => 'Number of threads for this analysis (Default 1)',
            'ciriSpanningDistance' => 'The maximum spanning distance used in CIRI (Default 500000)',
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
