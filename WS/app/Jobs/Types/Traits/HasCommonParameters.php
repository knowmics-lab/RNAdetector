<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Models\Job;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

trait HasCommonParameters
{

    abstract protected function validateFileParameter(string $parameter): bool;

    /**
     * Get a list of common parameters specs
     *
     * @return array
     */
    private static function commonParametersSpec(): array
    {
        return [
            'paired'          => 'A boolean value to indicate whether sequencing strategy is paired-ended or not (Default false)',
            'firstInputFile'  => 'Required, input file for the analysis. FASTQ or BAM',
            'secondInputFile' => 'Required if paired is true and inputType is fastq. The second reads file',
            'inputType'       => 'Required, type of the input file (fastq, bam)',
            'convertBam'      => 'If inputType is bam converts input in another format: fastq.',
            'trimGalore'      => [
                'enable'  => 'A boolean value to indicate whether trim galore should run (This parameter works only for fastq files)',
                'quality' => 'Minimal PHREAD quality for trimming (Default 20)',
                'length'  => 'Minimal reads length (Default 14)',
            ],
        ];
    }

    /**
     * Get the common parameters laravel validation array
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return array
     */
    private static function commonParametersValidation(Request $request): array
    {
        $parameters = (array)$request->get('parameters', []);

        return [
            'paired'             => ['filled', 'boolean'],
            'firstInputFile'     => ['required', 'string'],
            'secondInputFile'    => [
                'nullable',
                Rule::requiredIf(
                    static function () use ($parameters) {
                        return $parameters['inputType'] === self::FASTQ && ((bool)($parameters['paired'] ?? false)) === true;
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
        ];
    }

    /**
     * Validate the common parameters of an input job
     *
     * @param \App\Models\Job $model
     * @param array           $validInputTypes
     * @param string          $fastQType
     *
     * @return bool
     */
    private function validateCommonParameters(Job $model, array $validInputTypes, string $fastQType): bool
    {
        $paired = (bool)$model->getParameter('paired', false);
        $inputType = $model->getParameter('inputType');
        if (!in_array($inputType, $validInputTypes, true) || !$this->validateFileParameter('firstInputFile')) {
            return false;
        }

        return !($paired && $inputType === $fastQType && !$this->validateFileParameter('secondInputFile'));
    }

}
