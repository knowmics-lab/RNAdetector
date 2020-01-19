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
use App\Models\Job;
use App\Models\Reference;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class SamplesGroupJobType extends AbstractJob
{
    use HasCommonParameters, ConvertsSamToBamTrait, RunTrimGaloreTrait, UseAlignmentTrait, UseCountingTrait;

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'jobs'        => 'A list of analysis job of the same type',
            'description' => 'An optional tsv file containing samples descriptions',
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
            'jobs'        => 'A list of the valid job identifiers',
            'codes'       => 'A list of the valid job sample codes',
            'description' => 'An optional path/url of the samples descriptions file (filtered by valid codes)',
            'metadata'    => 'An optional list of available metadata for this sample group (columns of the description file)',
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
            'jobs'        => ['required', 'array', Rule::exists('jobs', 'id')],
            'description' => ['filled', 'string'],
        ];
    }

    /**
     * @inheritDoc
     */
    public function isInputValid(): bool
    {
        $description = $this->getParameter('description', null);
        $disk = Storage::disk('public');
        $dir = $this->getJobDirectory() . '/';

        return !(!empty($description) && !$disk->exists($dir . $description));
    }

    /**
     * Find all valid jobs
     *
     * @param array $jobs
     *
     * @return \App\Models\Job[]
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function processValidJobs(array $jobs): array
    {
        /** @var Job[] $models */
        $models = array_filter(
            array_map(
                static function ($job) {
                    return Job::whereId($job)->first();
                },
                $jobs
            )
        );
        if (!count($models)) {
            throw new ProcessingJobException('No jobs have been specified');
        }

        return $models;
    }

    /**
     * Checks if all job types are the same
     *
     * @param \App\Models\Job[] $models
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkJobsTypes(array $models): void
    {
        $firstType = $models[0]->job_type;
        foreach ($models as $model) {
            if ($model->job_type !== $firstType) {
                throw new ProcessingJobException('All jobs must be of the same type');
            }
        }
    }

    /**
     * Pull a property from an array of models
     *
     * @param \App\Models\Job[] $models
     * @param string            $property
     *
     * @return array
     */
    private function pullProperty(array $models, string $property): array
    {
        return array_map(
            static function ($job) use ($property) {
                return $job->$property;
            },
            $models
        );
    }

    /**
     * Build a description file
     *
     * @param \App\Models\Job[] $models
     *
     * @return array
     */
    private function makeDescriptionFile(array $models): array
    {
        $descriptionRelative = $this->getJobFile('description_', '.tsv');
        $descriptionFile = $this->absoluteJobPath($descriptionRelative);
        $descriptionUrl = Storage::disk('public')->url($descriptionRelative);
        $content = "SampleId\tSampleGroup" . PHP_EOL . implode(
                PHP_EOL,
                array_map(
                    function ($job) {
                        return $job->sample_code . "\t" . $this->model->sample_code;
                    },
                    $models
                )
            );
        @file_put_contents($descriptionFile, $content);
        @chmod($descriptionFile, 0777);

        return [$descriptionRelative, $descriptionUrl, ['SampleGroup']];
    }

    /**
     * @param \App\Models\Job[] $models
     * @param string[]          $validCodes
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function filterDescriptionFile(array $models, array $validCodes): array
    {
        $inputDescription = $this->model->getParameter('description', null);
        $dir = $this->getJobDirectory() . '/';
        if (empty($inputDescription) || !Storage::disk('public')->exists($dir . $inputDescription)) {
            return $this->makeDescriptionFile($models);
        }
        $descriptionRelative = $this->getJobFile('description_', '.tsv');
        $descriptionFile = $this->absoluteJobPath($descriptionRelative);
        $descriptionUrl = Storage::disk('public')->url($descriptionRelative);
        $metas = [];
        $inputFp = fopen($inputDescription, 'rb');
        $outputFp = fopen($descriptionFile, 'wb');
        if (!$inputFp) {
            throw new ProcessingJobException('Unable to open description file.');
        }
        if (!$outputFp) {
            fclose($inputFp);
            throw new ProcessingJobException('Unable to open output file.');
        }
        $firstLine = true;
        while (($data = fgetcsv($inputFp, 0, "\t")) !== false) {
            $firstElement = array_shift($data); // The first element must be always the sample identifier
            array_unshift($data, ($firstLine ? 'SampleGroup' : $this->model->sample_code));
            if ($firstLine) {
                $metas = $data;
                $firstLine = false;
            } elseif (!in_array($firstElement, $validCodes, true)) {
                continue;
            }
            array_unshift($data, $firstElement);
            fputcsv($outputFp, $data, "\t");
        }
        fclose($outputFp);
        fclose($inputFp);

        return [$descriptionRelative, $descriptionUrl, $metas];
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
        $jobs = $this->getParameter('jobs', []);
        $models = $this->processValidJobs($jobs);
        $this->checkJobsTypes($models);
        /** @var int[] $validJobs */
        $validJobs = $this->pullProperty($models, 'id');
        /** @var string[] $validCodes */
        $validCodes = $this->pullProperty($models, 'sample_code');
        [$descriptionRelative, $descriptionUrl, $metadata] = $this->filterDescriptionFile($models, $validCodes);
        $this->model->setOutput(
            [
                'jobs'        => $validJobs,
                'codes'       => $validCodes,
                'description' => ['path' => $descriptionRelative, 'url' => $descriptionUrl],
                'metadata'    => $metadata,
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
        return 'Groups several jobs of the same type together (ideally a group of samples)';
    }
}
