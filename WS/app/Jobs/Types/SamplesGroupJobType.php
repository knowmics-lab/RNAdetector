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
            'jobs' => 'A list of analysis job of the same type',
        ];
    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [];
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
            'jobs' => ['required', 'array', Rule::exists('jobs', 'id')],
        ];
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
        /** @var Job[] $models */
        $models = array_filter(array_map(static function ($job) {
            return Job::whereId($job)->first();
        }, $jobs));
        if (!count($models)) {
            throw new ProcessingJobException('No jobs have been specified');
        }
        $firstType = $models[0]->job_type;
        foreach ($models as $model) {
            if ($model->job_type !== $firstType) {
                throw new ProcessingJobException('All jobs must be of the same type');
            }
        }
        $this->model->setOutput([]);
        $this->model->save();
    }


    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Groups several jobs together';
    }
}
