<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Controllers\Api;

use App\Http\Controllers\Controller;
use App\Http\Resources\Job as JobResource;
use App\Http\Resources\JobCollection;
use App\Http\Resources\UserCollection;
use App\Jobs\Request as JobRequest;
use App\Jobs\Types\Factory;
use App\Models\Job;
use App\Models\User;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;
use Illuminate\Support\Arr;
use Illuminate\Support\Str;
use Illuminate\Validation\Rule;

class JobController extends Controller
{

    /**
     * JobController constructor.
     */
    public function __construct()
    {
        $this->authorizeResource(Job::class, 'job');
    }


    /**
     * Display a listing of the resource.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return \App\Http\Resources\JobCollection
     */
    public function index(Request $request): JobCollection
    {
        $perPage = (int)($request->get('per_page') ?? 15);
        if ($perPage < 0) {
            $perPage = 15;
        }
        /** @var \App\Models\User $user */
        $user = \Auth::guard('api')->user();
        if ($user->admin) {
            return new JobCollection(Job::paginate($perPage));
        }

        return new JobCollection(Job::whereUserId($user->id)->paginate($perPage));
    }

    /**
     * Prepare array for nested validation
     *
     * @param array $specs
     *
     * @return array
     */
    private function _prepareNestedValidation(array $specs): array
    {
        $nestedSpecs = [];
        foreach ($specs as $field => $rules) {
            if (!Str::startsWith($field, 'parameters.')) {
                $field = 'parameters.' . $field;
            }
            $nestedSpecs[$field] = $rules;
        }

        return $nestedSpecs;
    }

    /**
     * Store a newly created resource in storage.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return \App\Http\Resources\Job
     * @throws \Illuminate\Validation\ValidationException
     */
    public function store(Request $request): JobResource
    {
        $jobTypes = Factory::listTypes();
        $validValues = $this->validate(
            $request,
            [
                'type'       => ['required', 'string', Rule::in($jobTypes->pluck('id'))],
                'parameters' => ['filled', 'array'],
            ]
        );
        $parametersValidation = $this->_prepareNestedValidation(
            Factory::validationSpec($validValues['type'], $request)
        );
        $validParameters = $this->validate($request, $parametersValidation);
        $type = $validValues['type'];
        $validParameters = $validParameters['parameters'] ?? [];
        $job = Job::create(
            [
                'job_type'       => $type,
                'status'         => Job::READY,
                'job_parameters' => [],
                'job_output'     => [],
                'log'            => '',
                'user_id'        => \Auth::guard('api')->id(),
            ]
        );
        $job->setParameters(Arr::dot($validParameters));
        $job->save();

        return new JobResource($job);
    }

    /**
     * Display the specified resource.
     *
     * @param \App\Models\Job $job
     *
     * @return \App\Http\Resources\Job
     */
    public function show(Job $job): JobResource
    {
        return new JobResource($job);
    }

    /**
     * Update the specified resource in storage.
     *
     * @param \Illuminate\Http\Request $request
     * @param \App\Models\Job          $job
     *
     * @return \App\Http\Resources\Job
     * @throws \Illuminate\Validation\ValidationException
     */
    public function update(Request $request, Job $job): JobResource
    {
        if (!$job->canBeModified()) {
            abort(400, 'Unable to modify a submitted, running, or completed job.');
        }
        $parametersValidation = $this->_prepareNestedValidation(Factory::validationSpec($job, $request));
        $validParameters = $this->validate($request, $parametersValidation);
        $validParameters = $validParameters['parameters'] ?? [];
        $job->addParameters(Arr::dot($validParameters));
        $job->save();

        return new JobResource($job);
    }

    /**
     * Remove the specified resource from storage.
     *
     * @param \App\Models\Job $job
     *
     * @return \Illuminate\Http\JsonResponse
     * @throws \Exception
     */
    public function destroy(Job $job): JsonResponse
    {
        if (!$job->canBeDeleted()) {
            abort(400, 'Unable to delete a queued or running job.');
        }
        $job->deleteJobDirectory();
        $job->delete();

        return response()->json(
            [
                'message' => 'Job deleted.',
                'errors'  => false,
            ]
        );
    }

    /**
     * Submit the specified resource for execution
     *
     * @param \App\Models\Job $job
     *
     * @return \App\Http\Resources\Job
     */
    public function submit(Job $job): JobResource
    {
        if (!$job->canBeModified()) {
            abort(400, 'Unable to submit a job that is already submitted.');
        }
        $job->setStatus(Job::QUEUED);
        JobRequest::dispatch($job);

        return new JobResource($job);
    }

    /**
     * Upload a file to the specified job
     *
     * @param \App\Models\Job $job
     *
     * @return mixed
     */
    public function upload(Job $job)
    {
        if (!$job->canBeModified()) {
            abort(400, 'Unable to upload a file for a job that is already submitted.');
        }
        set_time_limit(0);
        /** @var \TusPhp\Tus\Server $server */
        $server = app('tus-server');
        $server->setApiPath(route('jobs.upload', $job, false))
               ->setUploadDir($job->getAbsoluteJobDirectory());
        $response = $server->serve();

        return $response->send();
    }
}
