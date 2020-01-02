<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs;

use App\Exceptions\ProcessingJobException;
use Auth;
use Illuminate\Bus\Queueable;
use Illuminate\Contracts\Queue\ShouldQueue;
use Illuminate\Foundation\Bus\Dispatchable;
use Illuminate\Queue\InteractsWithQueue;
use Illuminate\Queue\SerializesModels;
use App\Models\Job as JobModel;
use Throwable;

class Request implements ShouldQueue
{
    use Dispatchable, InteractsWithQueue, Queueable, SerializesModels;

    /**
     * Delete the job if its models no longer exist.
     *
     * @var bool
     */
    public $deleteWhenMissingModels = true;

    public $timeout = 0;

    /**
     * @var \App\Models\Job
     */
    protected $model;

    /**
     * Create a new job instance.
     *
     * @param \App\Models\Job $model
     */
    public function __construct(JobModel $model)
    {
        $this->model = $model;
    }

    /**
     * Execute the job.
     *
     * @return void
     */
    public function handle(): void
    {
        $jobProcessor = null;
        try {
            if ($this->model->shouldNotRun()) {  // job is being processed (or has been processed) by another thread.
                $this->delete();

                return;
            }
            $this->model->log = '';
            $this->model->setStatus(JobModel::PROCESSING);
            $this->delete();
            $jobProcessor = Types\Factory::get($this->model);
            if (!$jobProcessor->isInputValid()) {
                throw new ProcessingJobException('Job input format is not valid');
            }
            Auth::login($this->model->user);
            $jobProcessor->handle();
            Auth::logout();
            $this->model->setStatus(JobModel::COMPLETED);
        } catch (Throwable $e) {
            $this->model->appendLog('Error: ' . $e);
            $this->model->setStatus(JobModel::FAILED);
            if ($jobProcessor !== null && ($jobProcessor instanceof Types\AbstractJob)) {
                $jobProcessor->cleanupOnFail();
            }
            $this->fail($e);
        }
    }
}
