<?php

namespace App\Jobs;

use App\Models\Job as JobModel;
use Illuminate\Bus\Queueable;
use Illuminate\Contracts\Queue\ShouldQueue;
use Illuminate\Foundation\Bus\Dispatchable;
use Illuminate\Queue\InteractsWithQueue;
use Illuminate\Queue\SerializesModels;
use Throwable;

class DeleteJobDirectory implements ShouldQueue
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
    public function handle()
    {
        try {
            if ($this->model->exists && $this->model->canBeDeleted()) {
                $this->model->deleteJobDirectory();
                $this->model->delete();
            }
        } catch (Throwable $e) {
            $this->fail($e);
        }
    }
}
