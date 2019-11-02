<?php

namespace App\Policies;

use App\Models\Job;
use App\Models\User;
use Illuminate\Auth\Access\HandlesAuthorization;

class JobPolicy
{
    use HandlesAuthorization;

    /**
     * Determine whether the user can view any jobs.
     *
     * @param \App\Models\User $user
     *
     * @return mixed
     */
    public function viewAny(User $user)
    {
        return true;
    }

    /**
     * Determine whether the user can view the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function view(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can create jobs.
     *
     * @param \App\Models\User $user
     *
     * @return mixed
     */
    public function create(User $user)
    {
        return true;
    }

    /**
     * Determine whether the user can update the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function update(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can delete the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function delete(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can restore the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function restore(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can permanently delete the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function forceDelete(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can submit the job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function submitJob(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }

    /**
     * Determine whether the user can upload files to a job.
     *
     * @param \App\Models\User $user
     * @param \App\Models\Job  $job
     *
     * @return mixed
     */
    public function uploadJob(User $user, Job $job)
    {
        return $user->admin || $job->user_id === $user->id;
    }
}
