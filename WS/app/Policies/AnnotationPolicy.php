<?php

namespace App\Policies;

use App\Models\Annotation;
use App\Models\User;
use Illuminate\Auth\Access\HandlesAuthorization;

class AnnotationPolicy
{
    use HandlesAuthorization;

    /**
     * Determine whether the user can view any annotations.
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
     * Determine whether the user can view the annotation.
     *
     * @param \App\Models\User       $user
     * @param \App\Models\Annotation $annotation
     *
     * @return mixed
     */
    public function view(User $user, Annotation $annotation)
    {
        return true;
    }

    /**
     * Determine whether the user can create annotations.
     *
     * @param \App\Models\User $user
     *
     * @return mixed
     */
    public function create(User $user)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can update the annotation.
     *
     * @param \App\Models\User       $user
     * @param \App\Models\Annotation $annotation
     *
     * @return mixed
     */
    public function update(User $user, Annotation $annotation)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can delete the annotation.
     *
     * @param \App\Models\User       $user
     * @param \App\Models\Annotation $annotation
     *
     * @return mixed
     */
    public function delete(User $user, Annotation $annotation)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can restore the annotation.
     *
     * @param \App\Models\User       $user
     * @param \App\Models\Annotation $annotation
     *
     * @return mixed
     */
    public function restore(User $user, Annotation $annotation)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can permanently delete the annotation.
     *
     * @param \App\Models\User       $user
     * @param \App\Models\Annotation $annotation
     *
     * @return mixed
     */
    public function forceDelete(User $user, Annotation $annotation)
    {
        return $user->admin;
    }
}
