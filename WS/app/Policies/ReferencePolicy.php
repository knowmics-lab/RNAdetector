<?php

namespace App\Policies;

use App\Models\User;
use App\Models\Reference;
use Illuminate\Auth\Access\HandlesAuthorization;

class ReferencePolicy
{
    use HandlesAuthorization;

    /**
     * Determine whether the user can view any references.
     *
     * @param  \App\Models\User  $user
     * @return mixed
     */
    public function viewAny(User $user)
    {
        return true;
    }

    /**
     * Determine whether the user can view the reference.
     *
     * @param  \App\Models\User  $user
     * @param  \App\Models\Reference  $reference
     * @return mixed
     */
    public function view(User $user, Reference $reference)
    {
        return  true;
    }

    /**
     * Determine whether the user can create references.
     *
     * @param  \App\Models\User  $user
     * @return mixed
     */
    public function create(User $user)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can update the reference.
     *
     * @param  \App\Models\User  $user
     * @param  \App\Models\Reference  $reference
     * @return mixed
     */
    public function update(User $user, Reference $reference)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can delete the reference.
     *
     * @param  \App\Models\User  $user
     * @param  \App\Models\Reference  $reference
     * @return mixed
     */
    public function delete(User $user, Reference $reference)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can restore the reference.
     *
     * @param  \App\Models\User  $user
     * @param  \App\Models\Reference  $reference
     * @return mixed
     */
    public function restore(User $user, Reference $reference)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can permanently delete the reference.
     *
     * @param  \App\Models\User  $user
     * @param  \App\Models\Reference  $reference
     * @return mixed
     */
    public function forceDelete(User $user, Reference $reference)
    {
        return $user->admin;
    }
}
