<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Policies;

use App\Models\User;
use Illuminate\Auth\Access\HandlesAuthorization;

class UserPolicy
{
    use HandlesAuthorization;

    /**
     * Determine whether the user can view any models.
     *
     * @param \App\Models\User $user
     * @return mixed
     */
    public function viewAny(User $user)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can view the model.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function view(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }

    /**
     * Determine whether the user can create models.
     *
     * @param \App\Models\User $user
     * @return mixed
     */
    public function create(User $user)
    {
        return $user->admin;
    }

    /**
     * Determine whether the user can update the model.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function update(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }

    /**
     * Determine whether the user can delete the model.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function delete(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }

    /**
     * Determine whether the user can restore the model.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function restore(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }

    /**
     * Determine whether the user can permanently delete the model.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function forceDelete(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }

    /**
     * Determine whether the user can generate the access token.
     *
     * @param \App\Models\User $user
     * @param \App\Models\User $model
     * @return mixed
     */
    public function generateToken(User $user, User $model)
    {
        return $user->admin || $user->id === $model->id;
    }
}
