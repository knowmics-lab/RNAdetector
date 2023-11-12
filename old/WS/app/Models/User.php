<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

use Illuminate\Database\Eloquent\Builder;
use Illuminate\Database\Eloquent\Relations\HasMany;
use Illuminate\Foundation\Auth\User as Authenticatable;
use Illuminate\Notifications\Notifiable;

class User extends Authenticatable
{
    use Notifiable;

    /**
     * The attributes that are mass assignable.
     *
     * @var array
     */
    protected $fillable = [
        'name',
        'email',
        'password',
        'api_token',
        'admin',
    ];

    /**
     * The attributes that should be hidden for arrays.
     *
     * @var array
     */
    protected $hidden = [
        'password',
        'remember_token',
        'api_token',
    ];

    /**
     * The attributes that should be cast to native types.
     *
     * @var array
     */
    protected $casts = [
        'email_verified_at' => 'datetime',
        'admin'             => 'boolean',
    ];

    /**
     * @return \Illuminate\Database\Eloquent\Relations\HasMany
     */
    public function jobs(): HasMany
    {
        return $this->hasMany(Job::class, 'user_id', 'id');
    }

    /**
     * Returns some statistics about the user or the system if the user is an administrator.
     *
     * @return array
     */
    public function statistics(): array
    {
        $stats = [];
        if ($this->admin) {
            $stats['jobs'] = [
                'all'        => Job::count(),
                'ready'      => Job::whereStatus(Job::READY)->count(),
                'queued'     => Job::whereStatus(Job::QUEUED)->count(),
                'processing' => Job::whereStatus(Job::PROCESSING)->count(),
                'failed'     => Job::whereStatus(Job::FAILED)->count(),
                'completed'  => Job::whereStatus(Job::COMPLETED)->count(),
            ];
        } else {
            $stats['jobs'] = [
                'all'        => Job::whereUserId($this->id)->count(),
                'ready'      => Job::whereUserId($this->id)->whereStatus(Job::READY)->count(),
                'queued'     => Job::whereUserId($this->id)->whereStatus(Job::QUEUED)->count(),
                'processing' => Job::whereUserId($this->id)->whereStatus(Job::PROCESSING)->count(),
                'failed'     => Job::whereUserId($this->id)->whereStatus(Job::FAILED)->count(),
                'completed'  => Job::whereUserId($this->id)->whereStatus(Job::COMPLETED)->count(),
            ];
        }

        return $stats;
    }

}
