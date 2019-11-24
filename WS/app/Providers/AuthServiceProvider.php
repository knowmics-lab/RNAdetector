<?php

namespace App\Providers;

use App\Models\Annotation;
use App\Models\Job;
use App\Models\Reference;
use App\Models\User;
use App\Policies\AnnotationPolicy;
use App\Policies\JobPolicy;
use App\Policies\ReferencePolicy;
use App\Policies\UserPolicy;
use Illuminate\Foundation\Support\Providers\AuthServiceProvider as ServiceProvider;

class AuthServiceProvider extends ServiceProvider
{
    /**
     * The policy mappings for the application.
     *
     * @var array
     */
    protected $policies = [
        User::class       => UserPolicy::class,
        Job::class        => JobPolicy::class,
        Reference::class  => ReferencePolicy::class,
        Annotation::class => AnnotationPolicy::class,
    ];

    /**
     * Register any authentication / authorization services.
     *
     * @return void
     */
    public function boot()
    {
        $this->registerPolicies();

        //
    }
}
