<?php

namespace App\Providers;

use Illuminate\Support\ServiceProvider;

class AppServiceProvider extends ServiceProvider
{
    /**
     * Register any application services.
     *
     * @return void
     */
    public function register()
    {
        //
    }

    /**
     * Bootstrap any application services.
     *
     * @return void
     */
    public function boot()
    {
        $bootedFile = storage_path('app/booted');
        if (!file_exists($bootedFile)) {
            @touch($bootedFile);
            @chmod($bootedFile, 0777);
        }
    }
}
