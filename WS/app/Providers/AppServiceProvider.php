<?php

namespace App\Providers;

use Illuminate\Support\ServiceProvider;
use Queue;

class AppServiceProvider extends ServiceProvider
{
    /**
     * Register any application services.
     *
     * @return void
     */
    public function register(): void
    {
        //
    }

    /**
     * Bootstrap any application services.
     *
     * @return void
     */
    public function boot(): void
    {
        Queue::looping(
            static function () {
                $bootedFile = storage_path('app/booted');
                if (!file_exists($bootedFile)) {
                    @touch($bootedFile);
                    @chmod($bootedFile, 0777);
                }
            }
        );
    }
}
