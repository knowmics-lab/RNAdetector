<?php

namespace App\Providers;

use TusPhp\Tus\Server as TusServer;
use Illuminate\Support\ServiceProvider;

class TusServiceProvider extends ServiceProvider
{
    /**
     * Register services.
     *
     * @return void
     */
    public function register()
    {
        $this->app->singleton(
            'tus-server',
            static function ($app) {
                $server = new TusServer();

                /*$server->setApiPath('/api/tus') // tus server endpoint.
                       ->setUploadDir(storage_path('app/public/uploads')); // uploads dir.*/

                return $server;
            }
        );
    }

    /**
     * Bootstrap services.
     *
     * @return void
     */
    public function boot()
    {
        //
    }
}
