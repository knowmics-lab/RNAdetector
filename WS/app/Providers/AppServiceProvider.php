<?php

namespace App\Providers;

use App\Utils;
use DB;
use Exception;
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
                $check = false;
                while (!$check) {
                    try {
                        DB::connection()->getPdo();
                        $check = true;
                    } catch (Exception $ignore) {
                        sleep(5);
                    }
                }
                $bootedFile = storage_path('app/booted');
                if (!file_exists($bootedFile)) {
                    @touch($bootedFile);
                    @chmod($bootedFile, 0777);
                }
                $versionNumberFile = storage_path('app/version_number');
                if (!file_exists($versionNumberFile)) {
                    @file_put_contents(
                        $versionNumberFile,
                        json_encode(
                            [
                                'version' => Utils::VERSION_NUMBER,
                            ]
                        )
                    );
                    @chmod($versionNumberFile, 0644);
                }
            }
        );
    }
}
