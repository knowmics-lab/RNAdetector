<?php
/**
 * Oncoreport Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Utils;
use Illuminate\Console\Command;
use Throwable;

class CheckForUpdate extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'update:check';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Checks if the update script should be started';


    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        $error = false;
        $message = '';
        $updateNeeded = false;
        $versionNumberFile = storage_path('app/version_number');
        if (file_exists($versionNumberFile)) {
            try {
                $content = json_decode(file_get_contents($versionNumberFile), true);
                $version = $content['version'] ?? Utils::DEFAULT_VERSION_NUMBER;
                $updateNeeded = Utils::VERSION_NUMBER > $version;
            } catch (Throwable $ex) {
                $error = true;
                $message = $ex->getMessage();
            }
        } else {
            $updateNeeded = true;
        }
        try {
            $this->line(
                json_encode(
                    [
                        'error'        => $error,
                        'message'      => $message,
                        'updateNeeded' => $updateNeeded,
                    ]
                )
            );
        } catch (Throwable $e) {
            $this->error($e->getMessage());

            return 100;
        }

        return 0;
    }
}
