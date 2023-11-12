<?php
/**
 * Oncoreport Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\SystemInfo;
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
        $updateNeeded = false;
        $message = '';
        try {
            $updateNeeded = (new SystemInfo())->isUpdateNeeded();
        } catch (Throwable $ex) {
            $error = true;
            $message = $ex->getMessage();
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
