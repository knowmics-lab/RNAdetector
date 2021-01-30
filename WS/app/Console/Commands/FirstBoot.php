<?php
/**
 * Oncoreport Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Utils;
use Illuminate\Console\Command;
use Symfony\Component\Process\Process;

class FirstBoot extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'first:boot';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Runs some commands that are required on first boot';

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
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

        return 0;
    }
}
