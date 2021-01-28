<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Packages;
use Cache;
use Illuminate\Console\Command;
use Throwable;

class ListPackages extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'packages:list';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'List all packages available in the main repository';

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        try {
            $packages = new Packages();
            $this->line(json_encode($packages->fetchNotInstalled()));
        } catch (Throwable $e) {
            $this->line(
                json_encode(
                    [
                        'error'   => true,
                        'message' => $e->getMessage(),
                    ]
                )
            );

            return 1;
        }

        return 0;
    }
}
