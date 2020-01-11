<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

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
     * Fetch available packages from the repository URL
     *
     * @return array
     */
    public static function fetchPackages(): array
    {
        $key = 'repo-packages';
        $url = 'https://rnadetector.atlas.dmi.unict.it/repo/packages.json';
        if (Cache::has($key)) {
            return Cache::get($key);
        }

        $content = json_decode(file_get_contents($url), true);
        Cache::put($key, $content, 24 * 60 * 60); // store for 1 day

        return $content;
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        try {
            $this->line(json_encode(self::fetchPackages()));;
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
