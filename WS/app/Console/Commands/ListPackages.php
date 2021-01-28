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

    private static function doesNeedUpdate(string $name, string $version): ?bool
    {
        $referenceDir = env('REFERENCES_PATH') . '/' . $name;
        $versionFile = $referenceDir . '/.version';
        if ($version === '' || !self::isPackageInstalled($name)) {
            return false;
        }
        if (!file_exists($versionFile)) {
            return true;
        }
        $currentVersion = file_get_contents($versionFile);

        return (strcasecmp($currentVersion, $version) !== 0);
    }

    private static function isPackageInstalled(string $name): bool
    {
        $referenceDir = env('REFERENCES_PATH') . '/' . $name;

        return file_exists($referenceDir) && is_dir($referenceDir) && file_exists($referenceDir . '/.installed');
    }

    /**
     * Filter packages removing installed ones
     *
     * @param array $packages
     *
     * @return array
     */
    public static function filterPackages(array $packages): array
    {
        if (isset($packages['packages'])) {
            $packages['packages'] = array_values(
                array_filter(
                    array_map(
                        static function ($pkg) {
                            $pkg['needsUpdate'] = self::doesNeedUpdate($pkg['name'], $pkg['version'] ?? '');

                            return $pkg;
                        },
                        $packages['packages']
                    ),
                    static function ($pkg) {
                        return !self::isPackageInstalled($pkg['name']) || $pkg['needsUpdate'];
                    }
                )
            );
        }

        return $packages;
    }

    /**
     * Fetch available packages from the repository URL
     *
     * @return array
     */
    public static function fetchPackages(): array
    {
        $url = 'https://rnadetector.atlas.dmi.unict.it/repo/packages.json';
        $content = json_decode(file_get_contents($url), true);

        return self::filterPackages($content);
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        try {
            $this->line(json_encode(self::fetchPackages()));
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
