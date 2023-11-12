<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App;

use App\Exceptions\CommandException;
use App\Exceptions\IgnoredException;
use App\Exceptions\ProcessingJobException;
use RecursiveDirectoryIterator;
use RecursiveIteratorIterator;
use SplFileInfo;
use Symfony\Component\Process\Exception\ProcessFailedException;
use Symfony\Component\Process\Process;
use Throwable;
use ZipArchive;

final class Packages
{

    public const REPO_URL = 'https://rnadetector.atlas.dmi.unict.it/repo/packages.json';

    /**
     * @var null|array
     */
    private $packages = null;

    private static function doesNeedUpdate(string $name, string $version): ?bool
    {
        $referenceDir = config('rnadetector.reference_path') . '/' . $name;
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

    public static function isPackageInstalled(string $name): bool
    {
        $referenceDir = config('rnadetector.reference_path') . '/' . $name;

        return file_exists($referenceDir) && is_dir($referenceDir) && file_exists($referenceDir . '/.installed');
    }

    private function fetch(): array
    {
        if ($this->packages === null) {
            $packages = json_decode(file_get_contents(self::REPO_URL), true);
            if (isset($packages['packages'])) {
                $packages['packages'] = array_combine(
                    array_map(
                        static function ($p) {
                            return $p['name'];
                        },
                        $packages['packages']
                    ),
                    $packages['packages']
                );
            }
            $this->packages = $packages;
        }

        return $this->packages;
    }

    /**
     * Filter packages removing installed ones
     *
     * @return array
     */
    public function fetchNotInstalled(): array
    {
        $packages = $this->fetch();
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
            usort(
                $packages['packages'],
                static function ($a, $b) {
                    // If one of the two packages needs an update, then it should be displayed before the others
                    if ($a['needsUpdate'] !== $b['needsUpdate']) {
                        return !$a['needsUpdate'] - !$b['needsUpdate'];
                    }

                    // If both packages need (or do not need) an update, then they will be sorted by name
                    return strcasecmp($a['name'], $b['name']);
                }
            );
        }

        return $packages;
    }

    public function fetchOne(string $name): ?array
    {
        $packages = $this->fetch();

        return $packages['packages'][$name] ?? null;
    }

    public function canBeUpdated(string $name): bool
    {
        $pkg = $this->fetchOne($name);

        return $pkg !== null && self::doesNeedUpdate($name, $pkg['version'] ?? '');
    }


}
