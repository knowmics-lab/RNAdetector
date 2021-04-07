<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App;

use App\Jobs\Types\Factory;
use App\Models\Job;
use Cache;

final class SystemInfo
{

    /**
     * Find how much memory this machine has
     *
     * @return int
     */
    public function maxMemory(): int
    {
        if (Cache::has('availableMemory')) {
            return Cache::get('availableMemory');
        }
        $fh = @fopen('/proc/meminfo', 'rb');
        if (!$fh) {
            return -1;
        }
        $mem = -1;
        while ($line = fgets($fh)) {
            $pieces = [];
            if (preg_match('/^MemTotal:\s+(\d+)\skB$/', $line, $pieces)) {
                $mem = (int)$pieces[1];
                break;
            }
        }
        @fclose($fh);
        Cache::put('availableMemory', $mem, now()->addDay());

        return $mem;
    }

    /**
     * Count the number of cores available for this machine
     *
     * @return int
     */
    public function numCores(): int
    {
        if (Cache::has('numCores')) {
            return Cache::get('numCores');
        }
        $cpuInfo = file_get_contents('/proc/cpuinfo');
        preg_match_all('/^processor/m', $cpuInfo, $matches);
        $count = count($matches[0]);
        Cache::put('numCores', $count, now()->addDay());

        return $count;
    }

    /**
     * Detect how much cpu cores will be used by the last two analysis.
     * For the detection we take the last 10 jobs to find the maximum number of user threads.
     * Then it multiplies this values by 2 to get a rough estimate.
     * If all jobs are within boundaries (< 1/3 * CPU cores) the number is a correct estimation.
     *
     * @return int
     */
    public function usedCores(): int
    {
        $jobs = Job::whereIn('status', [Job::QUEUED, Job::PROCESSING])->latest()->take(10)->get();

        return (2 * $jobs->map(
                static function (Job $job) {
                    return Factory::threads($job);
                }
            )->max());
    }

    /**
     * Transforms this object in an array
     *
     * @return array
     */
    public function toArray(): array
    {
        return [
            'data' => [
                'containerVersion'       => Utils::VERSION,
                'containerVersionNumber' => Utils::VERSION_NUMBER,
                'maxMemory'              => $this->maxMemory(),
                'numCores'               => $this->numCores(),
                'usedCores'              => $this->usedCores(),
            ],
        ];
    }
}
