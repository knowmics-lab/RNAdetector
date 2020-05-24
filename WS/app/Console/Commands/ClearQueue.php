<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Models\Job;
use App\Utils;
use Illuminate\Console\Command;
use Queue;
use Symfony\Component\Process\Process;

class ClearQueue extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'queue:clear';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Clears the queue before closing the docker container';


    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        while (($j = Queue::pop()) !== null) {
            $j->delete();
        }
        $this->call('queue:flush');
        foreach (Job::whereStatus(Job::QUEUED)->get() as $job) {
            $job->status = Job::READY;
            $job->save();
        }
        foreach (Job::whereStatus(Job::PROCESSING)->get() as $job) {
            $job->status = Job::FAILED;
            $job->appendLog('Job failed since queue was cleared!', true, false);
            $job->save();
        }

        return 0;
    }
}
