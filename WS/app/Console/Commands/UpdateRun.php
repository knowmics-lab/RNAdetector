<?php
/**
 * Oncoreport Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use Illuminate\Console\Command;
use Symfony\Component\Process\Process;

class UpdateRun extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'update:run';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Runs the update script';

    /**
     * Execute an update script
     *
     * @param string $script
     * @param bool   $warn
     *
     * @return int
     */
    protected function runScript(string $script, bool $warn = false): int
    {
        $updateScript = AbstractJob::scriptPath($script);
        if (file_exists($updateScript)) {
            try {
                AbstractJob::runCommand(
                    [
                        'bash',
                        $updateScript,
                    ],
                    null,
                    null,
                    function ($type, $buffer) {
                        if ($type === Process::ERR) {
                            $this->error($buffer);
                        } else {
                            $this->line($buffer);
                        }
                    }
                );
            } catch (ProcessingJobException $e) {
                $this->error($e->getMessage());

                return 100;
            }
        } elseif ($warn) {
            $this->warn('No update script to run');
        }

        return 0;
    }


    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        $retCode = $this->runScript('pre_update.bash');
        if ($retCode !== 0) {
            return $retCode;
        }
        $retCode = $this->call(
            'migrate',
            [
                '--force'          => true,
                '--no-interaction' => true,
            ]
        );
        if ($retCode !== 0) {
            return $retCode;
        }

        return $this->runScript('post_update.bash', true);
    }
}
