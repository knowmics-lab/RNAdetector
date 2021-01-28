<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Packages;
use App\Utils;
use Illuminate\Console\Command;

class InstallPackage extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'packages:install {name} {--local : install a local package}';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Install a package from the main repository';


    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        $name = $this->argument('name');
        $local = $this->option('local');
        if ($local) {
            $packageFile = env('REFERENCES_PATH') . '/' . $name . '.tar.bz2';
            $packageUrl = $packageFile;
            $packageMd5 = env('REFERENCES_PATH') . '/' . $name . '.tar.bz2.md5';
            if (!file_exists($packageFile)) {
                $this->error('The requested package does not exist locally.');

                return 1;
            }
            if (!file_exists($packageMd5)) {
                $this->error('Unable to find package checksum file.');

                return 2;
            }
        } else {
            $this->info('Fetching package details...');
            $packages = new Packages();
            $package = $packages->fetchOne($name);
            if ($package === null) {
                $this->error('Package not found in the repository.');

                return 3;
            }
            $packageUrl = $package['url'];
            $packageMd5 = $package['md5'];
        }
        if (!$packageUrl || !$packageMd5) {
            $this->error('Unknown error.');

            return 4;
        }
        Utils::runCommand(
            [
                'bash',
                realpath(env('BASH_SCRIPT_PATH') . '/install.package.sh'),
                '-n',
                $name,
                '-u',
                $packageUrl,
                '-m',
                $packageMd5,
            ],
            env('REFERENCES_PATH'),
            null,
            function ($type, $buffer) {
                $this->output->write('<info>' . $buffer . '</info>');
            }
        );

        return 0;
    }
}
