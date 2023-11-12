<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Models\Reference;
use App\Packages;
use App\Utils;
use Illuminate\Http\Request;
use Storage;

class InstallPackagesJobType extends AbstractJob
{

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'names' => 'A list of package names that will be installed',
        ];
    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [
            'done' => 'A boolean that is true if the job has been successfully completed',
        ];
    }

    /**
     * Returns an array containing rules for input validation.
     *
     * @param  \Illuminate\Http\Request  $request
     *
     * @return array
     */
    public static function validationSpec(Request $request): array
    {
        return [
            'names'   => ['required', 'array'],
            'names.*' => ['required', 'alpha_dash', 'max:255'],
        ];
    }

    /**
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \Throwable
     */
    public function handle(): void
    {
        $this->log('Starting job.');
        $names = $this->getParameter('names', []);

        foreach ($names as $name) {
            $this->log('Fetching details for package ' . $name . '...');
            $packages = new Packages();
            $package = $packages->fetchOne($name);
            throw_if($package === null, ProcessingJobException::class, 'Package ' . $name . ' was not found in the repository.');
            $packageUrl = $package['url'];
            $packageMd5 = $package['md5'];
            throw_if(!$packageUrl || !$packageMd5, ProcessingJobException::class, 'Unknown error');
            $this->log('Installing ' . $name . '...');
            self::runCommand(
                [
                    'bash',
                    self::scriptPath('install.package.sh'),
                    '-n',
                    $name,
                    '-u',
                    $packageUrl,
                    '-m',
                    $packageMd5,
                ],
                $this->getAbsoluteJobDirectory(),
                null,
                function ($type, $buffer) {
                    $this->log($buffer, false);
                },
                [
                    3 => 'Package url is required',
                    4 => 'Package name is required',
                    5 => 'Unable to download the package',
                    6 => 'Unable to download the package checksum',
                    7 => 'Checksum control failed',
                    8 => 'Unable to install the package',
                ]
            );
            $this->log('Package ' . $name . ' installed!');
        }

        $this->model->setOutput(['type' => self::OUT_TYPE_CONFIRMATION, 'done' => true]);
        $this->log('Job completed!');
    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Install a package from the repository';
    }

    /**
     * @inheritDoc
     */
    public static function displayName(): string
    {
        return 'Install Package';
    }
}
