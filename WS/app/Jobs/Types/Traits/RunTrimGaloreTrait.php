<?php


namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Models\Job;

trait RunTrimGaloreTrait
{
    /**
     * Call trimGalore using input parameters
     *
     * @param \App\Models\Job $model
     * @param bool            $paired
     * @param string          $firstInputFile
     * @param string|null     $secondInputFile
     * @param int             $quality
     * @param int             $length
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private static function runTrimGalore(
        Job $model,
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile = null,
        int $quality = 20,
        $length = 14
    ): array {
        $outputDirectory = $model->getJobTempFileAbsolute('trim_galore_');
        $command = [
            'bash',
            env('BASH_SCRIPT_PATH') . '/trim_galore.bash',
            '-q',
            $quality,
            '-l',
            $length,
            '-o',
            $outputDirectory,
            '-f',
            $firstInputFile,
        ];
        if ($paired) {
            $command[] = '-s';
            $command[] = $secondInputFile;
        }
        $output = AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Input file does not exist.',
                4 => 'Second input file does not exist.',
                5 => 'Output directory must be specified.',
                6 => 'Output directory is not writable.',
            ]
        );
        if (!file_exists($outputDirectory) || !is_dir($outputDirectory)) {
            throw new ProcessingJobException('Unable to create trimGalore output folder');
        }
        if ($paired) {
            $firstOutput = $outputDirectory . '/' . basename($firstInputFile, '.fastq') . '_val_1.fq';
            $secondOutput = $outputDirectory . '/' . basename($secondInputFile, '.fastq') . '_val_2.fq';
            if (!file_exists($firstOutput) || !file_exists($secondOutput)) {
                throw new ProcessingJobException('Unable to create output files');
            }
        } else {
            $firstOutput = $outputDirectory . '/' . basename($firstInputFile, '.fastq') . '_trimmed.fq';
            $secondOutput = null;
            if (!file_exists($firstOutput)) {
                throw new ProcessingJobException('Unable to create output files');
            }
        }

        return [$firstOutput, $secondOutput, $output];
    }
}
