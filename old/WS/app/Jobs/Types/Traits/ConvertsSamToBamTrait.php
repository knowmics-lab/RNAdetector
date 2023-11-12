<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Models\Job;
use Symfony\Component\Process\Process;

trait ConvertsSamToBamTrait
{

    /**
     * Converts BAM file to fastq file
     *
     * @param \App\Models\Job $model
     * @param string          $samFile
     * @param bool            $unique
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private static function convertSamToBam(
        Job $model,
        string $samFile,
        bool $unique = true
    ): string {
        $model->appendLog('Converting SAM to BAM.');
        $bamFile = ($unique) ? $model->getJobTempFileAbsolute('sam2bam_', '.bam') : $model->getJobFileAbsolute('sam2bam_', '.bam');
        $output = AbstractJob::runCommand(
            [
                'bash',
                AbstractJob::scriptPath('sam2bam.sh'),
                '-i',
                $samFile,
                '-o',
                $bamFile,
            ],
            $model->getAbsoluteJobDirectory(),
            null,
            static function ($type, $buffer) use ($model) {
                $model->appendLog($buffer, false);
            },
            [
                3 => 'Input file does not exist.',
                4 => 'Output file must be specified.',
                5 => 'Output directory is not writable.',
                6 => 'An error occurred during samtools execution.',
            ]
        );
        if (!file_exists($bamFile)) {
            throw new ProcessingJobException('Unable to convert sam to bam.');
        }
        // $model->appendLog($output);
        $model->appendLog('SAM converted to BAM.');

        return $bamFile;
    }

}
