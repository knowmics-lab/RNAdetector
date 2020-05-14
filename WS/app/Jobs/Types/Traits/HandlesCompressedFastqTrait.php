<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\IgnoredException;
use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Models\Job;
use App\Utils;
use Storage;

trait HandlesCompressedFastqTrait
{

    /**
     * Checks if a FASTQ file is compressed and extracts its content
     *
     * @param \App\Models\Job $model
     * @param string          $fastqFile
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private static function checksForCompression(
        Job $model,
        ?string $fastqFile
    ): string {
        if ($fastqFile === null) {
            return $fastqFile;
        }
        try {
            AbstractJob::runCommand(
                [
                    'bash',
                    AbstractJob::scriptPath('prepare_fq.sh'),
                    $fastqFile,
                ],
                $model->getAbsoluteJobDirectory(),
                null,
                static function ($type, $buffer) use ($model) {
                    $model->appendLog($buffer, false);
                },
                [
                    1 => 'Input file is required.',
                    2 => 'Input file does not exist.',
                    3 => 'Extracted output file does not exist.',
                    4 => Utils::IGNORED_ERROR_CODE,
                ]
            );
        } catch (IgnoredException $e) {
            if ($e->getCode() === 4) {
                $disk = Storage::disk('public');
                $path = pathinfo($fastqFile, PATHINFO_FILENAME);
                if (!$disk->exists($model->getJobDirectory() . '/' . $path)) {
                    throw new ProcessingJobException('Unable to find extracted output file (' . $path . ').');
                }

                return $path;
            }
        }

        return $fastqFile;
    }

}
