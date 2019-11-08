<?php


namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Models\Job;

trait ConvertsBamToFastqTrait
{

    /**
     * Converts BAM file to fastq file
     *
     * @param \App\Models\Job $model
     * @param bool            $paired
     * @param string          $firstInputFile
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private static function convertBamToFastq(Job $model, bool $paired, string $firstInputFile): array
    {
        $firstFastQ = $model->getJobTempFileAbsolute('bam2fastq_', '.fastq');
        $secondFastQ = ($paired) ? $model->getJobTempFileAbsolute('bam2fastq_', '.fastq') : null;
        //call conversion script
        if (!$paired && !file_exists($firstFastQ)) {
            throw new ProcessingJobException('Unable to convert bam to fastq.');
        }
        if ($paired && !file_exists($firstFastQ) && !file_exists($secondFastQ)) {
            throw new ProcessingJobException('Unable to convert bam to fastq.');
        }

        return [$firstFastQ, $secondFastQ];
    }

}
