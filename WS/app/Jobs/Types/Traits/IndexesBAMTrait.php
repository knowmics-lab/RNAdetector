<?php


namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Models\Job;

trait IndexesBAMTrait
{

    abstract protected function getFilePaths(string $filename): array;

    abstract protected function getFileNameWithoutExtension(string $filename): string;

    abstract public function log(string $text, bool $appendNewLine = true, bool $commit = true): void;

    /**
     * Indexes a BAM file and compute coverage
     *
     * @param \App\Models\Job $model
     * @param string          $bamFile
     * @param boolean         $unsorted
     * @param int             $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexBAM(
        Job $model,
        string $bamFile,
        bool $unsorted = false,
        int $threads = 1
    ): array {
        $outputFile = $this->getFilePaths(($unsorted) ? $this->getFileNameWithoutExtension($bamFile) . '_sorted.bam' : $bamFile);
        $command = [
            'bash',
            AbstractJob::scriptPath('prepare_bam.sh'),
            '-f',
            $outputFile[1],
        ];
        if ($unsorted) {
            $command[] = '-s';
            $command[] = '-u';
            $command[] = $bamFile;
            $command[] = '-t';
            $command[] = $threads;
        }
        AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                101 => 'Unsorted BAM file does not exist.',
                102 => 'Unable to sort BAM file',
                103 => 'BAM file does not exist.',
                104 => 'Unable to write BAM index file.',
                105 => 'Unable to compute BAM coverage.',
            ]
        );
        if (!file_exists($outputFile[1])) {
            throw new ProcessingJobException('Unable to validate output file');
        }

        return $outputFile;
    }

}
