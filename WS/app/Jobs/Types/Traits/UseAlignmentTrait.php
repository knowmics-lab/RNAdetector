<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\AbstractJob;
use App\Models\Annotation;
use App\Models\Job;
use App\Models\Reference;
use Storage;

trait UseAlignmentTrait
{

    /**
     * Runs TopHat
     *
     * @param \App\Models\Job        $model
     * @param bool                   $paired
     * @param string                 $firstInputFile
     * @param string|null            $secondInputFile
     * @param \App\Models\Reference  $genome
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runTophat(
        Job $model,
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        Reference $genome,
        Annotation $annotation,
        int $threads = 1
    ): string {
        if ($genome->isAvailableFor('tophat')) {
            throw new ProcessingJobException('The selected reference has not been indexed for tophat.');
        }
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $model->appendLog('Aligning reads using TopHat.');
        $bamOutput = $model->getJobTempFileAbsolute('tophat_output', '.bam');
        $command = [
            'bash',
            AbstractJob::scriptPath('tophat.bash'),
            '-a',
            $annotation->path,
            '-g',
            $genome->basename(),
            '-t',
            $threads,
            '-f',
            $firstInputFile,
            '-o',
            $bamOutput,
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
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Second input file does not exist.',
                6 => 'Output file must be specified.',
                7 => 'Output directory is not writable.',
                8 => 'An error occurred during tophat2 execution.',
                9 => 'Unable to find output bam file.',
            ]
        );
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create TopHat output file');
        }
        $model->appendLog($output);
        $model->appendLog('Alignment completed.');

        return $bamOutput;
    }

    /**
     * Runs TopHat
     *
     * @param \App\Models\Job       $model
     * @param bool                  $paired
     * @param string                $firstInputFile
     * @param string|null           $secondInputFile
     * @param \App\Models\Reference $genome
     * @param int                   $threads
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runHisat(
        Job $model,
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        Reference $genome,
        int $threads = 1
    ): string {
        if ($genome->isAvailableFor('hisat')) {
            throw new ProcessingJobException('The selected reference has not been indexed for HISAT.');
        }
        $model->appendLog('Aligning with HISAT.');
        $bamOutput = $model->getJobTempFileAbsolute('hisat_output', '.bam');
        $command = [
            'bash',
            AbstractJob::scriptPath('hisat.sh'),
            '-g',
            $genome->basename(),
            '-t',
            $threads,
            '-f',
            $firstInputFile,
            '-o',
            $bamOutput,
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
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
                7 => 'An error occurred during HISAT 2 execution.',
                8 => 'Unable to find output sam file.',
            ]
        );
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create HISAT output file');
        }
        $model->appendLog($output);
        $model->appendLog('Alignment completed.');

        return $bamOutput;
    }

    /**
     * Runs Salmon
     *
     * @param \App\Models\Job       $model
     * @param bool                  $paired
     * @param string                $firstInputFile
     * @param string|null           $secondInputFile
     * @param string                $inputType
     * @param \App\Models\Reference $transcriptome
     * @param int                   $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runSalmon(
        Job $model,
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        string $inputType,
        Reference $transcriptome,
        int $threads = 1
    ): array {
        if (!$transcriptome->isAvailableFor('salmon')) {
            throw new ProcessingJobException('The specified reference sequence is not indexed for salmon.');
        }
        $model->appendLog('Quantifying with Salmon.');
        $salmonOutputRelative = $model->getJobTempFile('salmon_output', '_sa.txt');
        $salmonOutput = $model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = Storage::disk('public')->url($salmonOutputRelative);
        switch ($inputType) {
            case 'fastq':
                $command = [
                    'bash',
                    AbstractJob::scriptPath('salmon_counting.sh'),
                    '-i',
                    $transcriptome->basename(),
                    '-f',
                    $firstInputFile,
                    '-t',
                    $threads,
                    '-o',
                    $salmonOutput,
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
                        5 => 'Output file must be specified.',
                        6 => 'Output directory is not writable.',
                        7 => 'Indexed transcriptome does not exist.',
                        8 => 'Unable to find output file.',
                        9 => 'An error occurred during salmon quant execution.',
                    ]
                );
                break;
            case 'BAM':
                $output = AbstractJob::runCommand(
                    [
                        'bash',
                        AbstractJob::scriptPath('salmon_counting_bam.sh'),
                        '-r',
                        $transcriptome->path,
                        '-i',
                        $firstInputFile,
                        '-t',
                        $threads,
                        '-o',
                        $salmonOutput,
                    ],
                    $model->getAbsoluteJobDirectory(),
                    null,
                    null,
                    [
                        3 => 'Input file does not exist.',
                        4 => 'FASTA transcriptome file does not exist.',
                        5 => 'Output directory must be specified.',
                        6 => 'Output directory is not writable.',
                        7 => 'Unable to find output file.',
                        8 => 'An error occurred during salmon quant execution.',
                    ]
                );
                break;
            default:
                throw new ProcessingJobException('Unsupported input type');
        }
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $model->appendLog($output);
        $model->appendLog('Transcripts quantification completed.');

        return [$salmonOutputRelative, $salmonOutputUrl];
    }

}
