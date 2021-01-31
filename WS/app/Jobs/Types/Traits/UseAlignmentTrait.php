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
     * Runs STAR
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
    private function runSTAR(
        Job $model,
        bool $paired,
        string $firstInputFile,
        ?string $secondInputFile,
        Reference $genome,
        Annotation $annotation,
        int $threads = 1
    ): string {
        if (!$genome->isAvailableFor('star')) {
            throw new ProcessingJobException('The selected reference has not been indexed for star.');
        }
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $model->appendLog('Aligning reads using STAR.');
        $bamOutput = $model->getJobFileAbsolute('star_output_', '.bam');
        $command = [
            'bash',
            AbstractJob::scriptPath('star.sh'),
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
        AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            static function ($type, $buffer) use ($model) {
                $model->appendLog($buffer, false);
            },
            [
                3   => 'Annotation file does not exist.',
                4   => 'Input file does not exist.',
                5   => 'Second input file does not exist.',
                6   => 'Output file must be specified.',
                7   => 'Output directory is not writable.',
                8   => 'Unable to create temporary directory.',
                9   => 'An error occurred during STAR alignment.',
                10  => 'Unable to find STAR output file.',
                11  => 'Unable to move output file.',
                12  => 'Unable to find output file.',
                101 => 'Unsorted BAM file does not exist.',
                102 => 'Unable to sort BAM file',
                103 => 'BAM file does not exist.',
                104 => 'Unable to write BAM index file.',
                105 => 'Unable to compute BAM coverage.',
            ]
        );
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create STAR output file');
        }
        $model->appendLog('Alignment completed.');

        return $bamOutput;
    }

    /**
     * Runs HISAT
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
        if (!$genome->isAvailableFor('hisat')) {
            throw new ProcessingJobException('The selected reference has not been indexed for HISAT.');
        }
        $model->appendLog('Aligning with HISAT.');
        $bamOutput = $model->getJobFileAbsolute('hisat_output_', '.bam');
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
        AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            static function ($type, $buffer) use ($model) {
                $model->appendLog($buffer, false);
            },
            [
                3   => 'Input file does not exist.',
                4   => 'Second input file does not exist.',
                5   => 'Output file must be specified.',
                6   => 'Output directory is not writable.',
                7   => 'An error occurred during HISAT 2 execution.',
                8   => 'Unable to find output sam file.',
                101 => 'Unsorted BAM file does not exist.',
                102 => 'Unable to sort BAM file',
                103 => 'BAM file does not exist.',
                104 => 'Unable to write BAM index file.',
                105 => 'Unable to compute BAM coverage.',
            ]
        );
        if (!file_exists($bamOutput)) {
            throw new ProcessingJobException('Unable to create HISAT output file');
        }
        // $model->appendLog($output);
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
        $salmonOutputRelative = $model->getJobFile('salmon_output_', '.txt');
        $salmonOutput = $model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = Storage::disk('public')->url($salmonOutputRelative);
        $harmonizedGeneRelative = $model->getJobFile('salmon_harmonized_genes_', '.txt');
        $harmonizedGene = $model->absoluteJobPath($harmonizedGeneRelative);
        $harmonizedGeneUrl = \Storage::disk('public')->url($harmonizedGeneRelative);
        $harmonizedTxRelative = $model->getJobFile('salmon_harmonized_transcripts_', '.txt');
        $harmonizedTx = $model->absoluteJobPath($harmonizedTxRelative);
        $harmonizedTxUrl = \Storage::disk('public')->url($harmonizedTxRelative);
        switch ($inputType) {
            case AbstractJob::FASTQ:
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
                    '-h',
                    $harmonizedGene,
                    '-n',
                    $harmonizedTx,
                ];
                if ($paired) {
                    $command[] = '-s';
                    $command[] = $secondInputFile;
                }
                AbstractJob::addMap($command, $transcriptome);
                AbstractJob::runCommand(
                    $command,
                    $model->getAbsoluteJobDirectory(),
                    null,
                    static function ($type, $buffer) use ($model) {
                        $model->appendLog($buffer, false);
                    },
                    [
                        3  => 'Input file does not exist.',
                        4  => 'Second input file does not exist.',
                        5  => 'Output file must be specified.',
                        6  => 'Output directory is not writable.',
                        7  => 'Indexed transcriptome does not exist.',
                        8  => 'Unable to find output file.',
                        9  => 'An error occurred during salmon quant execution.',
                        10 => 'Harmonized transcripts output file is required.',
                        11 => 'Unable to harmonize output file.',
                    ]
                );
                break;
            case AbstractJob::BAM:
                $command = [
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
                    '-h',
                    $harmonizedGene,
                    '-n',
                    $harmonizedTx,
                ];
                AbstractJob::addMap($command, $transcriptome);
                AbstractJob::runCommand(
                    $command,
                    $model->getAbsoluteJobDirectory(),
                    null,
                    static function ($type, $buffer) use ($model) {
                        $model->appendLog($buffer, false);
                    },
                    [
                        3  => 'Input file does not exist.',
                        4  => 'FASTA transcriptome file does not exist.',
                        5  => 'Output directory must be specified.',
                        6  => 'Output directory is not writable.',
                        7  => 'Unable to find output file.',
                        8  => 'An error occurred during salmon quant execution.',
                        9  => 'Harmonized transcripts output file is required.',
                        10 => 'Unable to harmonize output file.',
                    ]
                );
                break;
            default:
                throw new ProcessingJobException('Unsupported input type');
        }
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        if (!file_exists($harmonizedGene) || !file_exists($harmonizedTx)) {
            throw new ProcessingJobException('Unable to create harmonized output files');
        }
        $model->appendLog('Transcripts quantification completed.');

        return [
            $salmonOutputRelative,
            $salmonOutputUrl,
            $harmonizedGeneRelative,
            $harmonizedGeneUrl,
            $harmonizedTxRelative,
            $harmonizedTxUrl,
        ];
    }

}
