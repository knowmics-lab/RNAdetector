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

trait UseCountingTrait
{

    use ConvertsBamToFastqTrait;

    /**
     * Runs HTseq-count
     *
     * @param \App\Models\Job        $model
     * @param string                 $countingInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runHTSEQ(Job $model, string $countingInputFile, Annotation $annotation, int $threads = 1): array
    {
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $htseqOutputRelative = $model->getJobTempFile('htseq_output', '_ht.txt');
        $htseqOutput = $model->absoluteJobPath($htseqOutputRelative);
        $htseqOutputUrl = \Storage::disk('public')->url($htseqOutputRelative);
        $output = AbstractJob::runCommand(
            [
                'bash',
                AbstractJob::scriptPath('htseqcount.bash'),
                '-a',
                $annotation->path,
                '-b',
                $countingInputFile,
                '-t',
                $threads,
                '-o',
                $htseqOutput,
            ],
            $model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
            ]
        );
        if (!file_exists($htseqOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }

        return [$htseqOutputRelative, $htseqOutputUrl, $output];
    }

    /**
     * Runs FeatureCount
     *
     * @param \App\Models\Job        $model
     * @param string                 $countingInputFile
     * @param \App\Models\Annotation $annotation
     * @param int                    $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runFeatureCount(Job $model, string $countingInputFile, Annotation $annotation, int $threads = 1): array
    {
        if (!$annotation->isGtf()) {
            throw new ProcessingJobException('The selected annotation must be in GTF format.');
        }
        $relativeOutput = $model->getJobTempFile('featurecount_output', '_fc.txt');
        $absoluteOutput = $model->absoluteJobPath($relativeOutput);
        $outputUrl = \Storage::disk('public')->url($relativeOutput);
        $output = AbstractJob::runCommand(
            [
                'bash',
                AbstractJob::scriptPath('htseqcount.bash'),
                '-a',
                $annotation->path,
                '-b',
                $countingInputFile,
                '-t',
                $threads,
                '-o',
                $absoluteOutput,
            ],
            $model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
            ]
        );
        if (!file_exists($absoluteOutput)) {
            throw new ProcessingJobException('Unable to create FeatureCount output file');
        }

        return [$relativeOutput, $outputUrl, $output];
    }

    /**
     * Run salmon for SmallRNA counting
     *
     * @param bool                  $paired
     * @param string                $topHatInputFile
     * @param \App\Models\Reference $transcriptome
     * @param int                   $threads
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function runSalmonCount(Job $model, bool $paired, string $topHatInputFile, Reference $transcriptome, int $threads = 1): array
    {
        if (!$transcriptome->isAvailableFor('salmon')) {
            throw new ProcessingJobException('The specified reference sequence is not indexed for salmon analysis.');
        }
        $model->appendLog('Converting TopHat output to FastQ');
        [$firstTempFastQ, $secondTempFastQ, $output] = self::convertBamToFastq($model, $paired, $topHatInputFile);
        $model->appendLog($output);
        $model->appendLog('TopHat output converted to FastQ');
        $model->appendLog('Running salmon');
        $salmonOutputRelative = $model->getJobTempFile('salmon_output', '_sa.txt');
        $salmonOutput = $model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = \Storage::disk('public')->url($salmonOutputRelative);
        $command = [
            'bash',
            AbstractJob::scriptPath('salmon_counting.sh'),
            '-i',
            $transcriptome->basename(),
            '-t',
            $threads,
            '-o',
            $salmonOutput,
            '-f',
            $firstTempFastQ,
        ];
        if ($paired) {
            $command[] = '-s';
            $command[] = $secondTempFastQ;
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
                7 => 'Indexed trascriptome does not exist.',
                8 => 'Unable to find output file.',
                9 => 'An error occurred during salmon quant execution.',
            ]
        );
        $model->appendLog($output);
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $model->appendLog('Count computation completed.');
        if (file_exists($firstTempFastQ)) {
            @unlink($firstTempFastQ);
        }
        if (file_exists($secondTempFastQ)) {
            @unlink($secondTempFastQ);
        }

        return [$salmonOutputRelative, $salmonOutputUrl];
    }

}
