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
        $model->appendLog('Counting with HTseq-count.');
        $htseqOutputRelative = $model->getJobFile('htseq_output_', '_ht.txt');
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
            static function ($type, $buffer) use ($model) {
                $model->appendLog(trim($buffer));
            },
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
                7 => 'Unable to find output file.',
                8 => 'Error running htseq-count.',
            ]
        );
        if (!file_exists($htseqOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }
        // $model->appendLog($output);
        $model->appendLog('Counting completed.');

        return [$htseqOutputRelative, $htseqOutputUrl];
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
        $model->appendLog('Counting with Feature-Count.');;
        $relativeOutput = $model->getJobFile('featurecount_output_', '_fc.txt');
        $absoluteOutput = $model->absoluteJobPath($relativeOutput);
        $outputUrl = \Storage::disk('public')->url($relativeOutput);
        $output = AbstractJob::runCommand(
            [
                'bash',
                AbstractJob::scriptPath('featurecounts.bash'),
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
            static function ($type, $buffer) use ($model) {
                $model->appendLog(trim($buffer));
            },
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
                7 => 'Unable to find output file.',
                8 => 'An error occurred during featureCounts execution.',
            ]
        );
        if (!file_exists($absoluteOutput)) {
            throw new ProcessingJobException('Unable to create FeatureCount output file');
        }
        // $model->appendLog($output);
        $model->appendLog('Counting completed.');


        return [$relativeOutput, $outputUrl];
    }

    /**
     * Run salmon for counting
     *
     * @param \App\Models\Job       $model
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
        [$firstTempFastQ, $secondTempFastQ] = self::convertBamToFastq($model, $paired, $topHatInputFile, true);
        $model->appendLog('Quantifying with salmon');
        $salmonOutputRelative = $model->getJobFile('salmon_output_', '_sa.txt');
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
            static function ($type, $buffer) use ($model) {
                $model->appendLog(trim($buffer));
            },
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
        // $model->appendLog($output);
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        $model->appendLog('Counts computation completed.');
        if (file_exists($firstTempFastQ)) {
            @unlink($firstTempFastQ);
        }
        if (file_exists($secondTempFastQ)) {
            @unlink($secondTempFastQ);
        }

        return [$salmonOutputRelative, $salmonOutputUrl];
    }

}
