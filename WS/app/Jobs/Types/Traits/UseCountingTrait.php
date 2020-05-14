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
        $htseqOutputRelative = $model->getJobFile('htseq_output_', '.txt');
        $htseqOutput = $model->absoluteJobPath($htseqOutputRelative);
        $htseqOutputUrl = \Storage::disk('public')->url($htseqOutputRelative);
        $harmonizedRelative = $model->getJobFile('htseq_harmonized_', '.txt');
        $harmonized = $model->absoluteJobPath($harmonizedRelative);
        $harmonizedUrl = \Storage::disk('public')->url($harmonizedRelative);
        $command = [
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
            '-h',
            $harmonized,
        ];
        AbstractJob::addMap($command, $annotation);
        AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            static function ($type, $buffer) use ($model) {
                $model->appendLog($buffer, false);
            },
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
                7 => 'Unable to find output file.',
                8 => 'Error running htseq-count.',
                9 => 'Unable to harmonize output file.',
            ]
        );
        if (!file_exists($htseqOutput)) {
            throw new ProcessingJobException('Unable to create HTseq-count output file');
        }
        if (!file_exists($harmonized)) {
            throw new ProcessingJobException('Unable to create harmonized output file');
        }
        // $model->appendLog($output);
        $model->appendLog('Counting completed.');

        return [$htseqOutputRelative, $htseqOutputUrl, $harmonizedRelative, $harmonizedUrl];
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
        $relativeOutput = $model->getJobFile('featurecount_output_', '.txt');
        $absoluteOutput = $model->absoluteJobPath($relativeOutput);
        $outputUrl = \Storage::disk('public')->url($relativeOutput);
        $harmonizedRelative = $model->getJobFile('featurecount_harmonized_', '.txt');
        $harmonized = $model->absoluteJobPath($harmonizedRelative);
        $harmonizedUrl = \Storage::disk('public')->url($harmonizedRelative);
        $command = [
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
            '-h',
            $harmonized,
        ];
        AbstractJob::addMap($command, $annotation);
        AbstractJob::runCommand(
            $command,
            $model->getAbsoluteJobDirectory(),
            null,
            static function ($type, $buffer) use ($model) {
                $model->appendLog($buffer, false);
            },
            [
                3 => 'Annotation file does not exist.',
                4 => 'Input file does not exist.',
                5 => 'Output file must be specified.',
                6 => 'Output directory is not writable.',
                7 => 'Unable to find output file.',
                8 => 'An error occurred during featureCounts execution.',
                9 => 'Unable to harmonize output file.',
            ]
        );
        if (!file_exists($absoluteOutput)) {
            throw new ProcessingJobException('Unable to create FeatureCount output file');
        }
        if (!file_exists($harmonized)) {
            throw new ProcessingJobException('Unable to create harmonized output file');
        }
        $model->appendLog('Counting completed.');


        return [$relativeOutput, $outputUrl, $harmonizedRelative, $harmonizedUrl];
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
        $salmonOutputRelative = $model->getJobFile('salmon_output_', '.txt');
        $salmonOutput = $model->absoluteJobPath($salmonOutputRelative);
        $salmonOutputUrl = \Storage::disk('public')->url($salmonOutputRelative);
        $harmonizedGeneRelative = $model->getJobFile('salmon_harmonized_genes_', '.txt');
        $harmonizedGene = $model->absoluteJobPath($harmonizedGeneRelative);
        $harmonizedGeneUrl = \Storage::disk('public')->url($harmonizedGeneRelative);
        $harmonizedTxRelative = $model->getJobFile('salmon_harmonized_transcripts_', '.txt');
        $harmonizedTx = $model->absoluteJobPath($harmonizedTxRelative);
        $harmonizedTxUrl = \Storage::disk('public')->url($harmonizedTxRelative);
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
            '-h',
            $harmonizedGene,
            '-n',
            $harmonizedTx,
        ];
        if ($paired) {
            $command[] = '-s';
            $command[] = $secondTempFastQ;
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
                7  => 'Indexed trascriptome does not exist.',
                8  => 'Unable to find output file.',
                9  => 'An error occurred during salmon quant execution.',
                10 => 'Harmonized transcripts output file is required.',
                11 => 'Unable to harmonize output file.',
            ]
        );
        // $model->appendLog($output);
        if (!file_exists($salmonOutput)) {
            throw new ProcessingJobException('Unable to create Salmon output file');
        }
        if (!file_exists($harmonizedGene) || !file_exists($harmonizedTx)) {
            throw new ProcessingJobException('Unable to create harmonized output files');
        }
        $model->appendLog('Counts computation completed.');
        if (file_exists($firstTempFastQ)) {
            @unlink($firstTempFastQ);
        }
        if (file_exists($secondTempFastQ)) {
            @unlink($secondTempFastQ);
        }

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
