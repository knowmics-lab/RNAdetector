<?php

/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\HasCommonParameters;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class DiffExprAnalysisJobType extends AbstractJob
{
    use HasCommonParameters;

    private const  LIMMA = "limma";
    private const  DESEQ = "deseq";
    private const  EDGER = "edger";
    private const VALID_DEA_METHODS = [self::LIMMA, self::DESEQ, self::EDGER];

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return array_merge(
            [
                'formatted_input'         => 'Directory with formatted input files',
                'sample_info'             => 'File with samples information',
                'DiffExprAnalysisEnable'  => 'Boolean if users want to perform differential expression analysis or not',
                'DEA_tool'                => 'Differentially expressed analysis tools: DESeq, edgeR, LIMMA',
            ]
        );

    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [
            'outputFile' => 'Normalized read counts or differentially expressed genes',
        ];
    }

    /**
     * Returns an array containing rules for input validation.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return array
     */
    public static function validationSpec(Request $request): array
    {
        $parameters = (array)$request->get('parameters', []);

        return array_merge(
            [
                'formatted_input'         => ['filled', 'alpha_dash'],
                'sample_info'             => ['filled', 'alpha_dash'],
                'DiffExprAnalysisEnable'  => ['filled', 'boolean'],
                'DEA_tool'                => ['filled', Rule::in(self::VALID_DEA_METHODS)],
            ]
        );
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        $formatted_input = $this->model->getParameter('formatted_input');
        $sample_info =$this->model->getParameter('sample_info');
        $DiffExprAnalysisEnable = $this->model->getParameter('DiffExprAnalysisEnable', false);
        $DEA_tool = $this->model->getParameter('DEA_tool', self::LIMMA);
        if (!in_array($DEA_tool, self::VALID_DEA_METHODS, true)) {
            return false;
        }
        return true;
    }

    private function runNormalization (string $formatted_input): array
    {
        $normalizationOutputRelative = $this->model->getJobFile('normalized_output', '.txt');
        $normalizationOutput = $this->model->absoluteJobPath($normalizationOutputRelative);
        $normalizationOutputUrl = \Storage::disk('public')->url($normalizationOutputRelative);
        $output = self::runCommand(
            [
                'Rscript',
                self::scriptPath('normalization.R'),
                'path_input',
                $formatted_input,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                // exit codes
            ]
        );
        if (!file_exists($normalizationOutput)) {
            throw new ProcessingJobException('Unable to create normalization output file');
        }
        $this->log($output);

        return [$normalizationOutputRelative, $normalizationOutputUrl];
    }

    /**
     * Runs DESeq2
     *
     * @param string $formatted_input
     * @param string $sample_info
     * @return array
     * @throws ProcessingJobException
     */
    private function runDeseq (string $formatted_input, string $sample_info):array
    {
        $deseqOutputRelative = $this->model->getJobFile('deseq_output', '.txt');
        $deseqOutput = $this->model->absoluteJobPath($deseqOutputRelative);
        $deseqOutputUrl = \Storage::disk('public')->url($deseqOutputRelative);
        $output = self::runCommand(
            [
                'Rscript',
                self::scriptPath('deseq.R'),
                'path_input',
                $formatted_input,
                'sample_info_path',
                $sample_info,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                // exit codes
            ]
        );
        if (!file_exists($deseqOutput)) {
            throw new ProcessingJobException('Unable to create DESeq2 output file');
        }
        $this->log($output);

        return [$deseqOutputRelative, $deseqOutputUrl];
    }

    /**
     * Runs edgeR
     *
     * @param string $formatted_input
     * @param string $sample_info
     * @return array
     * @throws ProcessingJobException
     */

    private function runEdger (string $formatted_input, string $sample_info):array
    {
        $edgerOutputRelative = $this->model->getJobFile('edger_output', '.txt');
        $edgerOutput = $this->model->absoluteJobPath($edgerOutputRelative);
        $edgerOutputUrl = \Storage::disk('public')->url($edgerOutputRelative);
        $output = self::runCommand(
            [
                'Rscript',
                self::scriptPath('edgeR.R'),
                'path_input',
                $formatted_input,
                'sample_info_path',
                $sample_info,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                // exit codes
            ]
        );
        if (!file_exists($edgerOutput)) {
            throw new ProcessingJobException('Unable to create edgeR output file');
        }
        $this->log($output);

        return [$edgerOutputRelative, $edgerOutputUrl];
    }

    private function runLimma (string $formatted_input, string $sample_info):array
    {
        $limmaOutputRelative = $this->model->getJobFile('limma_output', '.txt');
        $limmaOutput = $this->model->absoluteJobPath($limmaOutputRelative);
        $limmaOutputUrl = \Storage::disk('public')->url($limmaOutputRelative);
        $output = self::runCommand(
            [
                'Rscript',
                self::scriptPath('limma.R'),
                'path_input',
                $formatted_input,
                'sample_info_path',
                $sample_info,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                // exit codes
            ]
        );
        if (!file_exists($limmaOutput)) {
            throw new ProcessingJobException('Unable to create LIMMA output file');
        }
        $this->log($output);

        return [$limmaOutputRelative, $limmaOutputUrl];
    }

    /**
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    public function handle(): void
    {
        $this->log('Starting diffential expression analysis.');
        $formatted_input = $this->model->getParameter('formatted_input');
        $sample_info = $this->model->getParameter('sample_info');
        $DiffExprAnalysisEnable = (bool)$this->model->getParameter('DiffExprAnalysisEnable', false);
        $DEA_tool = $this->model->getParameter('DEA_tool', self::LIMMA);
        if ($DiffExprAnalysisEnable){
            switch ($DEA_tool){
                case self::DESEQ:
                    $this->log('Starting differential expression analysis with DESeq 2');
                    [$outputFile, $outputUrl] = $this->runDeseq($formatted_input, $sample_info);
                    $this->log('Differential expression analysis completed');
                    break;
                case self::EDGER:
                    $this->log('Starting differential expression analysis with edgeR');
                    [$outputFile, $outputUrl] = $this->runEdger($formatted_input, $sample_info);
                    $this->log('Differential expression analysis completed');
                    break;
                case self::LIMMA:
                    $this->log('Starting differential expression analysis with LIMMA');
                    [$outputFile, $outputUrl] = $this->runLimma($formatted_input, $sample_info);
                    $this->log('Differential expression analysis completed');
                    break;
                default:
                    throw new ProcessingJobException("Invalid differential expression analysis tool");
            }
        } else {
            $this->log('Starting normalization of raw reads counts');
            [$outputFile, $outputUrl] = $this->runNormalization($formatted_input);
            $this->log('Normalization completed');
        }
        $this->model->setOutput(
            [
                'outputFile' => ['path' => $outputFile, 'url' => $outputUrl],
            ]
        );
        $this->log('Analysis completed.');
        $this->model->save();
    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Runs differential expression analysis';
    }

}
