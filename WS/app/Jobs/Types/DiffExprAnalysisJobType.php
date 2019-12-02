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
use App\Models\Annotation;
use App\Models\Reference;
use App\Utils;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

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
                'DEA_tool'          => 'Differentially expressed analysis tools: DESeq, edgeR, LIMMA',
                'formatted_input'   => 'Directory with formatted input files',
                'sample_info'       => 'File with samples information',
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
            'outputFile' => 'Differentially expressed genes',
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
                'DEA_tool'          => ['filled', Rule::in(self::VALID_DEA_METHODS)],
                'formatted_input'   => ['filled', 'alpha_dash'],
                'sample_info'       => ['filled', 'alpha_dash'],
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
        $DEA_tool = $this->model->getParameter('DEA_tool', self::LIMMA);
        if (!in_array($DEA_tool, self::VALID_DEA_METHODS, true)) {
            return false;
        }

        return true;
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
        $deseqOutputRelative = $this->model->getJobTempFile('deseq_output', '.txt');
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
            []
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
        $edgerOutputRelative = $this->model->getJobTempFile('edger_output', '.txt');
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
            []
        );
        if (!file_exists($edgerOutput)) {
            throw new ProcessingJobException('Unable to create edgeR output file');
        }
        $this->log($output);

        return [$edgerOutputRelative, $edgerOutputUrl];
    }

    private function runLimma (string $formatted_input, string $sample_info):array
    {
        $limmaOutputRelative = $this->model->getJobTempFile('limma_output', '.txt');
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
            []
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
