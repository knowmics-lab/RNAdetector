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
use App\Models\Job;
use App\Utils;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;

class PathwayAnalysisJobType extends AbstractJob
{
    use HasCommonParameters;

    private const VALID_ORGANISMS = ['hsa', 'mmu', 'rno'];

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'degs_analysis' => 'The DEGs analysis from which data are gathered',
            'degs'          => [
                'p_cutoff'      => 'a p-value cutoff for exporting differentially genes (default: 0.05)',
                'p_use_fdr'     => 'boolean indicating whether p-value cutoff is applied to FDR p-values (default: TRUE)',
                'lfc_threshold' => 'minimum absolute Log-Fold-Change cutoff for exporting differentially genes (default: 0)',
            ],
            'pathways'      => [
                'organism'  => 'the organism to use for pathway analysis (default: hsa). Supported values are: hsa, mmu, rno.',
                'p_cutoff'  => 'a p-value cutoff for exporting significantly impacted pathways (default: 0.05)',
                'p_use_fdr' => 'boolean indicating whether p-value cutoff is applied to FDR p-values (default: TRUE)',
            ],
        ];

    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [
            'outputFile' => 'Archive containing report file.',
            'reportFile' => 'Path of the report index file.',
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
        return [
            'degs_analysis'      => ['required', Rule::exists('jobs', 'id')],
            'degs_analysis'      => ['required', Rule::exists('jobs', 'id')],
            'degs'               => ['filled', 'array'],
            'degs.p_cutoff'      => ['filled', 'numeric', 'min:0', 'max:1'],
            'degs.p_use_fdr'     => ['filled', 'boolean'],
            'degs.lfc_threshold' => ['filled', 'numeric'],
            'pathways'           => ['filled', 'array'],
            'pathways.organism'  => ['filled', Rule::in(self::VALID_ORGANISMS)],
            'pathways.p_cutoff'  => ['filled', 'numeric', 'min:0', 'max:1'],
            'pathways.p_use_fdr' => ['filled', 'boolean'],
        ];
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
        $this->log('Starting pathway analysis.');
        $degsAnalysis = Job::whereId($this->getParameter('degs_analysis'))->firstOrFail();
        $degsParameters = $this->getParameter('degs', []);
        $pathwayParameters = $this->getParameter('pathways', []);
        $degsParameters['p_cutoff'] = $degsParameters['p_cutoff'] ?? 0.05;
        $degsParameters['p_use_fdr'] = $degsParameters['p_use_fdr'] ?? true;
        $degsParameters['lfc_threshold'] = $degsParameters['lfc_threshold'] ?? 0;
        $pathwayParameters['p_cutoff'] = $pathwayParameters['p_cutoff'] ?? 0.05;
        $pathwayParameters['p_use_fdr'] = $pathwayParameters['p_use_fdr'] ?? true;
        $pathwayParameters['organism'] = $pathwayParameters['organism'] ?? self::VALID_ORGANISMS[0];
        if ($degsAnalysis->job_type !== 'diff_expr_analysis_job_type') {
            throw new ProcessingJobException('DEGs analysis error: job type invalid.');
        }
        $degsOutput = $degsAnalysis->job_output;
        $degReportFile = $degsOutput['reportFile']['path'] ?? null;
        if (!$degReportFile) {
            throw new ProcessingJobException('The selected DEGs analysis does not contain any result.');
        }
        $degReportFileAbsolute = $this->model->absoluteJobPath($degReportFile);
        if (!file_exists($degReportFileAbsolute)) {
            throw new ProcessingJobException('The selected DEGs analysis does not contain any result.');
        }
        $degReportDirectory = dirname($degReportFileAbsolute);
        $pathReportDirectory = $this->model->getJobFile('pathway_report_');
        $pathReport = $this->model->absoluteJobPath($pathReportDirectory);
        $pathReportUrl = \Storage::disk('public')->url($pathReportDirectory);
        $command = [
            'Rscript',
            self::scriptPath('pathway_analysis.R'),
            '-i',
            $degReportDirectory,
            '-o',
            $pathReport,
            '--degs-p',
            $degsParameters['p_cutoff'],
            '--degs-lfc',
            $degsParameters['lfc_threshold'],
            '--path-org',
            $pathwayParameters['organism'],
            '--path-p',
            $pathwayParameters['p_cutoff'],
        ];
        if (!$degsParameters['p_use_fdr']) {
            $command[] = '--degs-no-fdr';
        }
        if (!$pathwayParameters['p_use_fdr']) {
            $command[] = '--path-no-fdr';
        }
        AbstractJob::runCommand(
            $command,
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            }
        );
        if (!file_exists($pathReport) && !is_dir($pathReport) && !file_exists($pathReport . '/index.html')) {
            throw new ProcessingJobException('Unable to create output report.');
        }
        Utils::recursiveChmod($pathReport, 0777);
        $this->log('Pathway Analysis completed.');
        $pathReportZip = $this->model->getJobFile('pathway_report_', '.zip');
        $pathReportZipAbsolute = $this->model->absoluteJobPath($pathReportZip);
        $pathReportZipUrl = \Storage::disk('public')->url($pathReportZip);
        $this->log('Building report archive.');
        if (!Utils::makeZipArchive($pathReport, $pathReportZipAbsolute)) {
            throw new ProcessingJobException('Unknown error during output archive creation.');
        }
        @chmod($pathReportZipAbsolute, 0777);
        if (!file_exists($pathReportZipAbsolute)) {
            throw new ProcessingJobException('Unable to create output archive.');
        }
        $this->log('Archive built.');
        $this->setOutput(
            [
                'type'       => self::OUT_TYPE_ANALYSIS_REPORT,
                'outputFile' => ['path' => $pathReportZip, 'url' => $pathReportZipUrl],
                'reportFile' => ['path' => $pathReportDirectory . '/index.html', 'url' => $pathReportUrl . '/index.html'],
            ]
        );
        $this->model->save();
    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Runs pathway analysis';
    }

    /**
     * @inheritDoc
     */
    public static function displayName(): string
    {
        return 'Pathway Analysis';
    }
}
