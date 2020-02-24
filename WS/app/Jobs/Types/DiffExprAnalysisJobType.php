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

    private const  LIMMA                        = 'limma';
    private const  DESEQ                        = 'deseq';
    private const  EDGER                        = 'edger';
    private const  VALID_DEGS_METHODS           = [self::LIMMA, self::DESEQ, self::EDGER];
    private const  VALID_WHEN_APPLY_FILTERS     = ['prenorm', 'postnorm'];
    private const  VALID_NORM_METHODS           = [self::EDGER, self::DESEQ];
    private const  VALID_NORM_ARGS_METHOD       = ['TMM', 'TMMwsp', 'RLE', 'upperquartile', 'none'];
    private const  VALID_NORM_ARGS_LOCFUNC      = ['median', 'shorth'];
    private const  VALID_DESEQ_FITTYPE          = ['parametric', 'local', 'mean'];
    private const  VALID_EDGER_MAIN_METHOD      = ['classic', 'glm'];
    private const  VALID_EDGER_TREND            = ['movingave', 'loess', 'none'];
    private const  VALID_EDGER_TAG_METHOD       = ['grid', 'optimize'];
    private const  VALID_EDGER_GLM_METHOD       = ['CoxReid', 'Pearson', 'deviance'];
    private const  VALID_EDGER_GLM_TREND_METHOD = ['auto', 'bin.spline', 'bin.loess', 'power', 'spline'];
    private const  VALID_LIMMA_NORMALIZE_METHOD = ['none', 'scale', 'quantile', 'cyclicloess'];
    private const  VALID_FIG_FORMATS            = ['png', 'jpg', 'tiff', 'bmp', 'pdf'];
    private const  DEFAULT_FIG_FORMATS          = ['png', 'pdf'];
    private const  VALID_ADJUST_METHOD          = ['qvalue', 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'none'];
    private const  VALID_METAP_METHOD           = [
        'simes',
        'bonferroni',
        'minp',
        'maxp',
        'weight',
        'pandora',
        'dperm.min',
        'dperm.max',
        'dperm.weight',
        'fisher',
        'fperm',
        'whitlock',
        'none',
    ];

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        $degs = implode(', ', self::VALID_DEGS_METHODS);

        return [
            'source_sample_group' => 'The source sample group used to get samples data',
            'sample_type'         => 'If supported, the type of analysis: gene/transcript',
            'condition_variables' => 'Variables in the source datset used to determine samples conditions',
            'contrasts'           => 'A list of contrasts to determine differential expression: an array containing arrays with two elements (control and case)',
            'parameters'          => [
                'pcut'              => 'a p-value cutoff for exporting differentially genes (default: 0.05)',
                'log_offset'        => 'an offset to be added to values during logarithmic transformations in order to avoid Infinity (default: 1)',
                'when.apply.filter' => 'when gene filters will be applied relative to normalization: "prenorm" (default), "postnorm"',
                'norm'              => 'normalization algorithm to be applied: "edger" (default), "deseq"',
                'norm.args'         => [
                    'method'  => 'If norm is edger, the normalization method to be used: "TMM" (default),"TMMwsp","RLE","upperquartile","none"',
                    'locfunc' => 'If norm is deseq, a function to compute a location for a sample: "median" (default), "shorth" (more reliable for low counts)',
                ],
                'stats'             => 'an array of algorithms to use for DEGs analysis. Valid options are: ' . $degs . '.',
                'stats.args'        => [
                    self::DESEQ => [
                        'fitType' => 'The type of fitting of dispersions to the mean intensity: "parametric" (default), "local", "mean".',
                    ],
                    self::EDGER => [
                        'main.method'   => 'Type of edger analysis: "classic" (default), "glm".',
                        'rowsum.filter' => 'If main.method is "classic", genes with total count (across all samples) below this value will be filtered out before estimating the dispersion (default: 5).',
                        'trend'         => 'If main.method is "classic", method for estimating dispersion trend: "movingave" (default), "loess", "none".',
                        'tag.method'    => 'If main.method is "classic", method for maximizing the posterior likelihood: "grid" (default), "optimize".',
                        'glm.method'    => 'If main.method is "glm", method for estimating the dispersion: "CoxReid" (default), "Pearson" or "deviance".',
                        'trend.method'  => 'If main.method is "glm", method used to estimate the trended dispersions: "auto" (default), "bin.spline", "bin.loess", "power", or "spline".',
                    ],
                    self::LIMMA => [
                        'normalize.method' => 'The microarray-style normalization method to be applied to the logCPM values: "none" (default), "scale", "quantile", "cyclicloess".',
                    ],
                ],
                'filters'           => [
                    'length'     => [
                        'length' => 'Genes/transcripts are accepted for further analysis if they are above specified length in kb. Use null to disable',
                    ],
                    'avg.reads'  => [
                        'average.per.bp' => 'A gene is accepted for further analysis if it has more average reads than the quantile of the average count distribution per "average.per.bp" base pairs.',
                        'quantile'       => 'A gene is accepted for further analysis if it has more average reads than the "quantile" of the average count distribution per average.per.bp base pairs.',
                    ],
                    'expression' => [
                        'median'   => 'Genes below the median of the overall count distribution are not accepted (boolean value, default: TRUE)',
                        'mean'     => 'Genes below the mean of the overall count distribution are not accepted (boolean value, default: FALSE)',
                        'quantile' => 'Genes below the "quantile" quantile of the overall count distribution are not accepted (number or null, default: NULL)',
                        'known'    => 'A set of known not-expressed genes in the system under investigation are used to estimate an expression cutoff (array of gene names or null, default: NULL)',
                    ],
                    'presence'   => [
                        'frac'          => 'A gene is considered for statistical testing if "frac"% of samples have more than "min.count" reads',
                        'min.count'     => 'A gene is considered for statistical testing if "frac"% of samples have more than "min.count" reads',
                        'per.condition' => 'The check is on all samples or group by group (boolean value, default: FALSE)',
                    ],
                ],
                'adjust.method'     => 'A p-value adjustment method: "qvalue" (DEFAULT), "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".',
                'meta.p.method'     => 'A method used for results meta-analysis if more than one DEGs algorithms are chosen: "simes" (default), "bonferroni", "minp", "maxp", "weight", "pandora", "dperm.min", "dperm.max", "dperm.weight", "fisher", "fperm", "whitlock" or "none"',
                'fig.formats'       => 'An array of figure formats used to export results: "png" (default), "jpg", "tiff", "bmp", "pdf" (default).',
                'num.cores'         => 'Number of cores for parallel processing.',
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
        /* @TODO */
        $parameters = (array)$request->get('parameters', []);

        return array_merge(
            [
                'formatted_input'        => ['filled', 'alpha_dash'],
                'sample_info'            => ['filled', 'alpha_dash'],
                'DiffExprAnalysisEnable' => ['filled', 'boolean'],
                'DEA_tool'               => ['filled', Rule::in(self::VALID_DEGS_METHODS)],
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
        /* @TODO */
        $formatted_input = $this->model->getParameter('formatted_input');
        $sample_info = $this->model->getParameter('sample_info');
        $DiffExprAnalysisEnable = $this->model->getParameter('DiffExprAnalysisEnable', false);
        $DEA_tool = $this->model->getParameter('DEA_tool', self::LIMMA);
        if (!in_array($DEA_tool, self::VALID_DEGS_METHODS, true)) {
            return false;
        }

        return true;
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
        /* @TODO */
        /*$this->log('Starting diffential expression analysis.');
        $formatted_input = $this->model->getParameter('formatted_input');
        $sample_info = $this->model->getParameter('sample_info');
        $DiffExprAnalysisEnable = (bool)$this->model->getParameter('DiffExprAnalysisEnable', false);
        $DEA_tool = $this->model->getParameter('DEA_tool', self::LIMMA);
        if ($DiffExprAnalysisEnable) {
            switch ($DEA_tool) {
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
        $this->log('Analysis completed.');*/
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
