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

class DiffExprAnalysisJobType extends AbstractJob
{
    use HasCommonParameters;

    private const  LIMMA                        = 'limma';
    private const  DESEQ                        = 'deseq';
    private const  EDGER                        = 'edger';
    private const  VALID_SAMPLE_TYPE            = ['gene', 'transcript'];
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

    private const FILTER_DEFAULTS = [
        'length'     => [
            'length' => null,
        ],
        'avg_reads'  => [
            'average_per_bp' => 100,
            'quantile'       => 0.75,
        ],
        'expression' => [
            'median'   => true,
            'mean'     => false,
            'quantile' => null,
            'known'    => null,
        ],
        'presence'   => [
            'frac'          => 0.25,
            'min_count'     => 10,
            'per_condition' => false,
        ],
    ];

    private const DEFAULT_PARAMETERS = [
        'pcut'              => 0.05,
        'log_offset'        => 1,
        'when_apply_filter' => self::VALID_WHEN_APPLY_FILTERS[0],
        'norm'              => self::VALID_NORM_METHODS[0],
        'norm_args'         => [
            'method'  => self::VALID_NORM_ARGS_METHOD[0],
            'locfunc' => self::VALID_NORM_ARGS_LOCFUNC[0],
        ],
        'stats'             => [self::VALID_DEGS_METHODS[0]],
        'stats_args'        => [
            self::DESEQ => [
                'fitType' => self::VALID_DESEQ_FITTYPE[0],
            ],
            self::EDGER => [
                'main_method'   => self::VALID_EDGER_MAIN_METHOD[0],
                'rowsum_filter' => 5,
                'trend'         => self::VALID_EDGER_TREND[0],
                'tag_method'    => self::VALID_EDGER_TAG_METHOD[0],
                'glm_method'    => self::VALID_EDGER_GLM_METHOD[0],
                'trend_method'  => self::VALID_EDGER_GLM_TREND_METHOD[0],
            ],
            self::LIMMA => [
                'normalize_method' => self::VALID_LIMMA_NORMALIZE_METHOD[0],
            ],
        ],
        'filters'           => [
            'length'     => [
                'length' => null,
            ],
            'avg_reads'  => [
                'average_per_bp' => 100,
                'quantile'       => 0.75,
            ],
            'expression' => [
                'median'   => true,
                'mean'     => false,
                'quantile' => null,
                'known'    => null,
            ],
            'presence'   => null,
        ],
        'adjust_method'     => self::VALID_ADJUST_METHOD[0],
        'meta_p_method'     => self::VALID_METAP_METHOD[0],
        'fig_formats'       => self::DEFAULT_FIG_FORMATS,
        'num_cores'         => 1,
    ];

    private const PARAMETERS_KEY_CONVERSION = [
        'log_offset'        => 'log.offset',
        'when_apply_filter' => 'when.apply.filter',
        'norm_args'         => 'norm.args',
        'stats_args'        => 'stats.args',
        'main_method'       => 'main.method',
        'rowsum_filter'     => 'rowsum.filter',
        'tag_method'        => 'tag.method',
        'glm_method'        => 'glm.method',
        'trend_method'      => 'trend.method',
        'normalize_method'  => 'normalize.method',
        'avg_reads'         => 'avg.reads',
        'average_per_bp'    => 'average.per.bp',
        'adjust_method'     => 'adjust.method',
        'meta_p_method'     => 'meta.p.method',
        'fig_formats'       => 'fig.formats',
        'num_cores'         => 'num.cores',
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
                'when_apply_filter' => 'when gene filters will be applied relative to normalization: "prenorm" (default), "postnorm"',
                'norm'              => 'normalization algorithm to be applied: "edger" (default), "deseq"',
                'norm_args'         => [
                    'method'  => 'If norm is edger, the normalization method to be used: "TMM" (default),"TMMwsp","RLE","upperquartile","none"',
                    'locfunc' => 'If norm is deseq, a function to compute a location for a sample: "median" (default), "shorth" (more reliable for low counts)',
                ],
                'stats'             => 'an array of algorithms to use for DEGs analysis. Valid options are: ' . $degs . '.',
                'stats_args'        => [
                    self::DESEQ => [
                        'fitType' => 'The type of fitting of dispersions to the mean intensity: "parametric" (default), "local", "mean".',
                    ],
                    self::EDGER => [
                        'main_method'   => 'Type of edger analysis: "classic" (default), "glm".',
                        'rowsum_filter' => 'If main_method is "classic", genes with total count (across all samples) below this value will be filtered out before estimating the dispersion (default: 5).',
                        'trend'         => 'If main_method is "classic", method for estimating dispersion trend: "movingave" (default), "loess", "none".',
                        'tag_method'    => 'If main_method is "classic", method for maximizing the posterior likelihood: "grid" (default), "optimize".',
                        'glm_method'    => 'If main_method is "glm", method for estimating the dispersion: "CoxReid" (default), "Pearson" or "deviance".',
                        'trend_method'  => 'If main_method is "glm", method used to estimate the trended dispersions: "auto" (default), "bin.spline", "bin.loess", "power", or "spline".',
                    ],
                    self::LIMMA => [
                        'normalize_method' => 'The microarray-style normalization method to be applied to the logCPM values: "none" (default), "scale", "quantile", "cyclicloess".',
                    ],
                ],
                'filters'           => [
                    'length'     => [
                        'length' => 'Genes/transcripts are accepted for further analysis if they are above specified length in kb. Use null to disable',
                    ],
                    'avg_reads'  => [
                        'average_per_bp' => 'A gene is accepted for further analysis if it has more average reads than the quantile of the average count distribution per "average.per.bp" base pairs.',
                        'quantile'       => 'A gene is accepted for further analysis if it has more average reads than the "quantile" of the average count distribution per average.per.bp base pairs.',
                    ],
                    'expression' => [
                        'median'   => 'Genes below the median of the overall count distribution are not accepted (boolean value, default: TRUE)',
                        'mean'     => 'Genes below the mean of the overall count distribution are not accepted (boolean value, default: FALSE)',
                        'quantile' => 'Genes below the "quantile" quantile of the overall count distribution are not accepted (number or null, default: NULL)',
                        'known'    => 'A set of known not-expressed genes in the system under investigation are used to estimate an expression cutoff (array of gene names or null, default: NULL)',
                    ],
                    'presence'   => [
                        'frac'          => 'A gene is considered for statistical testing if "frac"% of samples have more than "min_count" reads',
                        'min_count'     => 'A gene is considered for statistical testing if "frac"% of samples have more than "min_count" reads',
                        'per_condition' => 'The check is on all samples or group by group (boolean value, default: FALSE)',
                    ],
                ],
                'adjust_method'     => 'A p-value adjustment method: "qvalue" (DEFAULT), "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none".',
                'meta_p_method'     => 'A method used for results meta-analysis if more than one DEGs algorithms are chosen: "simes" (default), "bonferroni", "minp", "maxp", "weight", "pandora", "dperm.min", "dperm.max", "dperm.weight", "fisher", "fperm", "whitlock" or "none"',
                'fig_formats'       => 'An array of figure formats used to export results: "png" (default), "jpg", "tiff", "bmp", "pdf" (default).',
                'num_cores'         => 'Number of cores for parallel processing.',
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
        $statDeseq = 'stats_args.' . self::DESEQ;
        $statEdger = 'stats_args.' . self::EDGER;
        $statLimma = 'stats_args.' . self::LIMMA;

        return [
            'source_sample_group'              => ['required', Rule::exists('jobs', 'id')],
            'sample_type'                      => ['required', Rule::in(self::VALID_SAMPLE_TYPE)],
            'condition_variables'              => ['required', 'array'],
            'contrasts'                        => ['required', 'array', 'min:1'],
            'contrasts.*.control'              => ['required', 'string'],
            'contrasts.*.case'                 => ['required', 'string'],
            'parameters'                       => ['filled', 'array'],
            'pcut'                             => ['filled', 'numeric', 'min:0', 'max:1'],
            'log_offset'                       => ['filled', 'numeric'],
            'when_apply_filter'                => ['filled', Rule::in(self::VALID_WHEN_APPLY_FILTERS)],
            'norm'                             => ['filled', Rule::in(self::VALID_NORM_METHODS)],
            'norm_args'                        => ['filled', 'array'],
            'norm_args.method'                 => ['filled', Rule::in(self::VALID_NORM_ARGS_METHOD)],
            'norm_args.locfunc'                => ['filled', Rule::in(self::VALID_NORM_ARGS_LOCFUNC)],
            'stats'                            => ['filled', 'array', Rule::in(self::VALID_DEGS_METHODS)],
            'stats_args'                       => ['filled', 'array'],
            $statDeseq                         => ['filled', 'array'],
            $statDeseq . '.fitType'            => ['filled', Rule::in(self::VALID_DESEQ_FITTYPE)],
            $statEdger                         => ['filled', 'array'],
            $statEdger . '.main_method'        => ['filled', Rule::in(self::VALID_EDGER_MAIN_METHOD)],
            $statEdger . '.rowsum_filter'      => ['filled', 'numeric'],
            $statEdger . '.trend'              => ['filled', Rule::in(self::VALID_EDGER_TREND)],
            $statEdger . '.tag_method'         => ['filled', Rule::in(self::VALID_EDGER_TAG_METHOD)],
            $statEdger . '.glm_method'         => ['filled', Rule::in(self::VALID_EDGER_GLM_METHOD)],
            $statEdger . '.trend_method'       => ['filled', Rule::in(self::VALID_EDGER_GLM_TREND_METHOD)],
            $statLimma                         => ['filled', 'array'],
            $statLimma . '.normalize_method'   => ['filled', Rule::in(self::VALID_LIMMA_NORMALIZE_METHOD)],
            'filters'                          => ['filled', 'array'],
            'filters.length'                   => ['filled', 'nullable', 'array'],
            'filters.length.length'            => ['filled', 'nullable', 'numeric'],
            'filters.avg_reads'                => ['filled', 'nullable', 'array'],
            'filters.avg_reads.average_per_bp' => ['filled', 'numeric'],
            'filters.avg_reads.quantile'       => ['filled', 'numeric', 'min:0', 'max:1'],
            'filters.expression'               => ['filled', 'nullable', 'array'],
            'filters.expression.median'        => ['filled', 'boolean'],
            'filters.expression.mean'          => ['filled', 'boolean'],
            'filters.expression.quantile'      => ['filled', 'nullable', 'numeric', 'min:0', 'max:1'],
            'filters.expression.known'         => ['filled', 'nullable', 'array'],
            'filters.presence'                 => ['filled', 'nullable', 'array'],
            'filters.presence.frac'            => ['filled', 'numeric', 'min:0', 'max:1'],
            'filters.presence.min_count'       => ['filled', 'numeric'],
            'filters.presence.per_condition'   => ['filled', 'boolean'],
            'adjust_method'                    => ['filled', Rule::in(self::VALID_ADJUST_METHOD)],
            'meta_p_method'                    => ['filled', Rule::in(self::VALID_METAP_METHOD)],
            'fig_formats'                      => ['filled', 'array', Rule::in(self::VALID_FIG_FORMATS)],
            'num_cores'                        => ['filled', 'integer'],
        ];
    }

    /**
     * @param \App\Models\Job $sourceSampleGroup
     *
     * @return array
     */
    private function prepareMetadata(Job $sourceSampleGroup): array
    {
        $contents = [];
        $groupOutput = $sourceSampleGroup->job_output;
        foreach ($groupOutput['metadata'] as $m) {
            if ($m['type'] === 'string') {
                $contents[$m['field']] = $m['content'];
            } /*else {
                // unsupported right now
            }*/
        }

        return $contents;
    }

    /**
     * @param array $conditionVariables
     * @param array $meta
     *
     * @return void
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkConditionVariables(array $conditionVariables, array $meta): void
    {
        $invalid = [];
        foreach ($conditionVariables as $var) {
            if (!isset($meta[$var])) {
                $invalid[] = $var;
            }
        }
        if (count($invalid) > 0) {
            throw new ProcessingJobException('Found invalid condition variables: ' . implode(', ', $invalid) . '.');
        }
    }

    /**
     * @param array  $variables
     * @param array  $meta
     * @param array  $result
     * @param string $prefix
     */
    private function recursiveConditionBuilder(array $variables, array $meta, array &$result, string $prefix = ''): void
    {
        $first = array_shift($variables);
        foreach ($meta[$first] as $value) {
            $elem = $prefix . $value;
            if (count($variables) > 0) {
                $this->recursiveConditionBuilder($variables, $meta, $result, $elem . '_');
            } else {
                $result[] = $elem;
            }
        }
    }

    /**
     * @param array $conditionVariables
     * @param array $meta
     *
     * @return array
     */
    private function prepareAvailableConditions(array $conditionVariables, array $meta): array
    {
        $result = [];
        $this->recursiveConditionBuilder($conditionVariables, $meta, $result);

        return $result;
    }

    /**
     * @param array $contrasts
     * @param array $conditions
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkAndPrepareContrasts(array $contrasts, array $conditions): array
    {
        $validContrasts = [];
        $invalid = [];
        foreach ($contrasts as $contrast) {
            $case = $contrast['case'];
            $control = $contrast['control'];
            if (in_array($case, $conditions, true) && in_array($control, $conditions, true)) {
                $validContrasts[] = [$control, $case];
            } else {
                $invalid[] = $contrast;
            }
        }
        if (count($validContrasts) <= 0) {
            throw new ProcessingJobException('No valid contrasts found');
        }
        if (count($invalid) > 0) {
            $invalidString = implode(
                ', ',
                array_map(
                    static function ($e) {
                        return $e['control'] . ' vs ' . $e['case'];
                    },
                    $invalid
                )
            );
            $this->log('Found invalid contrasts: ' . $invalidString . '.');
        }

        return $validContrasts;
    }

    /**
     * @param array $arr
     *
     * @return array
     */
    private static function mapKeys(array $arr): array
    {
        $res = [];
        foreach ($arr as $key => $value) {
            $mapKey = self::PARAMETERS_KEY_CONVERSION[$key] ?? $key;
            if (is_array($value)) {
                $value = self::mapKeys($value);
            }
            $res[$mapKey] = $value;
        }

        return $res;
    }

    /**
     * @param array $parameters
     *
     * @return array
     */
    private function prepareParameters(array $parameters): array
    {
        $finalParameters = [];
        foreach (self::DEFAULT_PARAMETERS as $par => $val) {
            $mapKey = self::PARAMETERS_KEY_CONVERSION[$par] ?? $par;
            switch ($par) {
                case 'norm_args':
                    $finalParameters[$mapKey] = self::mapKeys(array_merge($val, $parameters[$par] ?? []));
                    break;
                case 'stats_args':
                    $args = $parameters[$par] ?? [];
                    $stats = $parameters['stats'] ?? self::DEFAULT_PARAMETERS['stats'];
                    $finalParameters[$mapKey] = [];
                    foreach ($stats as $s) {
                        $finalParameters[$mapKey][$s] = self::mapKeys(array_merge($val[$s], $args[$s] ?? []));
                    }
                    break;
                case 'filters':
                    $args = $parameters[$par] ?? [];
                    $finalParameters[$mapKey] = [];
                    foreach (array_keys($val) as $f) {
                        $filterArgs = $args[$f] ?? $val[$f];
                        $finalParameters[$mapKey][$f] = ($filterArgs !== null) ? self::mapKeys(
                            array_merge(self::FILTER_DEFAULTS[$f], $filterArgs)
                        ) : null;
                    }
                    break;
                default:
                    $finalParameters[$mapKey] = $parameters[$par] ?? $val;
                    break;
            }
        }

        return $finalParameters;
    }

    /**
     * @param \App\Models\Job $sourceSampleGroup
     * @param array           $groupOutput
     * @param string          $sampleType
     * @param array           $conditionVariables
     * @param array           $validContrasts
     * @param array           $parameters
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function prepareConfigFile(
        Job $sourceSampleGroup,
        array $groupOutput,
        string $sampleType,
        array $conditionVariables,
        array $validContrasts,
        array $parameters
    ): array {
        $configFile = $this->getJobFileAbsolute('deg_config_', '.json');
        $degReportDirectory = $this->model->getJobFile('deg_report_');
        $degReport = $this->model->absoluteJobPath($degReportDirectory);
        $degReportUrl = \Storage::disk('public')->url($degReportDirectory);
        $dataField = $sampleType === self::VALID_SAMPLE_TYPE[1] ? 'harmonizedTranscriptsFile' : 'harmonizedFile';
        $jsonContent = [
            'description.file'     => $sourceSampleGroup->absoluteJobPath($groupOutput['description']['path']),
            'data.file'            => $sourceSampleGroup->absoluteJobPath($groupOutput[$dataField]['path']),
            'data.type'            => $sampleType,
            'conditions.variables' => $conditionVariables,
            'contrasts'            => $validContrasts,
            'output.directory'     => $degReport,
            'parameters'           => $parameters,
        ];
        file_put_contents($configFile, json_encode($jsonContent));
        if (!file_exists($configFile)) {
            throw new ProcessingJobException('Unable to create analysis config file.');
        }

        return [$configFile, $degReportDirectory, $degReport, $degReportUrl];
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
        $sourceSampleGroup = Job::whereId($this->getParameter('source_sample_group'))->firstOrFail();
        $sampleType = $this->getParameter('sample_type', self::VALID_SAMPLE_TYPE[0]);
        $conditionVariables = (array)$this->getParameter('condition_variables');
        $contrasts = (array)$this->getParameter('contrasts');
        $parameters = (array)$this->getParameter('parameters');
        if ($sourceSampleGroup->job_type !== 'samples_group_job_type') {
            throw new ProcessingJobException('Source sample group error: job type invalid.');
        }
        $groupOutput = $sourceSampleGroup->job_output;
        if ($sampleType === self::VALID_SAMPLE_TYPE[1] && !isset($groupOutput['harmonizedTranscriptsFile'])) {
            throw new ProcessingJobException('The selected sample group does not support transcripts differential expression analysis.');
        }
        $metadata = $this->prepareMetadata($sourceSampleGroup);
        $this->checkConditionVariables($conditionVariables, $metadata);
        $conditions = $this->prepareAvailableConditions($conditionVariables, $metadata);
        $validContrasts = $this->checkAndPrepareContrasts($contrasts, $conditions);
        $parameters = $this->prepareParameters($parameters);
        [$configFile, $degReportDirectory, $degReport, $degReportUrl] = $this->prepareConfigFile(
            $sourceSampleGroup,
            $groupOutput,
            $sampleType,
            $conditionVariables,
            $validContrasts,
            $parameters
        );
        AbstractJob::runCommand(
            [
                'Rscript',
                self::scriptPath('de_analysis.R'),
                '-c',
                $configFile,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log(trim($buffer));
            }
        );
        if (!file_exists($degReport) && !is_dir($degReport) && !file_exists($degReport . '/index.html')) {
            throw new ProcessingJobException('Unable to create output report.');
        }
        $this->log('DEGs Analysis completed.');
        $degReportZip = $this->model->getJobFile('deg_report_', '.zip');
        $degReportZipAbsolute = $this->model->absoluteJobPath($degReportZip);
        $degReportZipUrl = \Storage::disk('public')->url($degReportZip);
        $this->log('Building report archive.');
        if (!Utils::makeZipArchive($degReport, $degReportZipAbsolute)) {
            throw new ProcessingJobException('Unknown error during output archive creation.');
        }
        if (!file_exists($degReportZipAbsolute)) {
            throw new ProcessingJobException('Unable to create output archive.');
        }
        $this->log('Archive built.');
        $this->setOutput(
            [
                'outputFile' => ['path' => $degReportZip, 'url' => $degReportZipUrl],
                'reportFile' => ['path' => $degReportDirectory . '/index.html', 'url' => $degReportUrl . '/index.html'],
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
        return 'Runs differential expression analysis';
    }

}
