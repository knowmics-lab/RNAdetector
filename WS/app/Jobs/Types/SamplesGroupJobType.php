<?php
/**
 * RNADetector Web Service
 *
 * @author A. La Ferlita, Ph.D. Student <alessandrolf90 at hotmail dot it>
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;

use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\ConvertsSamToBamTrait;
use App\Jobs\Types\Traits\HasCommonParameters;
use App\Jobs\Types\Traits\RunTrimGaloreTrait;
use App\Jobs\Types\Traits\UseAlignmentTrait;
use App\Jobs\Types\Traits\UseCountingTrait;
use App\Models\Job;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;
use ZipArchive;

class SamplesGroupJobType extends AbstractJob
{
    use HasCommonParameters, ConvertsSamToBamTrait, RunTrimGaloreTrait, UseAlignmentTrait, UseCountingTrait;

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'jobs'        => 'A list of analysis job of the same type',
            'description' => 'An optional tsv file containing samples descriptions',
            'de_novo'     => 'An optional boolean to ignore all pre-built samples descriptions',
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
            'jobs'                      => 'A list of the valid job identifiers',
            'codes'                     => 'A list of the valid job sample codes',
            'description'               => 'An optional path/url of the samples descriptions file (filtered by valid codes)',
            'metadata'                  => 'An optional list of available metadata for this sample group (columns of the description file)',
            'outputFile'                => 'The raw output files',
            'harmonizedFile'            => 'The harmonized output file',
            'harmonizedTranscriptsFile' => 'An optional transcripts harmonized output file',
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
            'jobs'        => ['required', 'array', Rule::exists('jobs', 'id')],
            'description' => ['filled', 'string'],
            'de_novo'     => ['filled', 'boolean'],
        ];
    }

    /**
     * @inheritDoc
     */
    public function isInputValid(): bool
    {
        $description = $this->getParameter('description', null);
        $disk = Storage::disk('public');
        $dir = $this->getJobDirectory() . '/';

        return !(!empty($description) && !$disk->exists($dir . $description));
    }

    /**
     * @param \App\Models\Job[] $jobs
     *
     * @return array
     */
    private function processSamplesGroups(array $jobs): array
    {
        $samplesGroups = [];
        $others = [];
        foreach ($jobs as $job) {
            if ($job->job_type === 'samples_group_job_type') {
                $samplesGroups[] = $job;
            } else {
                $others[$job->id] = $job;
            }
        }
        $samplesGroupsJobs = [];
        $descriptionFiles = [];
        foreach ($samplesGroups as $job) {
            $containedJobs = $job->getOutput('jobs');
            $descriptionFiles[] = $job->absoluteJobPath($job->getOutput('description')['path']);
            if ($jobs !== null && is_array($jobs)) {
                foreach ($containedJobs as $cjId) {
                    if (isset($others[$cjId])) {
                        $samplesGroupsJobs[$cjId] = $others[$cjId];
                        unset($others[$cjId]);
                    } else {
                        $samplesGroupsJobs[$cjId] = Job::whereId($cjId)->first();
                    }
                }
            }
        }
        $finalJobsArray = array_merge(array_values($others), array_filter(array_values($samplesGroupsJobs)));

        return [$finalJobsArray, $descriptionFiles];
    }

    /**
     * Find all valid jobs
     *
     * @param array $jobs
     *
     * @return \App\Models\Job[]
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function processValidJobs(array $jobs): array
    {
        /** @var Job[] $models */
        $models = array_filter(
            array_map(
                static function ($jobId) {
                    return Job::whereId($jobId)->first();
                },
                $jobs
            )
        );
        if (!count($models)) {
            throw new ProcessingJobException('No jobs have been specified');
        }

        return $models;
    }

    /**
     * Checks if all job types are the same
     *
     * @param \App\Models\Job[] $models
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkJobsTypes(array $models): string
    {
        $firstType = $models[0]->job_type;
        foreach ($models as $model) {
            if ($model->job_type !== $firstType) {
                throw new ProcessingJobException('All jobs must be of the same type');
            }
        }

        return $firstType;
    }

    /**
     * Pull a property from an array of models
     *
     * @param \App\Models\Job[] $models
     * @param string            $property
     *
     * @return array
     */
    private function pullProperty(array $models, string $property): array
    {
        return array_map(
            static function ($job) use ($property) {
                return $job->$property;
            },
            $models
        );
    }

    /**
     * Prepare the meta array using the output of the make_descriptions.R script
     *
     * @param string $descriptor
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function prepareMetaArray(string $descriptor): array
    {
        $fp = @fopen($descriptor, 'rb');
        if (!$fp) {
            throw new ProcessingJobException('Unable to open descriptor file');
        }
        $meta = [];
        while (($data = fgetcsv($fp, 0, "\t")) !== false) {
            $meta[] = [
                'field'   => $data[0],
                'type'    => $data[1],
                'content' => explode(';', $data[2]),
            ];
        }
        @fclose($fp);

        return $meta;
    }

    /**
     * Make the description file using the make_descriptions.R script
     *
     * @param string[] $validCodes
     * @param array    $preBuiltDescriptions
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function makeDescriptionsFile(array $validCodes, array $preBuiltDescriptions): array
    {
        $inputDescription = $this->model->getParameter('description', null);
        $deNovo = $this->getParameter('de_novo', false);
        if ($deNovo) {
            $preBuiltDescriptions = [];
        }
        $dir = $this->getJobDirectory() . '/';
        if (!empty($inputDescription) && Storage::disk('public')->exists($dir . $inputDescription)) {
            $preBuiltDescriptions[] = realpath($this->getAbsoluteJobDirectory() . '/' . $inputDescription);
        }
        $descriptionsListFile = $this->getJobFileAbsolute('descriptions_list_', '.txt');
        @chmod($descriptionsListFile, 0777);
        @file_put_contents($descriptionsListFile, implode(PHP_EOL, $preBuiltDescriptions) . PHP_EOL);
        $samplesListFile = $this->getJobFileAbsolute('samples_list_', '.txt');
        @chmod($samplesListFile, 0777);
        @file_put_contents($samplesListFile, implode(PHP_EOL, $validCodes) . PHP_EOL);
        $descriptionRelative = $this->getJobFile('description_', '.tsv');
        $descriptionFile = $this->absoluteJobPath($descriptionRelative);
        $descriptionUrl = Storage::disk('public')->url($descriptionRelative);
        $descriptorFile = $this->getJobFileAbsolute('description_descriptor_', '.tsv');
        self::runCommand(
            [
                'Rscript',
                self::scriptPath('make_descriptions.R'),
                '-d',
                $descriptionsListFile,
                '-g',
                $this->model->sample_code,
                '-s',
                $samplesListFile,
                '-o',
                $descriptionFile,
                '-p',
                $descriptorFile,
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            }
        );
        if (!file_exists($descriptionFile)) {
            throw new ProcessingJobException('Unable to create samples description file');
        }
        if (!file_exists($descriptorFile)) {
            throw new ProcessingJobException('Unable to create samples descriptor file');
        }
        @chmod($descriptionFile, 0777);
        $metas = $this->prepareMetaArray($descriptorFile);
        @unlink($descriptorFile);
        @unlink($descriptionsListFile);
        @unlink($samplesListFile);

        return [$descriptionRelative, $descriptionUrl, $metas];
    }

    /**
     * Checks if all jobs have completed (COMPLETED or FAILED state)
     *
     * @param \App\Models\Job[] $jobs
     *
     * @return bool
     */
    private function checksForCompletion(array $jobs): bool
    {
        $completed = true;
        foreach ($jobs as $job) {
            $completed = $completed && $job->hasCompleted();
        }

        return $completed;
    }

    /**
     * Checks if all jobs have completed (COMPLETED or FAILED state)
     * and returns only completed jobs.
     *
     * @param \App\Models\Job[] $jobs
     *
     * @return \App\Models\Job[]
     */
    private function waitForCompletion(array $jobs): array
    {
        $this->log('Waiting for grouped jobs to complete.', false);
        while (!$this->checksForCompletion($jobs)) {
            sleep(600); // Wait for 10 minutes
            $this->log('.', false);
            // Refresh all jobs
            foreach ($jobs as $job) {
                $job->refresh();
            }
        }
        $this->log('');

        return array_filter(
            $jobs,
            static function (Job $job) {
                return $job->status === Job::COMPLETED;
            }
        );
    }

    /**
     * Make sample compose file
     *
     * @param \App\Models\Job[] $jobs
     *
     * @return array
     */
    private function makeSampleComposeFile(array $jobs): array
    {

        $sampleComposeData = array_filter(
            array_map(
                static function (Job $job) {
                    $res = Factory::sampleGroupFunctions($job);
                    if ($res === null) {
                        return null;
                    }

                    return array_map(
                        static function (callable $fn) use ($job) {
                            return $fn($job) ?? 'NA';
                        },
                        $res
                    );
                },
                $jobs
            )
        );
        $content = implode(
                PHP_EOL,
                array_map(
                    static function ($data) {
                        return implode("\t", $data);
                    },
                    $sampleComposeData
                )
            ) . PHP_EOL;
        $sampleComposeFile = $this->getJobFileAbsolute('job_compose_', '.txt');
        @file_put_contents($sampleComposeFile, $content);
        @chmod($sampleComposeFile, 0777);
        $hasTranscripts = true;
        foreach ($sampleComposeData as $datum) {
            if ($datum[3] === 'NA') {
                $hasTranscripts = false;
                break;
            }
        }

        return [$sampleComposeFile, $sampleComposeData, $hasTranscripts];
    }

    /**
     * Make raw output zip file
     *
     * @param array $sampleComposeContent
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function createRawZip(array $sampleComposeContent): array
    {
        $outputRelative = $this->getJobFile('raw_output_', '.zip');
        $output = $this->absoluteJobPath($outputRelative);
        $outputUrl = \Storage::disk('public')->url($outputRelative);
        $zip = new ZipArchive();
        if ($zip->open($output, ZipArchive::CREATE | ZipArchive::OVERWRITE) === true) {
            foreach ($sampleComposeContent as $content) {
                if ($content[2] !== null) {
                    $zip->addFile($content[2], basename($content[2]));
                }
            }
            $zip->close();
        } else {
            throw new ProcessingJobException('Unable to create raw output zip file.');
        }
        @chmod($output, 0777);

        return [$outputRelative, $outputUrl];
    }

    /**
     * @param string $sampleComposeFile
     * @param bool   $ciri
     * @param bool   $transcripts
     *
     * @return array
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function composeOutputFile(string $sampleComposeFile, bool $ciri, bool $transcripts): array
    {
        $outputRelative = $this->getJobFile('harmonized_output_', '.txt');
        $output = $this->absoluteJobPath($outputRelative);
        $outputUrl = \Storage::disk('public')->url($outputRelative);
        $txRelative = null;
        $txFile = null;
        $txUrl = null;
        $command = [
            'Rscript',
            self::scriptPath('compose.R'),
            '-i',
            $sampleComposeFile,
            '-o',
            $output,
        ];
        if ($ciri) {
            $command[] = '-c';
        } elseif ($transcripts) {
            $txRelative = $this->getJobFile('harmonized_transcripts_output_', '.txt');
            $txFile = $this->absoluteJobPath($outputRelative);
            $txUrl = \Storage::disk('public')->url($txRelative);
            $command[] = '-t';
            $command[] = '-s';
            $command[] = $txFile;
        }
        self::runCommand(
            $command,
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            }
        );
        @chmod($output, 0777);
        if (!$ciri && $transcripts) {
            @chmod($txFile, 0777);
        }
        if (!file_exists($output)) {
            throw new ProcessingJobException('Unable to create harmonized output file');
        }
        if (!$ciri && $transcripts && !file_exists($txFile)) {
            throw new ProcessingJobException('Unable to create harmonized transcripts output file');
        }

        return [$outputRelative, $outputUrl, $txRelative, $txUrl];
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
        $jobs = $this->getParameter('jobs', []);
        $models = $this->processValidJobs($jobs);
        $models = $this->waitForCompletion($models);
        $this->log('All jobs have been completed.');
        [$models, $preBuiltDescriptions] = $this->processSamplesGroups($models);
        $type = $this->checkJobsTypes($models);
        /** @var int[] $validJobs */
        $validJobs = $this->pullProperty($models, 'id');
        /** @var string[] $validCodes */
        $validCodes = $this->pullProperty($models, 'sample_code');
        $this->log('Processing description file.');
        [$descriptionRelative, $descriptionUrl, $metadata] = $this->makeDescriptionsFile(
            $validCodes,
            $preBuiltDescriptions
        );
        [$sampleComposeFile, $sampleComposeContent, $hasTranscripts] = $this->makeSampleComposeFile($models);
        $this->log('Creating raw output zip file.');
        [$rawPath, $rawUrl] = $this->createRawZip($sampleComposeContent);
        $this->log('Creating harmonized output file.');
        $isCiri = $type === 'circ_rna_job_type';
        [$harmonizedPath, $harmonizedUrl, $txPath, $txUrl] = $this->composeOutputFile(
            $sampleComposeFile,
            $isCiri,
            !$isCiri && $hasTranscripts
        );
        $output = [
            'type'           => self::OUT_TYPE_ANALYSIS_HARMONIZED_DESCRIPTION,
            'jobs'           => $validJobs,
            'codes'          => $validCodes,
            'description'    => ['path' => $descriptionRelative, 'url' => $descriptionUrl],
            'metadata'       => $metadata,
            'outputFile'     => ['path' => $rawPath, 'url' => $rawUrl],
            'harmonizedFile' => ['path' => $harmonizedPath, 'url' => $harmonizedUrl],
        ];
        if ($txPath !== null && $txUrl !== null) {
            $output['type'] = self::OUT_TYPE_ANALYSIS_HARMONIZED_TRANSCRIPTS_DESCRIPTION;
            $output['harmonizedTranscriptsFile'] = ['path' => $txPath, 'url' => $txUrl];
        }
        $this->model->setOutput($output);
        $this->model->save();
        $this->log('Sample group has been created.');
    }


    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Groups several jobs of the same type together (ideally a group of samples)';
    }
}
