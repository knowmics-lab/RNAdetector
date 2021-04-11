<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ResubmitException;
use App\Models\Job;
use App\Models\Job as JobModel;
use App\Utils;
use Illuminate\Contracts\Filesystem\Filesystem;
use Illuminate\Http\Request;
use Storage;
use Symfony\Component\Process\Exception\ProcessFailedException;

/**
 * Class AbstractJob
 *
 * @method string getJobTempFile(string $prefix = '', string $suffix = '')
 * @method string getJobTempFileAbsolute(string $prefix = '', string $suffix = '')
 * @method string getJobFile(string $prefix = '', string $suffix = '')
 * @method string getJobFileAbsolute(string $prefix = '', string $suffix = '')
 * @method string absoluteJobPath(string $path)
 * @method string getAbsoluteJobDirectory()
 * @method string getJobDirectory()
 * @method mixed getParameter($parameter = null, $default = null)
 * @method \App\Models\Job setOutput($parameter, $value = null)
 * @package App\Jobs\Types
 */
abstract class AbstractJob
{

    public const FASTQ             = 'fastq';
    public const BAM               = 'BAM';
    public const SAM               = 'SAM';
    public const VALID_INPUT_TYPES = [self::FASTQ, self::BAM, self::SAM];

    public const HTSEQ_COUNTS           = 'htseq';
    public const FEATURECOUNTS_COUNTS   = 'feature-counts';
    public const SALMON                 = 'salmon';
    public const VALID_COUNTING_METHODS = [self::HTSEQ_COUNTS, self::FEATURECOUNTS_COUNTS, self::SALMON];

    public const HISAT2                    = 'hisat2';
    public const STAR                      = 'star';
    public const VALID_ALIGN_QUANT_METHODS = [self::SALMON, self::HISAT2, self::STAR];

    public const OUT_TYPE_CONFIRMATION                                = 'confirmation';
    public const OUT_TYPE_ANALYSIS                                    = 'analysis';
    public const OUT_TYPE_ANALYSIS_HARMONIZED                         = 'analysis-harmonized';
    public const OUT_TYPE_ANALYSIS_HARMONIZED_TRANSCRIPTS             = 'analysis-harmonized-transcripts';
    public const OUT_TYPE_ANALYSIS_HARMONIZED_DESCRIPTION             = 'analysis-harmonized-description';
    public const OUT_TYPE_ANALYSIS_HARMONIZED_TRANSCRIPTS_DESCRIPTION = 'analysis-harmonized-transcripts-description';
    public const OUT_TYPE_ANALYSIS_REPORT                             = 'analysis-report';

    /**
     * @var int
     */
    protected $defaultNumberOfThreads = 1;

    /**
     * @var \App\Models\Job
     */
    protected $model;

    /**
     * @var \Illuminate\Contracts\Filesystem\Filesystem
     */
    protected $publicFolder;

    /**
     * AbstractJob constructor.
     *
     * @param  \App\Models\Job  $model
     */
    public function __construct(JobModel $model)
    {
        $this->model = $model;
        $this->publicFolder = Storage::disk('public');
    }


    /**
     * Get the model containing all the data for this job
     *
     * @return \App\Models\Job
     */
    public function getModel(): JobModel
    {
        return $this->model;
    }

    /**
     * Set the model containing all the data for this job
     *
     * @param  \App\Models\Job  $model
     *
     * @return $this
     */
    public function setModel(JobModel $model): self
    {
        $this->model = $model;

        return $this;
    }

    /**
     * Returns the name that will be displayed for this job type
     *
     * @return string
     */
    abstract public static function displayName(): string;

    /**
     * Returns a description for this job
     *
     * @return string
     */
    abstract public static function description(): string;

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    abstract public static function parametersSpec(): array;

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    abstract public static function outputSpec(): array;

    /**
     * Returns an array containing rules for input validation.
     *
     * @param  \Illuminate\Http\Request  $request
     *
     * @return array
     */
    abstract public static function validationSpec(Request $request): array;

    /**
     * Returns an array of 5 functions one for each field of the sample group composition script input
     *
     * @return callable[]|null
     */
    public static function sampleGroupFunctions(): ?array
    {
        return null;
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        return true;
    }

    /**
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    abstract public function handle(): void;

    /**
     * Returns the number of threads used by this job
     *
     * @return int
     */
    public function threads(): int
    {
        return (int)$this->getParameter('threads', $this->defaultNumberOfThreads);
    }

    /**
     * Get the output of this job and modifies it so that it can be returned as an output to an API call
     * The default implementation does not mutate the output array
     *
     * @return array
     */
    public function mutateOutput(): array
    {
        return $this->model->job_output;
    }

    /**
     * Actions to clean up everything if the job failed.
     * The default implementation does nothing.
     */
    public function cleanupOnFail(): void
    {
    }

    /**
     * Append text to the log of this job
     *
     * @param  string  $text
     * @param  bool  $appendNewLine
     * @param  bool  $commit
     */
    public function log(string $text, bool $appendNewLine = true, bool $commit = true): void
    {
        $this->model->appendLog($text, $appendNewLine, $commit);
    }

    /**
     * Move a file from a job directory to any destination
     *
     * @param  string|null  $source
     * @param  string  $destination
     */
    public function moveFile(?string $source, string $destination): void
    {
        if (!$source) {
            return;
        }
        $absoluteSourceFilename = realpath($this->model->getAbsoluteJobDirectory() . '/' . $source);
        if ($absoluteSourceFilename === false) {
            return;
        }
        @rename($absoluteSourceFilename, $destination);
        @chmod($destination, 0777);
    }

    /**
     * Add the map path of a reference
     *
     * @param  array  $command
     * @param  \App\Models\Annotation|\App\Models\Reference  $annotation
     *
     * @return bool
     */
    public static function addMap(array &$command, $annotation): bool
    {
        if ($annotation->map_path !== null) {
            $command[] = '-x';
            $command[] = $annotation->map_path;

            return true;
        }

        return false;
    }

    /**
     * Checks if a file is exists in the current job directory
     *
     * @param  string  $file
     *
     * @return bool
     */
    protected function fileExistsRelative(string $file): bool
    {
        return $this->publicFolder->exists($file);
    }

    /**
     * Checks if a file is exists in the current job directory
     *
     * @param  string|null  $file
     *
     * @return bool
     */
    protected function fileExists(?string $file): bool
    {
        if (empty($file)) {
            return false;
        }

        return $this->fileExistsRelative($this->model->getJobDirectory() . '/' . $file);
    }

    /**
     * Checks if a parameter contains the name of a valid file.
     * A file is considered valid only if it exists in the current job directory
     *
     * @param  string  $parameter
     *
     * @return bool
     */
    protected function validateFileParameter(string $parameter): bool
    {
        $file = $this->model->getParameter($parameter);

        return $this->fileExists($file);
    }

    /**
     * Returns the real path of a script
     *
     * @param  string  $script
     *
     * @return string
     */
    public static function scriptPath(string $script): string
    {
        return realpath(config('rnadetector.scripts_path') . '/' . $script);
    }

    /**
     * Helper function used to generate absolute path and url of a relative path
     *
     * @param  string  $pathRelative
     *
     * @return array
     */
    private function pathHelper(string $pathRelative): array
    {
        $pathAbsolute = $this->model->absoluteJobPath($pathRelative);
        $url = $this->publicFolder->url($pathRelative);

        return [$pathRelative, $pathAbsolute, $url];
    }

    /**
     * Returns the paths and url of a new file in the job directory
     * The result is an array with 3 elements: [0] relative path; [1] absolute path; [2] url.
     *
     * @param  string  $prefix
     * @param  string  $suffix
     *
     * @return array
     */
    protected function getJobFilePaths(string $prefix = '', string $suffix = ''): array
    {
        return $this->pathHelper($this->model->getJobFile($prefix, $suffix));
    }

    /**
     * This function returns the relative path, absolute path, and url of a file.
     * The result is an array with 3 elements: [0] relative path; [1] absolute path; [2] url.
     *
     * @param  string  $filename
     *
     * @return array
     */
    protected function getFilePaths(string $filename): array
    {
        $filename = str_ireplace($this->model->getAbsoluteJobDirectory(), '', $filename);

        return $this->pathHelper($this->model->getJobDirectory() . '/' . $filename);
    }

    /**
     * Returns the filename without any extension
     *
     * @param  string  $filename
     *
     * @return string
     */
    protected function getFileNameWithoutExtension(string $filename): string
    {
        $parts = pathinfo($filename);

        return $parts['dirname'] . '/' . $parts['filename'];
    }

    /**
     * This function returns the relative path and url of a file to be saved as output of a job.
     *
     * @param  string|array|null  $filename
     *
     * @return array|null
     */
    protected function getFilePathsForOutput($filename): ?array
    {
        if (!$filename) {
            return null;
        }
        if (is_array($filename) && count($filename) === 3) {
            if (!$this->fileExistsRelative($filename[0])) {
                return null;
            }

            return ['path' => $filename[0], 'url' => $filename[2]];
        }
        [$path, , $url] = $this->pathHelper($filename);
        if (!$this->fileExistsRelative($path)) {
            return null;
        }

        return compact('path', 'url');
    }

    /**
     * Runs a shell command and checks for successful completion of execution
     *
     * @param  array  $command
     * @param  string|null  $cwd
     * @param  int|null  $timeout
     * @param  callable|null  $callback
     * @param  array  $errorCodeMap
     *
     * @return string
     * @throws \App\Exceptions\ProcessingJobException
     */
    public static function runCommand(
        array $command,
        ?string $cwd = null,
        ?int $timeout = null,
        ?callable $callback = null,
        array $errorCodeMap = []
    ): string {
        try {
            return Utils::runCommand($command, $cwd, $timeout, $callback);
        } catch (ProcessFailedException $e) {
            throw Utils::mapCommandException($e, $errorCodeMap);
        }
    }

    /**
     * Stop this job and resubmit after the specified amount of minutes
     *
     * @param  int  $minutes
     *
     * @return void
     */
    protected function resubmit($minutes = 5): void
    {
        throw (new ResubmitException())->setAfter($minutes);
    }

    /**
     * @param  string  $name
     * @param  array  $arguments
     *
     * @return mixed
     */
    public function __call(string $name, array $arguments)
    {
        if (method_exists($this->model, $name)) {
            return $this->model->$name(...$arguments);
        }
        throw new \RuntimeException('Undefined method ' . $name);
    }


}

