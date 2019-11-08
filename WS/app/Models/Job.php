<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

use Illuminate\Database\Eloquent\Model;
use Illuminate\Database\Eloquent\Relations\BelongsTo;
use Storage;

/**
 * App\Models\Job
 *
 * @property int                             $id
 * @property string                          $job_type
 * @property string                          $status
 * @property array                           $job_parameters
 * @property array                           $job_output
 * @property string                          $log
 * @property int                             $user_id
 * @property \Illuminate\Support\Carbon|null $created_at
 * @property \Illuminate\Support\Carbon|null $updated_at
 * @property-read \App\Models\User           $user
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job newModelQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job newQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job query()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereCreatedAt($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereId($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereJobOutput($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereJobParameters($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereJobType($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereLog($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereStatus($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereUpdatedAt($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Job whereUserId($value)
 * @mixin \Eloquent
 */
class Job extends Model
{

    public const READY      = 'ready';
    public const QUEUED     = 'queued';
    public const PROCESSING = 'processing';
    public const COMPLETED  = 'completed';
    public const FAILED     = 'failed';

    /**
     * The attributes that are mass assignable.
     *
     * @var array
     */
    protected $fillable = [
        'job_type',
        'status',
        'job_parameters',
        'job_output',
        'log',
        'user_id',
    ];

    /**
     * The attributes that should be cast to native types.
     *
     * @var array
     */
    protected $casts = [
        'job_parameters' => 'array',
        'job_output'     => 'array',
    ];

    /**
     * @return \Illuminate\Database\Eloquent\Relations\BelongsTo
     */
    public function user(): BelongsTo
    {
        return $this->belongsTo(User::class, 'user_id', 'id');
    }

    /**
     * Set the status attribute.
     *
     * @param string $value
     */
    public function setStatusAttribute($value): void
    {
        if (!in_array(
            $value,
            [
                self::READY,
                self::QUEUED,
                self::PROCESSING,
                self::COMPLETED,
                self::FAILED,
            ],
            true
        )) {
            $value = self::READY;
        }
        $this->attributes['status'] = $value;
    }

    /**
     * Returns the path of the job storage directory
     *
     * @return string
     */
    public function getJobDirectory(): string
    {
        $path = 'jobs/' . $this->id;
        $disk = Storage::disk('public');
        if (!$disk->exists($path)) {
            $disk->makeDirectory($path);
        }

        return $path;
    }

    /**
     * Returns the absolute path of the job storage directory
     *
     * @return string
     */
    public function getAbsoluteJobDirectory(): string
    {
        return storage_path('app/public/' . $this->getJobDirectory());
    }

    /**
     * Returns the path of a temporary file in the job directory
     *
     * @param string $prefix
     * @param string $suffix
     *
     * @return string
     */
    public function getJobTempFile(string $prefix = '', string $suffix = ''): string
    {
        $filename = preg_replace('/[^\w]+/', '', uniqid($prefix, true)) . $suffix;

        return $this->getJobDirectory() . '/' . $filename;
    }

    /**
     * Returns the absolute path of a temporary file in the job directory
     *
     * @param string $prefix
     * @param string $suffix
     *
     * @return string
     */
    public function getJobTempFileAbsolute(string $prefix = '', string $suffix = ''): string
    {
        return storage_path('app/public/' . $this->getJobTempFile($prefix, $suffix));
    }

    /**
     * Delete the job directory
     *
     * @return bool
     */
    public function deleteJobDirectory(): bool
    {
        return Storage::disk('public')->deleteDirectory($this->getJobDirectory());
    }

    /**
     * Checks if the current job can be modified
     *
     * @return bool
     */
    public function canBeModified(): bool
    {
        return $this->status === self::READY;
    }

    /**
     * Checks if the current job can be deleted
     *
     * @return bool
     */
    public function canBeDeleted(): bool
    {
        return in_array($this->status, [self::READY, self::COMPLETED, self::FAILED], true);
    }

    /**
     * Checks if the current job should run or not.
     * Only queued jobs should run.
     *
     * @return bool
     */
    public function shouldNotRun(): bool
    {
        return in_array($this->status, [self::PROCESSING, self::COMPLETED, self::FAILED], true);
    }

    /**
     * Set a new status value and save this model
     *
     * @param $newStatus
     *
     * @return $this
     */
    public function setStatus($newStatus): self
    {
        $this->status = $newStatus;
        $this->save();

        return $this;
    }

    /**
     * Set the value of one or more output data.
     * If $parameter is an associative array sets multiple parameters at the same time.
     *
     * @param array|string $parameter
     * @param null|mixed   $value
     *
     * @return $this
     */
    public function setOutput($parameter, $value = null): self
    {
        $tmp = $this->job_output;
        if (!is_array($tmp)) {
            $tmp = [];
        }
        if ($value === null && is_array($parameter)) {
            foreach ($parameter as $p => $v) {
                data_set($tmp, $p, $v);
            }
        } else {
            data_set($tmp, $parameter, $value);
        }
        $this->job_output = $tmp;

        return $this;
    }

    /**
     * Get the value of an output data
     *
     * @param string|array|null $parameter
     * @param mixed             $default
     *
     * @return mixed
     */
    public function getOutput($parameter = null, $default = null)
    {
        if ($parameter === null) {
            return $this->job_output;
        }
        if (is_array($parameter)) {
            $slice = [];
            foreach ($parameter as $key) {
                $slice[$key] = data_get($this->job_output, $key, $default);
            }

            return $slice;
        }

        return data_get($this->job_output, $parameter, $default);
    }

    /**
     * Set the value of a parameter
     *
     * @param string $parameter
     * @param mixed  $value
     *
     * @return $this
     */
    public function setParameter(string $parameter, $value): self
    {
        $tmp = $this->job_parameters;
        data_set($tmp, $parameter, $value);
        $this->job_parameters = $tmp;

        return $this;
    }

    /**
     * Add parameters to this job
     *
     * @param array $parameters
     *
     * @return $this
     */
    public function addParameters(array $parameters): self
    {
        $tmp = $this->job_parameters;
        foreach ($parameters as $param => $value) {
            data_set($tmp, $param, $value);
        }
        $this->job_parameters = $tmp;

        return $this;
    }

    /**
     * Set parameters of this job
     *
     * @param array $parameters
     *
     * @return $this
     */
    public function setParameters(array $parameters): self
    {
        $this->job_parameters = [];

        return $this->addParameters($parameters);
    }

    /**
     * Get the value of a parameter
     *
     * @param string|array|null $parameter
     * @param mixed             $default
     *
     * @return mixed
     */
    public function getParameter($parameter = null, $default = null)
    {
        if ($parameter === null) {
            return $this->job_parameters;
        }
        if (is_array($parameter)) {
            $slice = [];
            foreach ($parameter as $key) {
                $slice[$key] = data_get($this->job_parameters, $key, $default);
            }

            return $slice;
        }

        return data_get($this->job_parameters, $parameter, $default);
    }

    /**
     * Append text to the log
     *
     * @param string $text
     * @param bool   $appendNewLine
     * @param bool   $commit
     */
    public function appendLog(string $text, bool $appendNewLine = true, bool $commit = true): void
    {
        if ($appendNewLine) {
            $text .= PHP_EOL;
        }
        //echo $text; // @TODO FOR DEBUG ONLY
        $this->log .= $text;
        if ($commit) {
            $this->save();
        }
    }


}
