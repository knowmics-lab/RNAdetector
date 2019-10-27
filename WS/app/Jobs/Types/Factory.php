<?php


namespace App\Jobs\Types;

use App\Exceptions\ProcessingJobException;
use App\Models\Job as JobModel;
use Cache;
use Illuminate\Support\Collection;
use Illuminate\Support\Str;
use Symfony\Component\Finder\Finder;

/**
 * @method static string description(\App\Models\Job|\App\Jobs\Types\AbstractJob|string $where)
 * @method static array parametersSpec(\App\Models\Job|\App\Jobs\Types\AbstractJob|string $where)
 * @method static array outputSpec(\App\Models\Job|\App\Jobs\Types\AbstractJob|string $where)
 * @method static array validationSpec(\App\Models\Job|\App\Jobs\Types\AbstractJob|string $where)
 */
class Factory
{
    /**
     * @var \Illuminate\Support\Collection|null
     */
    private static $_cachedList = null;

    /**
     * Given a Job Model, this method builds an object that will be able to process this job
     *
     * @param \App\Models\Job $jobModel
     * @return \App\Jobs\Types\AbstractJob
     * @throws \App\Exceptions\ProcessingJobException
     */
    public static function get(JobModel $jobModel): AbstractJob
    {
        $jobClass = '\App\Jobs\Types\\' . Str::studly($jobModel->job_type);
        if (!class_exists($jobClass)) {
            throw new ProcessingJobException('Unable to find a class suitable for this job.');
        }
        try {
            $r = new \ReflectionClass($jobClass);
            if ($r->isSubclassOf(AbstractJob::class) && !$r->isAbstract() && $r->getConstructor()->getNumberOfRequiredParameters() === 1) {
                /** @var \App\Jobs\Types\AbstractJob $obj */
                $obj = $r->newInstance($jobModel);
                return $obj;
            }
        } catch (\ReflectionException $e) {
            throw new ProcessingJobException('An error occurred during type class instantiation', 0, $e);
        }
        throw new ProcessingJobException('The type of job ' . $jobModel->id . ' is not valid!');
    }

    /**
     * Implementation of virtual static methods
     *
     * @param string $name
     * @param array  $arguments
     * @return mixed
     * @throws \App\Exceptions\ProcessingJobException
     * @throws \ReflectionException
     */
    public static function __callStatic($name, $arguments)
    {
        if (in_array($name, [
            'description',
            'parametersSpec',
            'outputSpec',
            'validationSpec',
        ])) {
            if (count($arguments) === 1) {
                $where    = $arguments[0];
                $jobClass = null;
                if ($where instanceof JobModel) {
                    $jobClass = '\App\Jobs\Types\\' . Str::studly($where->job_type);
                } elseif ($where instanceof AbstractJob) {
                    $jobClass = get_class($where);
                } elseif (is_string($where) && class_exists($where)) {
                    $r = new \ReflectionClass($where);
                    if ($r->isSubclassOf(AbstractJob::class) && !$r->isAbstract()) {
                        $jobClass = $where;
                    }
                } elseif (is_string($where) && !class_exists($where)) {
                    $jobClass = '\App\Jobs\Types\\' . Str::studly($where);
                }
                if (!class_exists($jobClass)) {
                    throw new ProcessingJobException('Invalid parameter where specified in static method ' . $name . ' from __callStatic.');
                }
                return call_user_func([$jobClass, $name]);
            }
            throw new ProcessingJobException('Undefined parameter in static method ' . $name . ' from __callStatic.');
        }
        throw new ProcessingJobException('Undefined static method ' . $name . ' from __callStatic.');
    }

    /**
     * Returns a list of all available job types
     *
     * @return \Illuminate\Support\Collection
     */
    public static function listTypes(): Collection
    {
        if (!Cache::has('types-list')) {
            $list   = [];
            $finder = new Finder();
            $finder->files()->name('*Type.php')->in(app_path('Jobs/Types/'));
            foreach ($finder as $file) {
                $ns = '\App\Jobs\Types';
                if ($relativePath = $file->getRelativePath()) {
                    $ns .= '\\' . str_replace('/', '\\', $relativePath);
                }
                $class = $ns . '\\' . $file->getBasename('.php');
                try {
                    $r = new \ReflectionClass($class);
                    if ($r->isSubclassOf(AbstractJob::class) && !$r->isAbstract() && $r->getConstructor()->getNumberOfRequiredParameters() === 1) {
                        $list[] = [
                            'id'          => Str::snake($file->getBasename('.php')),
                            'description' => call_user_func([$class, 'description']),
                        ];
                    }
                } catch (\ReflectionException $ignore) {
                    continue;
                }
            }
            $list = collect($list);
            Cache::put('types-list', $list, now()->addMinutes(10)); //List is computed every 10 minutes
            return $list;
        }
        return Cache::get('types-list');
    }

}
