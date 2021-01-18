<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Models\Annotation;
use Illuminate\Http\Request;
use Illuminate\Validation\Rule;
use Storage;

class AnnotationUploadJobType extends AbstractJob
{

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'name'     => 'A name for this annotation',
            'type'     => 'The type of this annotation: gtf or bed (Default gtf)',
            'file'     => 'The file for this annotation',
            'map_file' => 'An optional map file (tab-separated with two columns) where each line contains the ID of an annotation and its Entrez Gene Id (Mirbase mature name for miRNAs).',
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
            'done' => 'A boolean that is true if the annotation has been successfully created',
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
            'name'     => ['required', 'alpha_dash', 'max:255'],
            'type'     => ['filled', Rule::in(['gtf', 'bed'])],
            'file'     => ['required', 'string'],
            'map_file' => ['nullable', 'string'],
        ];
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        $file = $this->model->getParameter('file');
        $mapFile = $this->model->getParameter('map_file');
        $disk = Storage::disk('public');
        $dir = $this->model->getJobDirectory() . '/';
        if (!$disk->exists($dir . $file)) {
            return false;
        }
        if ($mapFile && !$disk->exists($dir . $mapFile)) {
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
        $this->log('Starting job.');
        $name = $this->model->getParameter('name');
        $type = $this->model->getParameter('type', 'gtf');
        $file = $this->model->getParameter('file');
        $mapFile = $this->model->getParameter('map_file');
        $annotationFileName = env('ANNOTATIONS_PATH') . '/' . $name . '.' . $type;
        $this->moveFile($file, $annotationFileName);
        if (!file_exists($annotationFileName)) {
            throw new ProcessingJobException('Unable to create annotation file.');
        }
        $mapFileName = null;
        if ($mapFile) {
            $mapFileName = env('ANNOTATIONS_PATH') . '/' . $name . '_map.tsv';
            $this->moveFile($mapFile, $mapFileName);
            if (!file_exists($mapFileName)) {
                throw new ProcessingJobException('Unable to create map file.');
            }
        }
        if ($type === 'gtf') {
            $this->log('Preparing gtf index for the genomic browser...');
            self::runCommand(
                [
                    'bash',
                    self::scriptPath('prepare_gtf.sh'),
                    '-f',
                    $annotationFileName,
                ],
                $this->model->getAbsoluteJobDirectory(),
                null,
                function ($type, $buffer) {
                    $this->log($buffer, false);
                },
                [
                    3 => 'Input file does not exist',
                    4 => 'Output directory is not writable',
                    5 => 'Unable to convert GTF to GFF3',
                    6 => 'Unable to sort GFF3 file',
                    7 => 'Unable to index GFF3 file',
                    8 => 'Unable to remove temporary files'
                ]
            );
        }
        Annotation::create(
            [
                'name'     => $name,
                'type'     => $type,
                'path'     => $annotationFileName,
                'map_path' => $mapFileName,
            ]
        )->save();
        $this->log('Job completed!');
        $this->model->setOutput(['type' => self::OUT_TYPE_CONFIRMATION, 'done' => true]);
        $this->model->save();
    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Upload a novel genome annotation';
    }

    /**
     * @inheritDoc
     */
    public static function displayName(): string
    {
        return 'Annotation Upload';
    }
}
