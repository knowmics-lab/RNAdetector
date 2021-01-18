<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */


namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Models\Annotation;
use App\Models\Job;
use App\Models\Reference;
use Illuminate\Support\Str;

trait UsesJBrowseTrait
{

    abstract protected function getJobFilePaths(string $prefix = '', string $suffix = ''): array;

    abstract public function log(string $text, bool $appendNewLine = true, bool $commit = true): void;

    /**
     * Make CIRIquant configuration file and returns its path
     *
     * @param \App\Models\Job        $model
     * @param array|null             $bamPath
     * @param \App\Models\Reference  $reference
     * @param \App\Models\Annotation $annotation
     *
     * @return array|null
     */
    private function makeJBrowseConfig(Job $model, ?array $bamPath, Reference $reference, Annotation $annotation): ?array
    {
        if (!$bamPath) {
            return null;
        }
        $this->log('Generating JBrowse config file...', false);
        if (!$reference->isFastaIndexed()) {
            $this->log('Unable to generate JBrowse config file since the reference genome has not been indexed.');

            return null;
        }
        if (!$annotation->hasGFF3()) {
            $this->log('Unable to generate JBrowse config file since the genome annotation has not been processed for JBrowse2.');

            return null;
        }
        $bamFileUri = $model->getJobUri(str_ireplace($model->getAbsoluteJobDirectory(), '', $bamPath[1]));
        $configFile = $this->getJobFilePaths('jbrowse_config_', '.json');
        $template = file_get_contents(resource_path('templates/jbrowse_config.json'));
        $template = str_replace(
            [
                '{GENOME_NAME}',
                '{GENOME_ID}',
                '{GENOME_PATH}',
                '{ANNOTATION_ID}',
                '{ANNOTATION_NAME}',
                '{ANNOTATION_GFF3_PATH}',
                '{SAMPLE_CODE}',
                '{SAMPLE_NAME}',
                '{SAMPLE_BAM_PATH}',
            ],
            [
                $reference->name,
                Str::slug($reference->name),
                $reference->getReferenceUri(),
                Str::slug($annotation->name),
                $annotation->name,
                $annotation->getGFF3Uri(),
                $model->sample_code,
                $model->name,
                $bamFileUri,
            ],
            $template
        );
        @file_put_contents($configFile[1], $template);
        @chmod($configFile[1], 0777);
        $this->log('Ok!');

        return $configFile;
    }

}
