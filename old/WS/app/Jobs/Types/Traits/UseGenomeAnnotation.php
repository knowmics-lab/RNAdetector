<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Models\Annotation;

trait UseGenomeAnnotation
{

    /**
     * @var Annotation|null
     */
    private $genomeAnnotation = null;

    /**
     * Checks if a valid Genome Annotation has been provided
     *
     * @param  \App\Models\Annotation|null  $genomeAnnotation
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkGenomeAnnotation(?Annotation $genomeAnnotation): void
    {
        if ($genomeAnnotation === null || !$genomeAnnotation->isGtf()) {
            throw new ProcessingJobException('An invalid genome annotation has been provided');
        }
    }

    /**
     * Get the current Genome Annotation
     *
     * @param  string  $defaults
     * @param  bool  $checks
     *
     * @return \App\Models\Annotation
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function getGenomeAnnotation(string $defaults = 'human_rna_annotation_name', bool $checks = true): Annotation
    {
        if ($this->genomeAnnotation === null) {
            $annotationName = $this->getParameter('annotation', config('rnadetector.' . $defaults));
            $this->genomeAnnotation = Annotation::whereName($annotationName)->first();
        }
        if ($checks) {
            $this->checkGenomeAnnotation($this->genomeAnnotation);
        }

        return $this->genomeAnnotation;
    }

}
