<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Models\Reference;

trait UseGenome
{

    /**
     * @var Reference|null
     */
    private $genome = null;

    /**
     * Checks if a valid genome has been provided
     *
     * @param \App\Models\Reference|null $genome
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkGenome(?Reference $genome): void
    {
        if ($genome === null) {
            throw new ProcessingJobException('An invalid genome has been provided');
        }
    }

    /**
     * Get the current genome
     *
     * @param string $defaults
     * @param bool   $checks
     *
     * @return \App\Models\Reference
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function getGenome(string $defaults = 'HUMAN_GENOME_NAME', bool $checks = true): Reference
    {
        if ($this->genome === null) {
            $genomeName = $this->getParameter('genome', env($defaults, $defaults));
            $this->genome = Reference::whereName($genomeName)->first();
        }
        if ($checks) {
            $this->checkGenome($this->genome);
        }

        return $this->genome;
    }

}
