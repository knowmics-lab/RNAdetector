<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types\Traits;


use App\Exceptions\ProcessingJobException;
use App\Models\Reference;

trait UseTranscriptome
{

    /**
     * @var Reference|null
     */
    private $transcriptome = null;

    /**
     * Checks if a valid transcriptome has been provided
     *
     * @param \App\Models\Reference|null $transcriptome
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function checkTranscriptome(?Reference $transcriptome): void
    {
        if ($transcriptome === null) {
            throw new ProcessingJobException('An invalid transcriptome has been provided');
        }
    }

    /**
     * Get the current transcriptome
     *
     * @param string $defaults
     * @param bool   $checks
     *
     * @return \App\Models\Reference
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function getTranscriptome(string $defaults = 'HUMAN_TRANSCRIPTOME_NAME', bool $checks = true): Reference
    {
        if ($this->transcriptome === null) {
            $transcriptomeName = $this->getParameter('transcriptome', env($defaults, $defaults));
            $this->transcriptome = Reference::whereName($transcriptomeName)->first();
        }
        if ($checks) {
            $this->checkTranscriptome($this->transcriptome);
        }

        return $this->transcriptome;
    }

}
