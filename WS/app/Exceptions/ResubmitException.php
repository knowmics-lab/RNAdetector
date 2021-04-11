<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Exceptions;

use RuntimeException;

class ResubmitException extends RuntimeException
{

    private $after = 5;

    /**
     * Get the number of minutes to wait for resubmission
     *
     * @return int
     */
    public function getAfter(): int
    {
        return $this->after;
    }

    /**
     * Number of minutes to wait for resubmission
     *
     * @param  int  $after
     *
     * @return $this
     */
    public function setAfter(int $after): self
    {
        $this->after = $after;

        return $this;
    }


}
