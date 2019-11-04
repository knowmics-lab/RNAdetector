<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\JsonResource;

class Job extends JsonResource
{
    /**
     * Transform the resource into an array.
     *
     * @param \Illuminate\Http\Request $request
     * @return array
     */
    public function toArray($request)
    {
        return [
            'data'  => [
                'id'         => $this->id,
                'type'       => $this->job_type,
                'status'     => $this->status,
                'parameters' => $this->job_parameters,
                'output'     => $this->job_output,
                'log'        => $this->log,
                'created_at' => $this->created_at,
                'updated_at' => $this->updated_at,
                'owner'      => new User($this->user),
            ],
            'links' => [
                'self'   => route('jobs.show', $this->resource),
                'owner'  => route('users.show', $this->user),
                'submit' => route('jobs.submit', $this->resource),
            ],
        ];
    }
}
