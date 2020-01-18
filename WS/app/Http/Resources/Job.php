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
     *
     * @return array
     */
    public function toArray($request)
    {
        return [
            'data'  => [
                'id'              => $this->id,
                'sample_code'     => $this->sample_code,
                'name'            => $this->name,
                'type'            => $this->job_type,
                'readable_type'   => $this->resource->readableJobType(),
                'status'          => $this->status,
                'parameters'      => $this->job_parameters,
                'output'          => $this->job_output,
                'log'             => $this->log,
                'created_at'      => $this->created_at,
                'created_at_diff' => $this->created_at->diffForHumans(),
                'updated_at'      => $this->updated_at,
                'updated_at_diff' => $this->updated_at->diffForHumans(),
                'owner'           => new User($this->user),
            ],
            'links' => [
                'self'   => route('jobs.show', $this->resource),
                'owner'  => route('users.show', $this->user),
                'upload' => route('jobs.upload', $this->resource),
                'submit' => route('jobs.submit', $this->resource),
            ],
        ];
    }
}
