<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\ResourceCollection;

class JobCollection extends ResourceCollection
{

    /**
     * Transform the resource collection into an array.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return array
     */
    public function toArray($request)
    {
        return $this->collection->map(
            static function (Job $item) {
                return [
                    'id'              => $item->id,
                    'name'            => $item->name,
                    'type'            => $item->job_type,
                    'readable_type'   => $item->resource->readableJobType(),
                    'status'          => $item->status,
                    'created_at'      => $item->created_at,
                    'created_at_diff' => $item->created_at->diffForHumans(),
                    'updated_at'      => $item->updated_at,
                    'updated_at_diff' => $item->updated_at->diffForHumans(),
                    'owner'           => new User($item->user),
                    'self'            => route('jobs.show', $item->resource),
                    'upload'          => route('jobs.upload', $item->resource),
                    'submit'          => route('jobs.submit', $item->resource),
                ];
            }
        )->all();
    }
}
