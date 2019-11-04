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
     * Indicates if the resource's collection keys should be preserved.
     *
     * @var bool
     */
    public $preserveKeys = true;

    /**
     * Transform the resource collection into an array.
     *
     * @param \Illuminate\Http\Request $request
     * @return array
     */
    public function toArray($request)
    {
        return $this->collection->map(static function ($item) {
            return [
                'id'         => $item->id,
                'type'       => $item->job_type,
                'status'     => $item->status,
                'created_at' => $item->created_at,
                'updated_at' => $item->updated_at,
                'owner'      => new User($item->user),
                'submit'     => route('jobs.submit', $item),
            ];
        })->keyBy('id')->all();
    }
}
