<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\JsonResource;

/**
 * Class Reference
 * @mixin \App\Models\Reference
 * @package App\Http\Resources
 */
class Reference extends JsonResource
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
                'name'            => $this->name,
                'available_for'   => $this->available_for,
                'has_map'         => $this->map_path !== null && file_exists($this->map_path),
                'created_at'      => $this->created_at,
                'created_at_diff' => $this->created_at->diffForHumans(),
                'updated_at'      => $this->updated_at,
                'updated_at_diff' => $this->updated_at->diffForHumans(),
            ],
            'links' => [
                'self' => route('reference.show', $this->resource),
            ],
        ];
    }
}
