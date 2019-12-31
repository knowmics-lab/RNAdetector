<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\JsonResource;

class Annotation extends JsonResource
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
                'type'            => $this->type,
                'created_at'      => $this->created_at,
                'created_at_diff' => $this->created_at->diffForHumans(),
                'updated_at'      => $this->updated_at,
                'updated_at_diff' => $this->updated_at->diffForHumans(),
            ],
            'links' => [
                'self' => route('annotation.show', $this->resource),
            ],
        ];
    }
}
