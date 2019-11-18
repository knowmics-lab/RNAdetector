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
                'id'   => $this->id,
                'name' => $this->name,
            ],
            'links' => [
                'self' => route('annotation.show', $this->resource),
            ],
        ];
    }
}
