<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\ResourceCollection;

class AnnotationCollection extends ResourceCollection
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
            static function ($item) {
                return [
                    'id'   => $item->id,
                    'name' => $item->name,
                ];
            }
        )->keyBy('id')->all();
    }
}
