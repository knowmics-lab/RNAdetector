<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\ResourceCollection;

class ReferenceCollection extends ResourceCollection
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
                /** @var Reference $item */
                return [
                    'id'              => $item->id,
                    'name'            => $item->name,
                    'available_for'   => $item->available_for,
                    'has_map'         => $item->map_path !== null && file_exists($item->map_path),
                    'created_at'      => $item->created_at,
                    'created_at_diff' => $item->created_at->diffForHumans(),
                    'updated_at'      => $item->updated_at,
                    'updated_at_diff' => $item->updated_at->diffForHumans(),
                ];
            }
        )->keyBy('id')->all();
    }
}
