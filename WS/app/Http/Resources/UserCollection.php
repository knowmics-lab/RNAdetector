<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\ResourceCollection;

class UserCollection extends ResourceCollection
{
    /**
     * Transform the resource collection into an array.
     *
     * @param \Illuminate\Http\Request $request
     * @return array
     */
    public function toArray($request)
    {
        return $this->collection->map(function ($item) {
            return [
                'id'         => $item->id,
                'name'       => $item->name,
                'email'      => $item->email,
                'admin'      => $item->admin,
                'created_at' => $item->created_at,
                'updated_at' => $item->updated_at,
            ];
        })->all();
    }
}
