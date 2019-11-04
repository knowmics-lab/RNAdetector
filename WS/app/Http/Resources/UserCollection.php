<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\ResourceCollection;

class UserCollection extends ResourceCollection
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
                'name'       => $item->name,
                'email'      => $item->email,
                'admin'      => $item->admin,
                'created_at' => $item->created_at,
                'updated_at' => $item->updated_at,
            ];
        })->keyBy('id')->all();
    }
}
