<?php

namespace App\Http\Resources;

use Illuminate\Http\Resources\Json\JsonResource;

class User extends JsonResource
{
    /**
     * Transform the resource into an array.
     *
     * @param \Illuminate\Http\Request $request
     * @return array
     */
    public function toArray($request): array
    {
        return [
            'data'  => [
                'id'         => $this->id,
                'name'       => $this->name,
                'email'      => $this->email,
                'admin'      => $this->admin,
                'created_at' => $this->created_at,
                'updated_at' => $this->updated_at,
                'jobs'       => new JobCollection($this->whenLoaded('jobs')),
            ],
            'links' => [
                'self' => route('users.show', $this->resource),
            ],
        ];
    }
}
