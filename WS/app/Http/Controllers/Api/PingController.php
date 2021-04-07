<?php

namespace App\Http\Controllers\Api;

use App\Http\Controllers\Controller;
use App\SystemInfo;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;

class PingController extends Controller
{
    /**
     * @return \Illuminate\Http\JsonResponse
     */
    public function ping(): JsonResponse
    {
        return response()->json(
            [
                'data' => 'pong',
            ]
        );
    }

    /**
     * @param  \Illuminate\Http\Request  $request
     *
     * @return mixed
     */
    public function user(Request $request)
    {
        return $request->user();
    }

    /**
     *
     * @return \Illuminate\Http\JsonResponse
     */
    public function sysInfo(): JsonResponse
    {
        return response()->json((new SystemInfo())->toArray());
    }
}
