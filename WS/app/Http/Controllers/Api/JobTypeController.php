<?php

namespace App\Http\Controllers\Api;

use App\Http\Controllers\Controller;
use App\Jobs\Types\Factory;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;

class JobTypeController extends Controller
{
    /**
     * Display a listing of the resource.
     *
     * @return \Illuminate\Http\JsonResponse
     */
    public function index(): JsonResponse
    {
        return response()->json([
            'data'  => Factory::listTypes()->keyBy('id'),
            'links' => [
                'self' => route('job-types.index'),
            ],
        ]);
    }

    /**
     * Display the specified resource.
     *
     * @param string $type
     * @return \Illuminate\Http\JsonResponse
     */
    public function show(string $type): JsonResponse
    {
        $types = Factory::listTypes();
        $res   = $types->where('id', '=', $type)->first();
        if (!$res) {
            abort(404, 'No query results for type ' . $type);
        }
        $id                = $res['id'];
        $res['parameters'] = Factory::parametersSpec($id);
        $res['output']     = Factory::outputSpec($id);
        return response()->json([
            'data'  => $res,
            'links' => [
                'self' => route('job-types.show', $id),
            ],
        ]);
    }

}
