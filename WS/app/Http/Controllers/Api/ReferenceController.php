<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Controllers\Api;

use App\Http\Controllers\Controller;
use App\Http\Resources\Reference as ReferenceResource;
use App\Http\Resources\ReferenceCollection;
use App\Models\Reference;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;

class ReferenceController extends Controller
{

    /**
     * JobController constructor.
     */
    public function __construct()
    {
        $this->authorizeResource(Reference::class, 'reference');
    }


    /**
     * Display a listing of the resource.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return \App\Http\Resources\ReferenceCollection
     */
    public function index(Request $request): ReferenceCollection
    {
        $perPage = (int)($request->get('per_page') ?? 15);
        if ($perPage < 0) {
            $perPage = 15;
        }

        return new ReferenceCollection(Reference::paginate($perPage));
    }

    /**
     * Display the specified resource.
     *
     * @param \App\Models\Reference $reference
     *
     * @return \App\Http\Resources\Reference
     */
    public function show(Reference $reference): ReferenceResource
    {
        return new ReferenceResource($reference);
    }

    /**
     * Remove the specified resource from storage.
     *
     * @param \App\Models\Reference $reference
     *
     * @return \Illuminate\Http\JsonResponse
     * @throws \Exception
     */
    public function destroy(Reference $reference): JsonResponse
    {
        $reference->delete();

        return response()->json(
            [
                'message' => 'Reference deleted.',
                'errors'  => false,
            ]
        );
    }

}
