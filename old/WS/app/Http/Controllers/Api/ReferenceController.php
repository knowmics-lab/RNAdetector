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
use App\Packages;
use Illuminate\Database\Eloquent\Builder;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;
use Throwable;

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
     * @param  \Illuminate\Http\Request  $request
     *
     * @return \App\Http\Resources\ReferenceCollection
     */
    public function index(Request $request): ReferenceCollection
    {
        return new ReferenceCollection(
            $this->handleBuilderRequest(
                $request,
                Reference::query(),
                static function (Builder $builder) use ($request) {
                    if ($request->has('indexed_for')) {
                        $algo = $request->get('indexed_for');
                        if ($algo) {
                            $builder->whereRaw('JSON_CONTAINS(`available_for`, "true", ?)', ['$.' . $algo]);
                        }
                    }

                    return $builder;
                }
            )
        );
    }

    /**
     * Display the specified resource.
     *
     * @param  \App\Models\Reference  $reference
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
     * @param  \App\Models\Reference  $reference
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

    /**
     * Lists all available packages
     *
     * @return \Illuminate\Http\JsonResponse
     */
    public function listPackages(): JsonResponse
    {
        $data = null;
        try {
            $packages = new Packages();

            $data = $packages->fetchNotInstalled();
        } catch (Throwable $e) {
            abort(500, $e->getMessage());
        }

        return response()->json($data);
    }

}
