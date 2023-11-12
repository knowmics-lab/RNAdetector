<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Controllers;

use Illuminate\Database\Eloquent\Builder;
use Illuminate\Foundation\Auth\Access\AuthorizesRequests;
use Illuminate\Foundation\Bus\DispatchesJobs;
use Illuminate\Foundation\Validation\ValidatesRequests;
use Illuminate\Http\Request;
use Illuminate\Routing\Controller as BaseController;

class Controller extends BaseController
{
    use AuthorizesRequests, DispatchesJobs, ValidatesRequests;

    /**
     * @param \Illuminate\Http\Request              $request
     * @param \Illuminate\Database\Eloquent\Builder $builder
     * @param callable|null                         $callback
     * @param string                                $defaultOrderField
     * @param string                                $defaultOrdering
     * @param int                                   $defaultPerPage
     *
     * @return \Illuminate\Contracts\Pagination\LengthAwarePaginator|\Illuminate\Database\Eloquent\Builder[]|\Illuminate\Database\Eloquent\Collection
     */
    protected function handleBuilderRequest(
        Request $request,
        Builder $builder,
        ?callable $callback = null,
        string $defaultOrderField = 'created_at',
        string $defaultOrdering = 'desc',
        int $defaultPerPage = 15
    ) {
        $filterBy = $request->get('filter_by');
        $filterValue = $request->get('filter_value');
        if ($filterBy && $filterValue) {
            $builder->where($filterBy, 'LIKE', '%' . $filterValue . '%');
        }
        $orderBy = (array)($request->get('order') ?? [$defaultOrderField]);
        $orderDirection = (array)($request->get('order_direction') ?? [$defaultOrdering]);
        if (!empty($orderBy)) {
            for ($i = 0, $count = count($orderBy); $i < $count; $i++) {
                if ($orderBy[$i]) {
                    $builder->orderBy($orderBy[$i], $orderDirection[$i] ?? $defaultOrdering);
                }
            }
        }
        if ($callback !== null && is_callable($callback)) {
            $builder = $callback($builder);
        }
        $perPage = (int)($request->get('per_page') ?? $defaultPerPage);
        if ($perPage > 0) {
            return $builder->paginate($perPage)->appends($request->input());
        }

        return $builder->get();
    }
}
