<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Http\Controllers\Api;

use App\Http\Controllers\Controller;
use App\Http\Resources\User as UserResource;
use App\Http\Resources\UserCollection;
use App\Models\User;
use Auth;
use Carbon\Carbon;
use Hash;
use Illuminate\Http\JsonResponse;
use Illuminate\Http\Request;
use Illuminate\Http\Response;
use Illuminate\Support\Str;
use Illuminate\Validation\Rule;

class UserController extends Controller
{
    /**
     * UserController constructor.
     */
    public function __construct()
    {
        $this->authorizeResource(User::class, 'user');
    }


    /**
     * Display a listing of the resource.
     *
     * @param \Illuminate\Http\Request $request
     * @return \App\Http\Resources\UserCollection
     */
    public function index(Request $request): UserCollection
    {
        $perPage = (int)($request->get('per_page') ?? 15);
        if ($perPage < 0) {
            $perPage = 15;
        }
        return new UserCollection(User::paginate($perPage));
    }

    /**
     * Store a newly created resource in storage.
     *
     * @param \Illuminate\Http\Request $request
     * @return \App\Http\Resources\User
     * @throws \Illuminate\Validation\ValidationException
     */
    public function store(Request $request): UserResource
    {
        $values                   = $this->validate($request, [
            'name'     => ['required', 'max:255'],
            'email'    => ['required', 'max:255', 'email', 'unique:users'],
            'password' => ['required', 'string', 'min:6', 'alpha_num'],
            'admin'    => ['filled', 'boolean'],
        ]);
        $model                    = User::create([
            'name'      => $values['name'],
            'email'     => $values['email'],
            'password'  => Hash::make($values['password']),
            'admin'     => $values['admin'] ?? false,
            'api_token' => null,
        ]);
        $model->email_verified_at = Carbon::now();
        $model->save();
        return new UserResource($model);
    }

    /**
     * Display the specified resource.
     *
     * @param \App\Models\User $user
     * @return \App\Http\Resources\User
     */
    public function show(User $user): UserResource
    {
        return new UserResource($user);
    }

    /**
     * Update the specified resource in storage.
     *
     * @param \Illuminate\Http\Request $request
     * @param \App\Models\User         $user
     * @return \App\Http\Resources\User
     * @throws \Illuminate\Validation\ValidationException
     */
    public function update(Request $request, User $user): UserResource
    {
        $rules = [
            'name'         => ['filled', 'max:255'],
            'email'        => ['filled', 'max:255', 'email', Rule::unique('users')->ignore($user->id)],
            'password'     => ['required_with_all:new_password', 'password:api'],
            'new_password' => ['filled', 'string', 'min:6', 'alpha_num'],
            'admin'        => ['filled', 'boolean'],
        ];
        if (Auth::guard('api')->user()->admin) {
            unset($rules['password']);
        }
        $values = $this->validate($request, $rules);
        if (isset($values['name'])) {
            $user->name = $values['name'];
        }
        if (isset($values['email'])) {
            $user->email             = $values['email'];
            $user->email_verified_at = Carbon::now();
        }
        if (isset($values['admin'])) {
            $user->admin = (bool)$values['admin'];
        }
        if (isset($values['new_password'])) {
            $user->password = Hash::make($values['new_password']);
        }
        $user->save();
        return new UserResource($user);
    }

    /**
     * Remove the specified resource from storage.
     *
     * @param \App\Models\User $user
     * @return \Illuminate\Http\JsonResponse
     * @throws \Exception
     */
    public function destroy(User $user): JsonResponse
    {
        $user->delete();
        return response()->json([
            'message' => 'User deleted.',
            'errors'  => false,
        ]);
    }

    public function token(User $user): JsonResponse
    {
        $token = Str::random(80);
        $user->forceFill([
            'api_token' => hash('sha256', $token),
        ])->save();
        return response()->json([
            'data'  => [
                'id'        => $user->id,
                'api_token' => $token,
            ],
            'links' => [
                'self' => route('users.show', $user),
            ],
        ]);
    }
}
