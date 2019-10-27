<?php

use Illuminate\Http\Request;

/*
|--------------------------------------------------------------------------
| API Routes
|--------------------------------------------------------------------------
|
| Here is where you can register API routes for your application. These
| routes are loaded by the RouteServiceProvider within a group which
| is assigned the "api" middleware group. Enjoy building your API!
|
*/

Route::apiResource('users', 'Api\\UserController')->names([
    'show' => 'users.show',
])->middleware('auth:api');

Route::middleware([
    'auth:api',
    'can:generate-token,user',
])->get('/users/{user}/token', 'Api\\UserController@token');

Route::middleware('auth:api')->get('/user', function (Request $request) {
    return $request->user();
});
