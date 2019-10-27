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

Route::apiResource('jobs', 'Api\\JobController')->names([
    'show' => 'jobs.show',
])->middleware('auth:api');

Route::middleware([
    'auth:api',
    'can:submit-job,job',
])->get('/jobs/{job}/submit', 'Api\\JobController@submit')->name('jobs.submit');

Route::middleware('auth:api')->get('/job-types', 'Api\\JobTypeController@index')->name('job-types.index');
Route::middleware('auth:api')->get('/job-types/{type}', 'Api\\JobTypeController@show')->name('job-types.show');
