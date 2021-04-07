<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

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

Route::get('/ping', 'Api\\PingController@ping');

Route::group(
    [
        'middleware' => 'auth:api',
    ],
    static function () {
        Route::get('/auth-ping', 'Api\\PingController@ping');
        Route::get('/sys-info', 'Api\\PingController@sysInfo');

        Route::apiResource('users', 'Api\\UserController')->names(
            [
                'show' => 'users.show',
            ]
        );

        Route::middleware('can:generate-token,user')->get('/users/{user}/token', 'Api\\UserController@token');

        Route::get('/user', 'Api\\PingController@user');

        Route::apiResource('jobs', 'Api\\JobController')->names(
            [
                'show' => 'jobs.show',
            ]
        );

        Route::middleware('can:submit-job,job')->get('/jobs/{job}/submit', 'Api\\JobController@submit')->name('jobs.submit');

        Route::middleware('can:upload-job,job')->any('/jobs/{job}/upload/{any?}', 'Api\\JobController@upload')->where('any', '.*')->name(
            'jobs.upload'
        );

        Route::get('/job-types', 'Api\\JobTypeController@index')->name('job-types.index');
        Route::get('/job-types/{type}', 'Api\\JobTypeController@show')->name('job-types.show');

        Route::apiResource('annotations', 'Api\\AnnotationController')->names(
            [
                'show' => 'annotation.show',
            ]
        )->except(['create', 'store', 'update']);

        Route::get('references/packages', 'Api\\ReferenceController@listPackages');

        Route::apiResource('references', 'Api\\ReferenceController')->names(
            [
                'show' => 'reference.show',
            ]
        )->except(['create', 'store', 'update']);
    }
);

