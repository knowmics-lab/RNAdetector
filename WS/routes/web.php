<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */
/*
|--------------------------------------------------------------------------
| Web Routes
|--------------------------------------------------------------------------
|
| Here is where you can register web routes for your application. These
| routes are loaded by the RouteServiceProvider within a group which
| contains the "web" middleware group. Now create something great!
|
*/

Route::get('/', 'IndexController@index');

Auth::routes(
    [
        'register' => false,
        'reset'    => false,
        'confirm'  => false,
        'verify'   => false,
    ]
);

Route::get('/home', 'HomeController@index')->name('home');
Route::get('/check-update', 'HomeController@doUpdate')->name('run-update')->middleware('can:view-any,App\\Models\\User');

Route::get('/user/reset-token', 'HomeController@resetToken')->name('reset-token');
Route::get('/user/change-password', 'HomeController@changePasswordForm')->name('change-password');
Route::post('/user/change-password', 'HomeController@doChangePassword')->name('do-change-password');
Route::get('/admin/users', 'UserController@index')->name('users-list')->middleware('can:view-any,App\\Models\\User');
Route::get('/admin/users/new', 'UserController@create')->name('users-create')->middleware(
    'can:create,App\\Models\\User'
);
Route::post('/admin/users/new', 'UserController@doCreate')->name('users-do-create')->middleware(
    'can:create,App\\Models\\User'
);
Route::get('/admin/users/{user}', 'UserController@show')->name('users-show')->middleware('can:view,user');
Route::post('/admin/users/{user}', 'UserController@update')->name('users-update')->middleware('can:update,user');
Route::get('/admin/users/{user}/delete', 'UserController@delete')->name('users-delete')->middleware('can:delete,user');
