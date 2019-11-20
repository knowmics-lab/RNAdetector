<?php

use App\Models\User;
use Carbon\Carbon;
use Illuminate\Database\Seeder;

class UsersTableSeeder extends Seeder
{
    /**
     * Run the database seeds.
     *
     * @return void
     */
    public function run()
    {
        $model                    = User::create([
            'name'      => 'admin',
            'email'     => 'admin@admin',
            'password'  => Hash::make('password'),
            'admin'     => true,
            'api_token' => null,
        ]);
        $model->email_verified_at = Carbon::now();
        $model->save();
    }
}
