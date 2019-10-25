<?php

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
        \App\Models\User::create([
            'name'      => 'admin',
            'email'     => 'admin@admin',
            'password'  => Hash::make('password'),
            'api_token' => '',
        ]);
    }
}
