<?php

use App\Models\User;
use Illuminate\Database\Seeder;

class UsersSeeder extends Seeder
{
    /**
     * Run the database seeds.
     *
     * @return void
     */
    public function run()
    {
        User::create(
            [
                'name'     => 'admin',
                'email'    => 'admin@admin',
                'password' => Hash::make('admin'),
            ]
        );
    }
}
