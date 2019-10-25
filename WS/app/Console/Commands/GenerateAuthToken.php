<?php

namespace App\Console\Commands;

use App\Models\User;
use Illuminate\Console\Command;
use Illuminate\Support\Str;

class GenerateAuthToken extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'auth:token {user} {--json}';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Generate a new authentication token for an user';

    /**
     * Create a new command instance.
     *
     * @return void
     */
    public function __construct()
    {
        parent::__construct();
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        $user       = $this->argument('user');
        $json       = $this->option('json');
        $userObject = User::whereEmail($user)->first();
        if ($userObject === null) {
            if ($json) {
                $this->line(json_encode([
                    'error' => 101,
                    'data'  => null,
                ]));
            } else {
                $this->error('User not found. Please specify a valid user.');
            }
        } else {
            $token = Str::random(80);
            $userObject->forceFill([
                'api_token' => hash('sha256', $token),
            ])->save();
            if ($json) {
                $this->line(json_encode([
                    'error' => 0,
                    'data'  => $token,
                ]));
            } else {
                $this->info('Token generated. The new token is: ' . $token);
            }
        }
    }
}
