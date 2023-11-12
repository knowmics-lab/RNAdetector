<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use Illuminate\Console\Command;

class MakeLinks extends Command
{
    /**
     * The console command signature.
     *
     * @var string
     */
    protected $signature = 'make:links';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Create symbolic links for RNADetector webservice';

    /**
     * Execute the console command.
     *
     * @return int
     * @throws \Illuminate\Contracts\Container\BindingResolutionException
     */
    public function handle(): int
    {
        $files = $this->laravel->make('files');
        if (!file_exists(public_path('references'))) {
            $files->link(
                storage_path('app/references'),
                public_path('references')
            );

            $this->info('The references directory has been linked.');
        }
        if (!file_exists(public_path('annotations'))) {
            $files->link(
                storage_path('app/annotations'),
                public_path('annotations')
            );

            $this->info('The annotations directory has been linked.');
        }

        return 0;
    }
}
