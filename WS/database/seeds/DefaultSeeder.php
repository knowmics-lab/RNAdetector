<?php

use App\Models\Annotation;
use App\Models\Reference;
use Illuminate\Database\Seeder;

class DefaultSeeder extends Seeder
{
    /**
     * Run the database seeds.
     *
     * @return void
     */
    public function run()
    {
        Reference::create(
            [
                'name'          => env('HUMAN_GENOME_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_GENOME_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => true,
                    'tophat' => true,
                    'salmon' => false,
                ],
            ]
        )->save();
        Reference::create(
            [
                'name'          => env('HUMAN_TRANSCRIPTOME_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_TRANSCRIPTOME_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => false,
                    'tophat' => false,
                    'salmon' => true,
                ],
            ]
        )->save();
        Reference::create(
            [
                'name'          => env('HUMAN_TRANSCRIPTOME_MRNA_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_TRANSCRIPTOME_MRNA_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => false,
                    'tophat' => false,
                    'salmon' => true,
                ],
            ]
        )->save();
        Reference::create(
            [
                'name'          => env('HUMAN_TRANSCRIPTOME_LNCRNA_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_TRANSCRIPTOME_LNCRNA_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => false,
                    'tophat' => false,
                    'salmon' => true,
                ],
            ]
        )->save();
        Reference::create(
            [
                'name'          => env('HUMAN_TRANSCRIPTOME_SNCRNA_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_TRANSCRIPTOME_SNCRNA_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => false,
                    'tophat' => false,
                    'salmon' => true,
                ],
            ]
        )->save();
        Annotation::create(
            [
                'name' => env('HUMAN_CIRI_ANNOTATION_NAME'),
                'path' => env('ANNOTATIONS_PATH') . '/' . env('HUMAN_CIRI_ANNOTATION_NAME') . '.gtf',
            ]
        )->save();
        Annotation::create(
            [
                'name' => env('HUMAN_SNCRNA_ANNOTATION_NAME'),
                'path' => env('ANNOTATIONS_PATH') . '/' . env('HUMAN_SNCRNA_ANNOTATION_NAME') . '.gtf',
            ]
        )->save();
    }
}
