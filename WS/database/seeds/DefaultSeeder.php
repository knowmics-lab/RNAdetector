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
    public function run(): void
    {
        Reference::create(
            [
                'name'          => env('HUMAN_GENOME_NAME'),
                'path'          => env('REFERENCES_PATH') . '/' . env('HUMAN_GENOME_NAME') . '/reference.fa',
                'available_for' => [
                    'bwa'    => true,
                    'tophat' => true,
                    'hisat'  => true,
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
                    'hisat'  => false,
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
                    'hisat'  => false,
                    'salmon' => true,
                ],
            ]
        )->save();
        Annotation::create(
            [
                'name' => env('HUMAN_CIRC_ANNOTATION_NAME'),
                'type' => 'gtf',
                'path' => env('ANNOTATIONS_PATH') . '/' . env('HUMAN_CIRC_ANNOTATION_NAME') . '.gtf',
            ]
        )->save();
        Annotation::create(
            [
                'name' => env('HUMAN_SNCRNA_ANNOTATION_NAME'),
                'type' => 'gtf',
                'path' => env('ANNOTATIONS_PATH') . '/' . env('HUMAN_SNCRNA_ANNOTATION_NAME') . '.gtf',
            ]
        )->save();
        Annotation::create(
            [
                'name' => env('HUMAN_RNA_ANNOTATION_NAME'),
                'type' => 'gtf',
                'path' => env('ANNOTATIONS_PATH') . '/' . env('HUMAN_RNA_ANNOTATION_NAME') . '.gtf',
            ]
        )->save();
    }
}
