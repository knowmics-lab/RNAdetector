<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Models\Annotation;
use App\Models\Reference;
use Illuminate\Console\Command;

class ImportReference extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'reference:import {name}';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Import a new indexed genome/transcriptome and its annotations';

    /**
     * Checks if a filename is valid
     *
     * @param $value
     *
     * @return bool
     */
    public function isValidFilename($value): bool
    {
        return !empty($value) && preg_match('/^[\pL\pM\pN_-]+$/u', $value) > 0;
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        $name = $this->argument('name');
        if (!$this->isValidFilename($name)) {
            $this->error('Reference sequence name is not valid!');

            return 1;
        }
        $genomePath = env('REFERENCES_PATH') . '/' . $name . '/';
        $configFile = $genomePath . '/spec.json';
        if (!file_exists($genomePath) && !is_dir($genomePath)) {
            $this->error('Reference sequence directory not found!');

            return 2;
        }
        if (!file_exists($configFile)) {
            $this->error('Invalid reference sequence: config file not found.');

            return 3;
        }
        if (!file_exists($genomePath . 'reference.fa')) {
            $this->error('Invalid reference sequence: fasta file not found.');

            return 4;
        }

        $config = json_decode(file_get_contents($configFile), true);

        $indexedFor = (array)($config['indexedFor'] ?? []);
        Reference::create(
            [
                'name'          => $name,
                'path'          => $genomePath . 'reference.fa',
                'available_for' => [
                    'bwa'    => $indexedFor['bwa'] ?? false,
                    'tophat' => $indexedFor['tophat'] ?? false,
                    'salmon' => $indexedFor['salmon'] ?? false,
                    'hisat'  => $indexedFor['hisat'] ?? false,
                ],
            ]
        )->save();

        $annotations = (array)($config['annotations'] ?? []);
        foreach ($annotations as $annotationSpec) {
            $annotation = $annotationSpec['name'];
            $type = $annotationSpec['type'] ?? Annotation::GTF;
            if ($type !== Annotation::BED && $type !== Annotation::GTF) {
                $this->warn('Found invalid annotation type for ' . $annotation . ': ' . $type . '.');
                continue;
            }
            if (!$this->isValidFilename($annotation)) {
                $this->warn('Found invalid annotation name: ' . $annotation . '.');
                continue;
            }
            $annotationSourceFile = $genomePath . $annotation . '.' . $type;
            if (!file_exists($annotationSourceFile)) {
                $this->warn('Source path of annotation ' . $annotation . ' not found.');
                continue;
            }
            $annotationDestinationFile = env('ANNOTATIONS_PATH') . '/' . $annotation . '.' . $type;
            @rename($annotationSourceFile, $annotationDestinationFile);
            if (!file_exists($annotationDestinationFile)) {
                $this->warn('Unable to write annotation file for ' . $annotation . '.');
                continue;
            }
            Annotation::create(
                [
                    'name' => $name,
                    'type' => $type,
                    'path' => $annotationDestinationFile,
                ]
            )->save();
        }
        @unlink($configFile);
        $this->info('Reference sequence imported correctly!');

        return 0;
    }
}
