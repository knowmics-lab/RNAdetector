<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Models\Annotation;
use App\Models\Reference;
use App\Utils;
use DirectoryIterator;
use Illuminate\Console\Command;
use Symfony\Component\Process\Process;
use ZipArchive;

class ImportReference extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'reference:import {filename : A zip archive containing the reference and its annotations}';

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
     * Import a reference genome and its annotations
     *
     * @param string $name
     *
     * @return int
     */
    private function import(string $name): int
    {
        if (!$this->isValidFilename($name)) {
            $this->error('Reference sequence name is not valid!');

            return 1;
        }
        if (Reference::whereName($name)->count() > 0) {
            $this->error('Another reference sequence with this name already exists!');

            return 2;
        }
        $genomePath = env('REFERENCES_PATH') . '/' . $name . '/';
        $configFile = $genomePath . 'spec.json';
        if (!file_exists($genomePath) && !is_dir($genomePath)) {
            $this->error('Reference sequence directory not found!');

            return 3;
        }
        if (!file_exists($configFile)) {
            $this->error('Invalid reference sequence: config file not found.');

            return 4;
        }
        if (!file_exists($genomePath . 'reference.fa')) {
            $this->error('Invalid reference sequence: fasta file not found.');

            return 5;
        }

        $config = json_decode(file_get_contents($configFile), true);

        $indexedFor = (array)($config['indexedFor'] ?? []);
        Reference::create(
            [
                'name'          => $name,
                'path'          => realpath($genomePath . 'reference.fa'),
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
            if (Annotation::whereName($annotation)->count() > 0) {
                $this->error('Another reference sequence with this name already exists!');

                return 2;
            }
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
                    'name' => $annotation,
                    'type' => $type,
                    'path' => realpath($annotationDestinationFile),
                ]
            )->save();
        }
        @unlink($configFile);
        $this->info('Reference sequence imported correctly!');

        return 0;
    }

    /**
     * Extract a ZIP file
     *
     * @param string $filename
     *
     * @return void
     */
    private function extract(string $filename): void
    {
        Utils::runCommand(
            [
                'tar',
                '-jxvf',
                $filename,
            ],
            env('REFERENCES_PATH'),
            null,
            function ($type, $buffer) {
                if ($type === Process::OUT) {
                    $this->info(trim($buffer));
                }
            }
        );
    }

    /**
     * Recursively apply chmod to a folder
     *
     * @param string $path
     * @param int    $mode
     */
    private function recursiveChmod(string $path, int $mode): void
    {
        $files = new DirectoryIterator(realpath($path));
        foreach ($files as $file) {
            @chmod($file->getRealPath(), $mode);
            if ($file->isDir() && !$file->isDot()) {
                $this->recursiveChmod($file->getRealPath(), $mode);
            }
        }
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        $filename = $this->argument('filename');
        $ext = pathinfo($filename, PATHINFO_EXTENSION);
        $name = basename($filename, '.' . $ext);
        $this->info('Extracting reference archive...');
        $this->extract(env('REFERENCES_PATH') . '/' . $filename);
        $this->recursiveChmod(env('REFERENCES_PATH') . '/' . $name, 0777);
        $this->info('Importing reference sequence...');

        return $this->import($name);
    }
}
