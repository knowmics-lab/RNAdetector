<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Models\Annotation;
use App\Models\Reference;
use App\Packages;
use App\Utils;
use DirectoryIterator;
use Illuminate\Console\Command;
use Symfony\Component\Process\Process;

class ImportReference extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'reference:import {name : The name of the reference sequence to import}';

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
     * @param string $version
     * @param bool   $isUpdating
     *
     * @return int
     */
    private function import(string $name, string $version, bool $isUpdating): int
    {
        if (!$this->isValidFilename($name)) {
            $this->error('Reference sequence name is not valid!');

            return 11;
        }
        if (!$isUpdating && Reference::whereName($name)->count() > 0) {
            $this->error('Another reference sequence with this name already exists!');

            return 12;
        }
        $genomePath = config('rnadetector.reference_path') . '/' . $name . '/';
        $configFile = $genomePath . 'spec.json';
        if (!file_exists($genomePath) && !is_dir($genomePath)) {
            $this->error('Reference sequence directory not found!');

            return 13;
        }
        if (!file_exists($configFile)) {
            $this->error('Invalid reference sequence: config file not found.');

            return 14;
        }
        $config = json_decode(file_get_contents($configFile), true);

        $hasReference = !((bool)($config['noReference'] ?? false));
        if ($hasReference) {
            if (!file_exists($genomePath . 'reference.fa')) {
                $this->error('Invalid reference sequence: fasta file not found.');

                return 15;
            }
            $indexedFor = (array)($config['indexedFor'] ?? []);
            $mapFile = (bool)($config['mapFile'] ?? false);
            $mapPath = null;
            if ($mapFile) {
                $realMapPath = realpath($genomePath . 'map_file.tsv');
                if ($realMapPath !== false && file_exists($realMapPath)) {
                    $mapPath = $realMapPath;
                }
            }
            Reference::updateOrCreate(
                [
                    'name' => $name,
                ],
                [
                    'path'          => realpath($genomePath . 'reference.fa'),
                    'available_for' => [
                        'bwa'    => $indexedFor['bwa'] ?? false,
                        'hisat'  => $indexedFor['hisat'] ?? false,
                        'salmon' => $indexedFor['salmon'] ?? false,
                        'star'   => $indexedFor['star'] ?? false,
                    ],
                    'map_path'      => $mapPath,
                ]
            )->save();
        }
        $this->info('Reference sequence imported correctly!');

        $this->info('Importing annotations...');
        $annotations = (array)($config['annotations'] ?? []);
        foreach ($annotations as $annotationSpec) {
            $annotation = $annotationSpec['name'];
            if (!$isUpdating && Annotation::whereName($annotation)->count() > 0) {
                $this->error('Another reference annotation with the name "' . $name . '" already exists!');

                return 16;
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
            $annotationDestinationFile = config('rnadetector.annotations_path') . '/' . $annotation . '.' . $type;
            @rename($annotationSourceFile, $annotationDestinationFile);
            if (!file_exists($annotationDestinationFile)) {
                $this->warn('Unable to write annotation file for ' . $annotation . '.');
                continue;
            }
            $annotationMapPath = null;
            $hasMapPath = $annotationSpec['mapFile'] ?? false;
            if ($hasMapPath) {
                $annotationMapPathSourceFile = $genomePath . $annotation . '_map.tsv';
                if (!file_exists($annotationMapPathSourceFile)) {
                    $this->warn('Source path of map file for ' . $annotation . ' not found.');
                } else {
                    $annotationMapPathDestinationFile = config('rnadetector.annotations_path') . '/' . $annotation . '_map.tsv';
                    @rename($annotationMapPathSourceFile, $annotationMapPathDestinationFile);
                    if (!file_exists($annotationMapPathDestinationFile)) {
                        $this->warn('Unable to write map file for ' . $annotation . '.');
                    } else {
                        $annotationMapPath = $annotationMapPathDestinationFile;
                    }
                }
            }
            $annModel = Annotation::updateOrCreate(
                [
                    'name' => $annotation,
                ],
                [
                    'type'     => $type,
                    'path'     => realpath($annotationDestinationFile),
                    'map_path' => $annotationMapPath,
                ]
            );
            $gff3File = $genomePath . $annotation . '.gff3.gz';
            if (file_exists($gff3File)) {
                @rename($gff3File, $annModel->getGFF3Path());
                if (!$annModel->hasGFF3()) {
                    $this->warn('Unable to write GFF3 file for ' . $annotation . '.');
                }
            }
            $tbiFile = $genomePath . $annotation . '.gff3.gz.tbi';
            if (file_exists($tbiFile)) {
                @rename($tbiFile, $annModel->getGFF3Path() . '.tbi');
                if (!file_exists($annModel->getGFF3Path() . '.tbi')) {
                    $this->warn('Unable to write GFF3 index file for ' . $annotation . '.');
                }
            }
            $annModel->save();
        }
        $this->info('Annotations imported correctly!');
        @unlink($configFile);
        @touch($genomePath . '/.installed');
        @chmod($genomePath . '/.installed', 0777);
        @file_put_contents($genomePath . '/.version', $version);
        @chmod($genomePath . '/.version', 0777);

        return 0;
    }

    /**
     * Extract a .tar.bz2 archive
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
            config('rnadetector.reference_path'),
            null,
            function ($type, $buffer) {
                if ($type === Process::OUT) {
                    $this->output->write('<info>' . $buffer . '</info>');
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
        $this->info('Changing permissions on ' . $path . '...', 'v');
        $files = new DirectoryIterator($path);
        foreach ($files as $file) {
            @chmod($file->getRealPath(), $mode);
            if ($file->isDir() && !$file->isDot()) {
                $this->recursiveChmod($file->getRealPath(), $mode);
            }
        }
    }

    /**
     * Recursively delete a folder
     *
     * @param string $dir
     */
    private function recursiveRmdir(string $dir): void
    {
        if (is_dir($dir)) {
            $objects = scandir($dir);
            foreach ($objects as $object) {
                if ($object !== "." && $object !== "..") {
                    if (is_dir($dir . DIRECTORY_SEPARATOR . $object) && !is_link($dir . DIRECTORY_SEPARATOR . $object)) {
                        $this->recursiveRmdir($dir . DIRECTORY_SEPARATOR . $object);
                    } else {
                        unlink($dir . DIRECTORY_SEPARATOR . $object);
                    }
                }
            }
            rmdir($dir);
        }
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
    {
        $name = $this->argument('name');
        $referenceDir = config('rnadetector.reference_path') . '/' . $name;
        $filename = $referenceDir . '.tar.bz2';
        if (!file_exists($filename)) {
            $this->error('Unable to find the reference archive.');

            return 1;
        }
        $packages = new Packages();
        $isUpdating = false;
        if (Packages::isPackageInstalled($name)) {
            if (!$packages->canBeUpdated($name)) {
                $this->error('A package with the same name has already been installed and no update is available.');

                return 2;
            }

            $this->info('Package already exists but an update has been found...installing update...');
            $this->info('Removing old package directory...');
            $this->recursiveRmdir($referenceDir);
            $isUpdating = true;
        }
        $this->info('Extracting reference archive...');
        $this->extract($filename);
        $this->recursiveChmod($referenceDir, 0777);
        $package = $packages->fetchOne($name);
        $version = $package['version'] ?? '';
        $this->info('Importing reference sequence...');

        return $this->import($name, $version, $isUpdating);
    }
}
