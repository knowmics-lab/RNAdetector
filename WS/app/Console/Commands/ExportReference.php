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
use RecursiveDirectoryIterator;
use RecursiveIteratorIterator;
use ZipArchive;

class ExportReference extends Command
{
    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'reference:export
                            {reference : The name of the reference sequence}
                            {--annotation=* : A list of annotation to include in the archive}
                            {--output-dir= : The output directory. It defaults to the current directory.}';

    /**
     * The console command description.
     *
     * @var string
     */
    protected $description = 'Export an indexed genome to a zip archive';

    /**
     * @param \App\Models\Reference $refModel
     * @param Annotation[]          $annModel
     * @param string                $baseDir
     *
     * @return string
     */
    private function makeConfigFile(Reference $refModel, array $annModel, string $baseDir): string
    {
        $configFile = $baseDir . '/spec.json';
        $configContent = [
            'indexedFor'  => $refModel->available_for,
            'annotations' => array_map(
                static function (Annotation $model) {
                    return [
                        'name' => $model->name,
                        'type' => $model->type,
                    ];
                },
                $annModel
            ),
        ];
        file_put_contents($configFile, json_encode($configContent));

        return $configFile;
    }

    /**
     * @param Annotation[] $annModel
     * @param string       $baseDir
     *
     * @return string[]
     */
    private function copyAnnotations(array $annModel, string $baseDir): array
    {
        $res = [];
        foreach ($annModel as $ann) {
            $file = $baseDir . '/' . $ann->name . '.' . $ann->type;
            if (!file_exists($ann->path)) {
                $this->warn('Unable to find annotation file for ' . $ann->name . '.');
                continue;
            }
            @copy($ann->path, $file);
            if (!file_exists($file)) {
                $this->warn('Unable to copy annotation file for ' . $ann->name . '.');
                continue;
            }
            $res[] = $file;
        }

        return $res;
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle()
    {
        $reference = $this->argument('reference');
        $annotations = $this->option('annotation');
        $outputDir = $this->option('output-dir') ?? getcwd();
        if (!is_writable($outputDir)) {
            $this->error('Output directory is not writable.');

            return 1;
        }
        $refModel = Reference::whereName($reference)->firstOrFail();
        $annModel = array_map(
            static function ($name) {
                return Annotation::whereName($name)->firstOrFail();
            },
            $annotations
        );
        $outputFile = $outputDir . '/' . $reference . '.zip';
        $baseDir = $refModel->basedir();
        $this->info('Creating config file...');
        $configFile = $this->makeConfigFile($refModel, $annModel, $baseDir);
        if (!file_exists($configFile)) {
            $this->error('Unable to write configuration file.');

            return 2;
        }
        $this->info('Copying annotations...');
        $annFiles = $this->copyAnnotations($annModel, $baseDir);
        $this->info('Building zip archive...');
        $zip = new ZipArchive();
        $zip->open($outputFile, ZipArchive::CREATE | ZipArchive::OVERWRITE);
        $files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($baseDir));
        foreach ($files as $name => $file) {
            if (!$file->isDir()) {
                $filePath = $file->getRealPath();
                $relativePath = substr($filePath, strlen($baseDir) + 1);
                $this->info('Adding ' . $relativePath . '...', 'v');
                $zip->addFile($filePath, $relativePath);
            }
        }
        $zip->close();
        $this->info('Cleaning up...');
        @unlink($configFile);
        foreach ($annFiles as $file) {
            @unlink($file);
        }
        $this->info('Completed! Results have been stored in ' . $outputFile);

        return 0;
    }
}
