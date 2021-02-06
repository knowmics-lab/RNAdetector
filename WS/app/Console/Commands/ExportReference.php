<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Console\Commands;

use App\Jobs\Types\AbstractJob;
use App\Models\Annotation;
use App\Models\Reference;
use App\Utils;
use Illuminate\Console\Command;
use Illuminate\Support\Str;
use Pheanstalk\Command\AbstractCommand;
use Symfony\Component\Process\Process;

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

    private $jsonTemplate = <<<TEMPLATE
{
  "name": "%s",
  "title": "%s",
  "description": "{DESCRIPTION}",
  "url": "{URL}/%s",
  "md5": "{URL}/%s",
  "version": "%s"
}
TEMPLATE;


    /**
     * Make config file
     *
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
            'noReference' => false,
            'indexedFor'  => $refModel->available_for,
            'mapFile'     => ($refModel->map_path !== null),
            'annotations' => array_map(
                static function (Annotation $model) {
                    return [
                        'name'    => $model->name,
                        'type'    => $model->type,
                        'mapFile' => $model->map_path !== null,
                    ];
                },
                $annModel
            ),
        ];
        file_put_contents($configFile, json_encode($configContent));

        return $configFile;
    }

    /**
     * Copy annotations to the reference folder
     *
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
            if ($ann->map_path) {
                if (!file_exists($ann->map_path)) {
                    $this->warn('Unable to find map file for ' . $ann->name . '.');
                    continue;
                }
                $mapFile = $baseDir . '/' . $ann->name . '_map.tsv';
                @copy($ann->map_path, $mapFile);
                if (!file_exists($mapFile)) {
                    $this->warn('Unable to copy map file for ' . $ann->name . '.');
                    continue;
                }
                $res[] = $mapFile;
            }
            if ($ann->hasGFF3()) {
                $gff3File = $baseDir . '/' . $ann->name . '.gff3.gz';
                @copy($ann->getGFF3Path(), $gff3File);
                if (!file_exists($gff3File)) {
                    $this->warn('Unable to copy GFF3 file for ' . $ann->name . '.');
                    continue;
                }
                $res[] = $gff3File;
            }
        }

        return $res;
    }

    /**
     * Execute the console command.
     *
     * @return mixed
     */
    public function handle(): int
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
        $outputFile = $outputDir . '/' . $reference . '.tar.bz2';
        $baseDir = realpath($refModel->basedir());
        $installedFile = $baseDir . '/.installed';
        $versionFile = $baseDir . '/.version';
        $hasInstalled = file_exists($installedFile);
        $hasVersion = file_exists($versionFile);
        $version = ($hasVersion) ? file_get_contents($versionFile) : '';
        if ($hasInstalled) {
            @unlink($installedFile);
        }
        if ($hasVersion) {
            @unlink($versionFile);
        }
        $this->info('Creating config file...');
        $configFile = $this->makeConfigFile($refModel, $annModel, $baseDir);
        if (!file_exists($configFile)) {
            $this->error('Unable to write configuration file.');

            return 2;
        }
        $this->info('Copying annotations...');
        $annFiles = $this->copyAnnotations($annModel, $baseDir);
        $this->info('Building archive...');
        Utils::runCommand(
            [
                'tar',
                '-jcvf',
                $outputFile,
                basename($baseDir),
            ],
            dirname($baseDir),
            null,
            function ($type, $buffer) {
                if ($type === Process::OUT) {
                    $this->output->write('<info>' . $buffer . '</info>');
                }
            }
        );
        $this->info('Cleaning up...');
        @unlink($configFile);
        foreach ($annFiles as $file) {
            @unlink($file);
        }
        if ($hasInstalled) {
            @touch($installedFile);
            @chmod($installedFile, 0777);
        }
        if ($hasVersion) {
            @file_put_contents($versionFile, $version);
            @chmod($versionFile, 0777);
        }
        $this->info('Computing MD5 checksum...');
        Utils::runCommand(
            [
                'bash',
                AbstractJob::scriptPath('compute_md5.sh'),
                '-i',
                $outputFile,
            ],
            dirname($baseDir),
            null,
            function ($type, $buffer) {
                if ($type === Process::ERR) {
                    $this->error($buffer);
                }
            }
        );
        $md5File = $outputFile . ".md5";
        if (!file_exists($md5File)) {
            $this->error('MD5 checksum does not exist!');

            return 3;
        }
        $hash = preg_split("/\\s+/", file_get_contents($md5File))[0];
        if (empty($hash)) {
            $this->error('MD5 checksum does not exist!');

            return 4;
        }
        $jsonFile = $outputFile . ".json";
        file_put_contents(
            $jsonFile,
            sprintf(
                $this->jsonTemplate,
                $reference,
                str_replace(['-', '_'], ' ', $reference),
                basename($outputFile),
                basename($md5File),
                $hash
            )
        );
        $this->info('Completed! Results have been stored in ' . $outputDir);

        return 0;
    }
}
