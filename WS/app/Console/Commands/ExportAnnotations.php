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
use Symfony\Component\Process\Process;

class ExportAnnotations extends Command
{

    /**
     * The name and signature of the console command.
     *
     * @var string
     */
    protected $signature = 'annotations:export
                            {name : The name of the final package}
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
     * @param Annotation[] $annModel
     * @param string       $baseDir
     *
     * @return string
     */
    private function makeConfigFile(array $annModel, string $baseDir): string
    {
        $configFile = $baseDir . '/spec.json';
        $configContent = [
            'noReference' => true,
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
        $name = $this->argument('name');
        $annotations = $this->option('annotation');
        $outputDir = $this->option('output-dir') ?? getcwd();
        if (!is_writable($outputDir)) {
            $this->error('Output directory is not writable.');

            return 1;
        }
        $referenceDirname = env('REFERENCES_PATH') . '/' . $name;
        if (!file_exists($referenceDirname) && !mkdir($referenceDirname, 0777) && !is_dir($referenceDirname)) {
            $this->warn(sprintf('Directory "%s" was not created', $referenceDirname));

            return 3;
        }
        $annModel = array_map(
            static function ($name) {
                return Annotation::whereName($name)->firstOrFail();
            },
            $annotations
        );
        $outputFile = $outputDir . '/' . $name . '.tar.bz2';
        $baseDir = realpath($referenceDirname);
        $installedFile = $baseDir . '/.installed';
        $hasInstalled = file_exists($installedFile);
        if ($hasInstalled) {
            @unlink($installedFile);
        }
        $this->info('Creating config file...');
        $configFile = $this->makeConfigFile($annModel, $baseDir);
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
        if (file_exists($referenceDirname)) {
            @rmdir($referenceDirname);
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
        $this->info('Completed! Results have been stored in ' . $outputFile);
        $this->info('JSON Template for the repository: ');
        $this->info(
            sprintf($this->jsonTemplate, $name, str_replace(['-', '_'], ' ', $name), basename($outputFile), basename($md5File), $hash)
        );

        return 0;
    }
}
