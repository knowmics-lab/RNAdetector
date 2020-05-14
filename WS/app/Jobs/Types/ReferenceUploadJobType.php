<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Models\Reference;
use Illuminate\Http\Request;
use Storage;

class ReferenceUploadJobType extends AbstractJob
{

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'name'      => 'A name for this reference sequence',
            'fastaFile' => 'The fasta file for this reference sequence',
            'index'     => [
                'bwa'    => 'A boolean for enabling bwa indexing',
                'tophat' => 'A boolean for enabling tophat indexing',
                'salmon' => 'A boolean for enabling salmon indexing',
                'hisat'  => 'A boolean for enabling hisat2 indexing',
            ],
            'map_file'  => 'An optional map file (tab-separated with two columns) where each line contains the ID of a transcript/gene and its Entrez Gene Id (Mirbase mature name for miRNAs).',
        ];
    }

    /**
     * Returns an array containing for each output value an help detailing its use.
     *
     * @return array
     */
    public static function outputSpec(): array
    {
        return [
            'done' => 'A boolean that is true if the annotation has been successfully created',
        ];
    }

    /**
     * Returns an array containing rules for input validation.
     *
     * @param \Illuminate\Http\Request $request
     *
     * @return array
     */
    public static function validationSpec(Request $request): array
    {
        return [
            'name'         => ['required', 'alpha_dash', 'max:255'],
            'fastaFile'    => ['required', 'string'],
            'index'        => ['required', 'array'],
            'index.bwa'    => ['filled', 'boolean'],
            'index.tophat' => ['filled', 'boolean'],
            'index.salmon' => ['filled', 'boolean'],
            'index.hisat'  => ['filled', 'boolean'],
            'map_file'     => ['filled', 'string'],
        ];
    }

    /**
     * Checks the input of this job and returns true iff the input contains valid data
     * The default implementation does nothing.
     *
     * @return bool
     */
    public function isInputValid(): bool
    {
        $fastaFile = $this->model->getParameter('fastaFile');
        $mapFile = $this->model->getParameter('map_file');
        $disk = Storage::disk('public');
        $dir = $this->model->getJobDirectory() . '/';
        if (!$disk->exists($dir . $fastaFile)) {
            return false;
        }
        if ($mapFile && !$disk->exists($dir . $mapFile)) {
            return false;
        }

        return true;
    }

    /**
     * @param string $referenceFilename
     * @param string $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexBWA(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for bwa.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('bwa_index.sh'),
                '-f',
                $referenceFilename,
                '-p',
                $referenceDirname . '/reference',
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                3 => 'Input file does not exist.',
                4 => 'Output prefix must be specified',
                5 => 'Output directory is not writable',
            ]
        );
        // $this->log($output);
        $this->log('Reference sequence indexed for bwa.');
    }

    /**
     * @param string $referenceFilename
     * @param string $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexTopHat(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for TopHat.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('bowtie2_index.sh'),
                '-f',
                $referenceFilename,
                '-p',
                $referenceDirname . '/reference',
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                3 => 'Input file does not exist.',
                4 => 'Output prefix must be specified.',
                5 => 'Output directory is not writable',
            ]
        );
        // $this->log($output);
        $this->log('Reference sequence indexed for TopHat.');
    }

    /**
     * @param string $referenceFilename
     * @param string $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexSalmon(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for Salmon.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('salmon_index_2.sh'),
                '-r',
                $referenceFilename,
                '-i',
                $referenceDirname . '/reference',
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                3 => 'FASTA file with transcripts does not exist.',
                4 => 'Indexed trascriptome folder does not exist.',
            ]
        );
        // $this->log($output);
        $this->log('Reference sequence indexed for Salmon.');
    }

    /**
     * @param string $referenceFilename
     * @param string $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexHisat(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for Hisat 2.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('hisat_index.sh'),
                '-f',
                $referenceFilename,
                '-p',
                $referenceDirname . '/reference',
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                3 => 'Input file does not exist.',
                4 => 'Output prefix must be specified.',
                5 => 'Output directory must be writeable.',
                6 => 'An unknown error occurred while running hisat2-build.',
            ]
        );
        // $this->log($output);
        $this->log('Reference sequence indexed for Hisat 2.');
    }


    /**
     * Handles all the computation for this job.
     * This function should throw a ProcessingJobException if something went wrong during the computation.
     * If no exceptions are thrown the job is considered as successfully completed.
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    public function handle(): void
    {
        $this->log('Starting job.');
        $name = $this->model->getParameter('name');
        $file = $this->model->getParameter('fastaFile');
        $index = (array)$this->model->getParameter('index', []);
        $mapFile = $this->model->getParameter('map_file');
        $absoluteSourceFilename = $this->model->getAbsoluteJobDirectory() . '/' . $file;
        $referenceDirname = env('REFERENCES_PATH') . '/' . $name;
        $referenceFilename = $referenceDirname . '/reference.fa';
        if (!mkdir($referenceDirname, 0777, true) && !is_dir($referenceDirname)) {
            throw new ProcessingJobException(sprintf('Directory "%s" was not created', $referenceDirname));
        }
        @chmod($referenceDirname, 0777);
        if (!file_exists($referenceDirname)) {
            throw new ProcessingJobException('Unable to create reference directory.');
        }
        rename($absoluteSourceFilename, $referenceFilename);
        @chmod($referenceFilename, 0777);
        if (!file_exists($referenceFilename)) {
            throw new ProcessingJobException('Unable to create source fasta file.');
        }
        $bwa = (bool)($index['bwa'] ?? false);
        $tophat = (bool)($index['tophat'] ?? false);
        $salmon = (bool)($index['salmon'] ?? false);
        $hisat = (bool)($index['hisat'] ?? false);
        if ($bwa) {
            $this->indexBWA($referenceFilename, $referenceDirname);
        }
        if ($tophat) {
            $this->indexTopHat($referenceFilename, $referenceDirname);
        }
        if ($salmon) {
            $this->indexSalmon($referenceFilename, $referenceDirname);
        }
        if ($hisat) {
            $this->indexHisat($referenceFilename, $referenceDirname);
        }
        $mapFileName = null;
        if ($mapFile) {
            $mapFileName = $referenceDirname . '/map_file.tsv';
            $this->moveFile($mapFile, $mapFileName);
            if (!file_exists($mapFileName)) {
                throw new ProcessingJobException('Unable to create map file.');
            }
        }
        Reference::create(
            [
                'name'          => $name,
                'path'          => $referenceFilename,
                'available_for' => [
                    'bwa'    => $bwa,
                    'tophat' => $tophat,
                    'salmon' => $salmon,
                    'hisat'  => $hisat,
                ],
                'map_path'      => $mapFileName,
            ]
        )->save();
        $this->log('Job completed!');
        $this->model->setOutput(['type' => self::OUT_TYPE_CONFIRMATION, 'done' => true]);
        $this->model->save();
    }

    /**
     * Returns a description for this job
     *
     * @return string
     */
    public static function description(): string
    {
        return 'Upload a novel genome annotation';
    }
}
