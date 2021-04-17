<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Jobs\Types\Traits\HandlesCompressedFastqTrait;
use App\Models\Reference;
use Illuminate\Http\Request;
use Storage;

class ReferenceUploadJobType extends AbstractJob
{

    use HandlesCompressedFastqTrait;

    /**
     * Returns an array containing for each input parameter an help detailing its content and use.
     *
     * @return array
     */
    public static function parametersSpec(): array
    {
        return [
            'name'             => 'A name for this reference sequence',
            'fastaFile'        => 'The fasta file for this reference sequence',
            'index'            => [
                'bwa'    => 'A boolean for enabling bwa indexing',
                'salmon' => 'A boolean for enabling salmon indexing',
                'hisat'  => 'A boolean for enabling hisat2 indexing',
                'star'   => 'A boolean for enabling STAR indexing',
            ],
            'custom_arguments' => [
                'bwa'    => 'An optional string containing custom arguments for bwa command line call',
                'salmon' => 'An optional string containing custom arguments for salmon command line call',
                'hisat'  => 'An optional string containing custom arguments for HISAT2 command line call',
                'star'   => 'An optional string containing custom arguments for STAR command line call',
            ],
            'map_file'         => 'An optional map file (tab-separated with two columns) where each line contains the ID of a transcript/gene and its Entrez Gene Id (Mirbase mature name for miRNAs).',
        ];
    }

    /**
     * @inheritDoc
     */
    public function threads(): int
    {
        return 1;
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
     * @param  \Illuminate\Http\Request  $request
     *
     * @return array
     */
    public static function validationSpec(Request $request): array
    {
        return [
            'name'                    => ['required', 'alpha_dash', 'max:255'],
            'fastaFile'               => ['required', 'string'],
            'index'                   => ['required', 'array'],
            'index.bwa'               => ['filled', 'boolean'],
            'index.salmon'            => ['filled', 'boolean'],
            'index.hisat'             => ['filled', 'boolean'],
            'index.star'              => ['filled', 'boolean'],
            'custom_arguments'        => ['filled', 'array'],
            'custom_arguments.bwa'    => ['filled', 'string'],
            'custom_arguments.salmon' => ['filled', 'string'],
            'custom_arguments.hisat'  => ['filled', 'string'],
            'custom_arguments.star'   => ['filled', 'string'],
            'map_file'                => ['nullable', 'string'],
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
     * @param  string  $referenceFilename
     * @param  string  $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexBWA(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for bwa.');
        self::runCommand(
            $this->appendCustomArguments(
                [
                    'bash',
                    self::scriptPath('bwa_index.sh'),
                    '-f',
                    $referenceFilename,
                    '-p',
                    $referenceDirname . '/reference',
                ],
                'custom_arguments.bwa'
            ),
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
     * @param  string  $referenceFilename
     * @param  string  $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexSalmon(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for Salmon.');
        self::runCommand(
            $this->appendCustomArguments(
                [
                    'bash',
                    self::scriptPath('salmon_index_2.sh'),
                    '-r',
                    $referenceFilename,
                    '-i',
                    $referenceDirname . '/reference',
                ],
                'custom_arguments.salmon'
            ),
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
     * @param  string  $referenceFilename
     * @param  string  $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexHisat(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for Hisat 2.');
        self::runCommand(
            $this->appendCustomArguments(
                [
                    'bash',
                    self::scriptPath('hisat_index.sh'),
                    '-f',
                    $referenceFilename,
                    '-p',
                    $referenceDirname . '/reference',
                ],
                'custom_arguments.hisat'
            ),
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
     * @param  string  $referenceFilename
     * @param  string  $referenceDirname
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexStar(string $referenceFilename, string $referenceDirname): void
    {
        $this->log('Indexing reference for STAR.');
        self::runCommand(
            $this->appendCustomArguments(
                [
                    'bash',
                    self::scriptPath('star_index.sh'),
                    '-f',
                    $referenceFilename,
                    '-p',
                    $referenceDirname . '/reference',
                ],
                'custom_arguments.star'
            ),
            $this->model->getAbsoluteJobDirectory(),
            null,
            function ($type, $buffer) {
                $this->log($buffer, false);
            },
            [
                3 => 'Input file does not exist.',
                4 => 'Output directory must be specified.',
                5 => 'Output directory must be writeable.',
                6 => 'An unknown error occurred while indexing.',
            ]
        );
        // $this->log($output);
        $this->log('Reference sequence indexed for STAR.');
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
        $absoluteSourceFilename = $this->model->getAbsoluteJobDirectory() . '/' . self::checksForCompression(
                $this->model,
                $absoluteSourceFilename
            );
        if ($absoluteSourceFilename === null) {
            throw new ProcessingJobException("An unknown error occurred!");
        }
        $referenceDirname = config('rnadetector.reference_path') . '/' . $name;
        $referenceFilename = $referenceDirname . '/reference.fa';
        if (!file_exists($referenceDirname) && !mkdir($referenceDirname, 0777, true) && !is_dir($referenceDirname)) {
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
        $salmon = (bool)($index['salmon'] ?? false);
        $hisat = (bool)($index['hisat'] ?? false);
        $star = (bool)($index['star'] ?? false);
        if ($bwa) {
            $this->indexBWA($referenceFilename, $referenceDirname);
        }
        if ($salmon) {
            $this->indexSalmon($referenceFilename, $referenceDirname);
        }
        if ($hisat) {
            $this->indexHisat($referenceFilename, $referenceDirname);
        }
        if ($star) {
            $this->indexStar($referenceFilename, $referenceDirname);
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
                    'salmon' => $salmon,
                    'hisat'  => $hisat,
                    'star'   => $star,
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

    /**
     * @inheritDoc
     */
    public static function displayName(): string
    {
        return 'Reference Upload';
    }
}
