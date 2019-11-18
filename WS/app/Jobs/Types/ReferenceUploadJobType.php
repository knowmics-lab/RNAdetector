<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Jobs\Types;


use App\Exceptions\ProcessingJobException;
use App\Models\Annotation;
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
            ],
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
            'index.bwa'    => ['required', 'boolean'],
            'index.tophat' => ['required', 'boolean'],
            'index.salmon' => ['required', 'boolean'],
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
        $disk = Storage::disk('public');
        $dir = $this->model->getJobDirectory() . '/';
        if (!$disk->exists($dir . $fastaFile)) {
            return false;
        }

        return true;
    }

    /**
     * @param string $referenceFilename
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexBWA(string $referenceFilename): void
    {
        $this->log('Indexing reference for bwa.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('bwa_index.sh'),
                '-f',
                $referenceFilename,
                '-p',
                basename($referenceFilename),
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Input file does not exist.',
                4 => 'Output prefix must be specified',
                5 => 'Output directory is not writable',
            ]
        );
        $this->log($output);
        $this->log('Reference sequence indexed for bwa.');
    }

    /**
     * @param string $referenceFilename
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexTopHat(string $referenceFilename): void
    {
        $this->log('Indexing reference for TopHat.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('bowtie2_index.sh'),
                '-f',
                $referenceFilename,
                '-p',
                dirname($referenceFilename),
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'Input file does not exist.',
                4 => 'Output prefix must be specified.',
                5 => 'Output directory is not writable',
            ]
        );
        $this->log($output);
        $this->log('Reference sequence indexed for TopHat.');
    }

    /**
     * @param string $referenceFilename
     *
     * @throws \App\Exceptions\ProcessingJobException
     */
    private function indexSalmon(string $referenceFilename): void
    {
        $this->log('Indexing reference for Salmon.');
        $output = self::runCommand(
            [
                'bash',
                self::scriptPath('salmon_index_2.sh'),
                '-r',
                $referenceFilename,
                '-i',
                dirname($referenceFilename),
            ],
            $this->model->getAbsoluteJobDirectory(),
            null,
            null,
            [
                3 => 'FASTA file with transcripts does not exist.',
                4 => 'Indexed trascriptome folder does not exist.',
            ]
        );
        $this->log($output);
        $this->log('Reference sequence indexed for Salmon.');
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
        $absoluteSourceFilename = $this->model->getAbsoluteJobDirectory() . '/' . $file;
        $referenceFilename = env('REFERENCES_PATH') . '/' . $name . '.fasta';
        rename($absoluteSourceFilename, $referenceFilename);
        if (!file_exists($referenceFilename)) {
            throw new ProcessingJobException('Unable to create source fasta file.');
        }
        $bwa = (bool)($index['bwa'] ?? false);
        $tophat = (bool)($index['tophat'] ?? false);
        $salmon = (bool)($index['salmon'] ?? false);
        if ($bwa) {
            $this->indexBWA($referenceFilename);
        }
        if ($tophat) {
            $this->indexTopHat($referenceFilename);
        }
        if ($salmon) {
            $this->indexSalmon($referenceFilename);
        }
        Reference::create(
            [
                'name'          => $name,
                'path'          => $referenceFilename,
                'available_for' => [
                    'bwa'    => $bwa,
                    'tophat' => $tophat,
                    'salmon' => $salmon,
                ],
            ]
        )->save();
        $this->log('Job completed!');
        $this->model->setOutput(['done' => true]);
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
