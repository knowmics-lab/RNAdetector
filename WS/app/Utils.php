<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App;

use App\Exceptions\CommandException;
use App\Exceptions\IgnoredException;
use App\Exceptions\ProcessingJobException;
use RecursiveDirectoryIterator;
use RecursiveIteratorIterator;
use SplFileInfo;
use Symfony\Component\Process\Exception\ProcessFailedException;
use Symfony\Component\Process\Process;
use Throwable;
use ZipArchive;

final class Utils
{

    public const IGNORED_ERROR_CODE = '===IGNORED===';

    /**
     * Runs a shell command and checks for successful completion of execution
     *
     * @param array         $command
     * @param string|null   $cwd
     * @param int|null      $timeout
     * @param callable|null $callback
     *
     * @return string|null
     */
    public static function runCommand(
        array $command,
        ?string $cwd = null,
        ?int $timeout = null,
        ?callable $callback = null
    ): ?string {
        $process = new Process($command, $cwd, null, null, $timeout);
        $process->run($callback);
        if (!$process->isSuccessful()) {
            throw new ProcessFailedException($process);
        }

        return $process->getOutput();
    }

    /**
     * Map command exception to message
     *
     * @param \Symfony\Component\Process\Exception\ProcessFailedException $e
     * @param array                                                       $errorCodeMap
     *
     * @return \App\Exceptions\ProcessingJobException|\App\Exceptions\IgnoredException
     */
    public static function mapCommandException(
        ProcessFailedException $e,
        array $errorCodeMap = []
    ) {
        $code = $e->getProcess()->getExitCode();
        if (isset($errorCodeMap[$code])) {
            if ($errorCodeMap[$code] === self::IGNORED_ERROR_CODE) {
                return new IgnoredException($code, $code);
            }

            return new ProcessingJobException($errorCodeMap[$code], $code, $e);
        }

        return new ProcessingJobException($e->getMessage(), $code, $e);
    }

    /**
     * Build a zip archive from a folder
     *
     * @param string $inputFolder
     * @param string $zipArchive
     *
     * @return bool
     */
    public static function makeZipArchive(string $inputFolder, string $zipArchive): bool
    {
        $rootPath = realpath($inputFolder);
        if (!file_exists($rootPath) && !is_dir($rootPath)) {
            return false;
        }
        $zip = new ZipArchive();
        $zip->open($zipArchive, ZipArchive::CREATE | ZipArchive::OVERWRITE);
        /** @var SplFileInfo[] $files */
        $files = new RecursiveIteratorIterator(new RecursiveDirectoryIterator($rootPath), RecursiveIteratorIterator::LEAVES_ONLY);
        foreach ($files as $name => $file) {
            if (!$file->isDir()) {
                $filePath = $file->getRealPath();
                $relativePath = substr($filePath, strlen($rootPath) + 1);
                $zip->addFile($filePath, $relativePath);
            }
        }
        $zip->close();

        return true;
    }

}
