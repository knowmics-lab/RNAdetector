<?php


namespace App;

use App\Exceptions\CommandException;
use App\Exceptions\ProcessingJobException;

final class Utils
{

    /**
     * Runs a shell command and checks for successful completion of execution
     *
     * @param string     $command
     * @param array|null $output
     *
     * @return boolean
     */
    public static function runCommand(string $command, array &$output = null): bool
    {
        $returnCode = -1;
        exec($command, $output, $returnCode);
        if ($returnCode !== 0) {
            throw new CommandException($returnCode);
        }
        return true;
    }

    /**
     * Map command exception to message
     *
     * @param string           $command
     * @param CommandException $e
     * @param array            $errorCodeMap
     *
     * @return \App\Exceptions\ProcessingJobException
     */
    public static function mapCommandException(string $command, CommandException $e, array $errorCodeMap = []): ProcessingJobException
    {
        $code = (int)$e->getMessage();
        if (isset($errorCodeMap[$code])) {
            return new ProcessingJobException($errorCodeMap[$code]);
        }
        return new ProcessingJobException('Execution of command "' . $command . '" returned error code ' . $code . '.');
    }

}
