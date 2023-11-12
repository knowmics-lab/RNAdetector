<?php

declare (strict_types=1);

namespace RnaDetector\Monorepo;

use MonorepoBuilderPrefix202310\Symplify\PackageBuilder\Parameter\ParameterProvider;
use MonorepoBuilderPrefix202310\Symplify\SmartFileSystem\SmartFileSystem;
use PharIo\Version\Version;
use Symplify\MonorepoBuilder\FileSystem\ComposerJsonProvider;
use Symplify\MonorepoBuilder\Release\Contract\ReleaseWorker\ReleaseWorkerInterface;
use Symplify\MonorepoBuilder\Release\Process\ProcessRunner;

use function array_filter;
use function array_map;
use function array_slice;
use function count;
use function explode;
use function file_exists;
use function getcwd;
use function sprintf;
use function str_replace;
use function trim;

final class FillChangelogFileWorker implements ReleaseWorkerInterface
{

    private const OPEN_COMMIT_SUBJECT    = 'open %s';
    private const UNRELEASED_PLACEHOLDER = '<!-- automatic release commit placeholder == DO NOT REMOVE == -->';
    private const GIT_LOG_COMMAND        = "git log --pretty=format:'Commit: %h%nSubject: %s%n<==>'";

    private readonly string $packageAliasFormat;

    public function __construct(
        private readonly ProcessRunner $processRunner,
        private readonly SmartFileSystem $smartFileSystem,
        private readonly ComposerJsonProvider $composerJsonProvider,
        ParameterProvider $parameterProvider
    ) {
        $this->packageAliasFormat = $parameterProvider->provideStringParameter('package_alias_format');
    }

    public function getDescription(Version $version): string
    {
        return "Fill the CHANGELOG.md file with the changes since the last release";
    }

    public function work(Version $version): void
    {
        $changelogFilePath = getcwd().'/CHANGELOG.md';
        if (!file_exists($changelogFilePath)) {
            return;
        }
        $gitHistory = $this->getGitCommitHistory($version);
        $gitHistoryString = $this->getGitHistoryString($gitHistory);
        $changelogFileContent = $this->smartFileSystem->readFile($changelogFilePath);
        $changelogFileContent = str_replace(self::UNRELEASED_PLACEHOLDER, $gitHistoryString, $changelogFileContent);
        $this->smartFileSystem->dumpFile($changelogFilePath, $changelogFileContent);
    }

    protected function getGitCommitHistory(Version $version): array
    {
        $logContent = $this->processRunner->run(self::GIT_LOG_COMMAND);
        $logData = array_filter(array_map("trim", explode("<==>", $logContent)));
        $logData = array_map(
            static function ($logEntry) {
                $logEntry = explode("\n", $logEntry);
                $logEntry = array_map(
                    static function ($logEntryLine) {
                        $logEntryData = explode(":", $logEntryLine, 2);

                        return trim($logEntryData[1]);
                    },
                    $logEntry
                );

                return [
                    "commit"  => $logEntry[0],
                    "subject" => $logEntry[1],
                ];
            },
            $logData
        );
        $lastCommitSubject = sprintf(self::OPEN_COMMIT_SUBJECT, $this->getVersionString($version));
        foreach ($logData as $key => $logEntry) {
            if ($logEntry["subject"] === $lastCommitSubject) {
                $logData = array_slice($logData, 0, $key);
                break;
            }
        }

        return $logData;
    }

    private function getGitHistoryString(array $gitHistory): string
    {
        $gitHistoryString = "";
        for ($i = count($gitHistory) - 1; $i >= 0; $i--) {
            $gitHistoryString .= sprintf("- %s [%s]\n", $gitHistory[$i]["subject"], $gitHistory[$i]["commit"]);
        }

        return trim($gitHistoryString);
    }

    private function getVersionString(Version $version): string
    {
        $rootComposerJson = $this->composerJsonProvider->getRootComposerJson();
        $rootVersion = new Version($rootComposerJson->getVersion());
        if ($rootVersion->isGreaterThan($version)) {
            return $this->getAliasFormat($rootVersion);
        }

        return $this->getAliasFormat($version);
    }

    private function getAliasFormat(Version $version): string
    {
        return str_replace(
            ['<major>', '<minor>'],
            [$version->getMajor()->getValue(), $version->getMinor()->getValue()],
            $this->packageAliasFormat
        );
    }
}