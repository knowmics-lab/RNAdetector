<?php

declare (strict_types=1);

namespace RnaDetector\Monorepo;

use MonorepoBuilderPrefix202310\Symplify\SmartFileSystem\SmartFileSystem;
use PharIo\Version\Version;
use Symplify\MonorepoBuilder\Release\Contract\ReleaseWorker\ReleaseWorkerInterface;

use function explode;
use function file_exists;
use function getcwd;
use function implode;

final class AddUnreleasedStubToChangelogFileWorker implements ReleaseWorkerInterface
{

    private const UNRELEASED_STUB = "\n## Unreleased\n\n<!-- automatic release commit placeholder == DO NOT REMOVE == -->\n";

    public function __construct(private readonly SmartFileSystem $smartFileSystem) {}

    public function getDescription(Version $version): string
    {
        return "Add the \"Unreleased\" section stub to the CHANGELOG.md file";
    }

    public function work(Version $version): void
    {
        $changelogFilePath = getcwd().'/CHANGELOG.md';
        if (!file_exists($changelogFilePath)) {
            return;
        }
        $changelogFileContent = $this->smartFileSystem->readFile($changelogFilePath);
        $sections = explode("\n##", $changelogFileContent, 2);
        $sections[0] .= self::UNRELEASED_STUB;
        $changelogFileContent = implode("\n##", $sections);
        $this->smartFileSystem->dumpFile($changelogFilePath, $changelogFileContent);
    }

}