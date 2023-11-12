<?php

declare (strict_types=1);

namespace RnaDetector\Monorepo;

use PharIo\Version\Version;
use Symplify\MonorepoBuilder\ComposerJsonManipulator\FileSystem\JsonFileManager;
use Symplify\MonorepoBuilder\ComposerJsonManipulator\ValueObject\ComposerJson;
use Symplify\MonorepoBuilder\FileSystem\ComposerJsonProvider;
use Symplify\MonorepoBuilder\Release\Contract\ReleaseWorker\ReleaseWorkerInterface;

class UpdatePackageVersionWorker implements ReleaseWorkerInterface
{

    protected ComposerJsonProvider $composerJsonProvider;
    protected JsonFileManager $jsonFileManager;

    public function __construct(
        ComposerJsonProvider $composerJsonProvider,
        JsonFileManager $jsonFileManager
    ) {
        $this->composerJsonProvider = $composerJsonProvider;
        $this->jsonFileManager = $jsonFileManager;
    }

    public function getDescription(Version $version): string
    {
        return "Update the version in all composer.json files to the current version";
    }

    public function work(Version $version): void
    {
        $this->updateRootComposer($version);
        $this->updatePackageComposerJsons($version);
    }

    protected function updateRootComposer(Version $version): void
    {
        $rootComposerJson = $this->composerJsonProvider->getRootComposerJson();
        $this->updateComposerJson($rootComposerJson, $version);
    }

    protected function updatePackageComposerJsons(Version $version): void
    {
        $packageComposerJsons = $this->composerJsonProvider->getPackageComposerJsons();
        foreach ($packageComposerJsons as $packageComposerJson) {
            $this->updateComposerJson($packageComposerJson, $version);
        }
    }

    protected function updateComposerJson(ComposerJson $composerJson, Version $version): void
    {
        $composerJsonFileInfo = $composerJson->getFileInfo();
        if ($composerJsonFileInfo === null) {
            return;
        }
        $composerJson->setVersion($version->getVersionString());
        $this->jsonFileManager->printJsonToFileInfo($composerJson->getJsonArray(), $composerJsonFileInfo);
    }
}