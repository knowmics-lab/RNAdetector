<?php

declare(strict_types=1);

require __DIR__.'/vendor/autoload.php';

use Symplify\MonorepoBuilder\ComposerJsonManipulator\ValueObject\ComposerJsonSection;
use Symplify\MonorepoBuilder\Config\MBConfig;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\AddTagToChangelogReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\PushNextDevReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\PushTagReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\SetCurrentMutualDependenciesReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\SetNextMutualDependenciesReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\TagVersionReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\UpdateBranchAliasReleaseWorker;
use Symplify\MonorepoBuilder\Release\ReleaseWorker\UpdateReplaceReleaseWorker;
use Symplify\MonorepoBuilder\ValueObject\Option;
use RnaDetector\Monorepo\AddUnreleasedStubToChangelogFileWorker;
use RnaDetector\Monorepo\FillChangelogFileWorker;
use RnaDetector\Monorepo\UpdatePackageNextVersionWorker;
use RnaDetector\Monorepo\UpdatePackageVersionWorker;

return static function (MBConfig $mbConfig): void {
    $mbConfig->packageDirectories([__DIR__.'/packages']);
    $mbConfig->defaultBranch('master');
    $mbConfig->dataToRemove(
        [
            ComposerJsonSection::REQUIRE_DEV  => [
                'phpunit/phpunit'             => '*',
                'pestphp/pest'                => '*',
                'pestphp/pest-plugin-laravel' => '*',
            ],
            ComposerJsonSection::REPOSITORIES => [
                Option::REMOVE_COMPLETELY,
            ],
        ]
    );

    // release workers - in order to execute
    $mbConfig->workers(
        [
            FillChangelogFileWorker::class,
            UpdatePackageVersionWorker::class,
            UpdateReplaceReleaseWorker::class,
            SetCurrentMutualDependenciesReleaseWorker::class,
            AddTagToChangelogReleaseWorker::class,
            TagVersionReleaseWorker::class,
            PushTagReleaseWorker::class,
            SetNextMutualDependenciesReleaseWorker::class,
            UpdateBranchAliasReleaseWorker::class,
            UpdatePackageNextVersionWorker::class,
            AddUnreleasedStubToChangelogFileWorker::class,
            PushNextDevReleaseWorker::class,
        ]
    );
};
