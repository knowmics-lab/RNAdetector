<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

use Illuminate\Database\Eloquent\Builder;
use Illuminate\Database\Eloquent\Model;
use RecursiveDirectoryIterator;
use RecursiveIteratorIterator;

/**
 * App\Models\Reference
 *
 * @property int                             $id
 * @property string                          $name
 * @property string                          $path
 * @property string|null                     $map_path
 * @property array                           $available_for
 * @property \Illuminate\Support\Carbon|null $created_at
 * @property \Illuminate\Support\Carbon|null $updated_at
 * @method static Builder|Reference newModelQuery()
 * @method static Builder|Reference newQuery()
 * @method static Builder|Reference query()
 * @method static Builder|Reference whereAvailableFor($value)
 * @method static Builder|Reference whereCreatedAt($value)
 * @method static Builder|Reference whereId($value)
 * @method static Builder|Reference whereName($value)
 * @method static Builder|Reference wherePath($value)
 * @method static Builder|Reference whereMapPath($value)
 * @method static Builder|Reference whereUpdatedAt($value)
 * @mixin \Eloquent
 */
class Reference extends Model
{

    /**
     * The attributes that are mass assignable.
     *
     * @var array
     */
    protected $fillable = [
        'name',
        'path',
        'map_path',
        'available_for',
    ];

    /**
     * The attributes that should be cast to native types.
     *
     * @var array
     */
    protected $casts = [
        'available_for' => 'array',
    ];

    /**
     * Delete a path
     *
     * @param $path
     *
     * @return void
     */
    private static function _deletePath($path): void
    {
        if (is_file($path)) {
            unlink($path);
        } elseif (is_dir($path)) {
            $files = new RecursiveIteratorIterator(
                new RecursiveDirectoryIterator($path, RecursiveDirectoryIterator::SKIP_DOTS),
                RecursiveIteratorIterator::CHILD_FIRST
            );

            foreach ($files as $fileinfo) {
                $todo = ($fileinfo->isDir() ? 'rmdir' : 'unlink');
                $todo($fileinfo->getRealPath());
            }

            rmdir($path);
        }
    }

    /**
     * Returns the base directory for this reference sequence
     *
     * @return string
     */
    public function basedir(): string
    {
        return dirname($this->path);
    }

    /**
     * Returns the basename for this reference sequence
     *
     * @return string
     */
    public function basename(): string
    {
        return $this->basedir() . '/reference';
    }

    /**
     * Checks if this reference sequence is available for a specific algorithm
     *
     * @param string $algorithm
     *
     * @return bool
     */
    public function isAvailableFor(string $algorithm): bool
    {
        return $this->available_for[$algorithm] ?? false;
    }

    /**
     * Checks if this reference sequence FASTA file is indexed
     *
     * @return bool
     */
    public function isFastaIndexed(): bool
    {
        return file_exists($this->path . '.fai');
    }

    /**
     * Returns the URI of the FASTA file of this reference
     *
     * @return string
     */
    public function getReferenceUri(): string
    {
        return '/references/' . str_ireplace(config('rnadetector.reference_path'), '', $this->path);
    }

    /**
     * @inheritDoc
     */
    public function delete(): ?bool
    {
        if (file_exists($this->basedir())) {
            self::_deletePath($this->basedir());
        }
        if ($this->map_path && file_exists($this->map_path)) {
            @unlink($this->map_path);
        }

        return parent::delete();
    }


}
