<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

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
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference newModelQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference newQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference query()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereAvailableFor($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereCreatedAt($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereId($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereName($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference wherePath($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereMapPath($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Reference whereUpdatedAt($value)
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
     * @inheritDoc
     */
    public function delete()
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
