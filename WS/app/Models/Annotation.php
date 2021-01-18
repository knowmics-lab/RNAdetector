<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

use Illuminate\Database\Eloquent\Model;

/**
 * App\Models\Annotation
 *
 * @property int                             $id
 * @property string                          $name
 * @property string                          $type
 * @property string                          $path
 * @property string|null                     $map_path
 * @property \Illuminate\Support\Carbon|null $created_at
 * @property \Illuminate\Support\Carbon|null $updated_at
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation newModelQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation newQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation query()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereCreatedAt($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereId($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereName($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereType($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation wherePath($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereMapPath($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereUpdatedAt($value)
 * @mixin \Eloquent
 */
class Annotation extends Model
{

    public const GTF = 'gtf';
    public const BED = 'bed';

    /**
     * The attributes that are mass assignable.
     *
     * @var array
     */
    protected $fillable = [
        'name',
        'type',
        'path',
        'map_path',
    ];

    /**
     * Checks if this annotation is in GTF format
     *
     * @return bool
     */
    public function isGtf(): bool
    {
        return $this->type === self::GTF;
    }

    /**
     * Checks if this annotation is in BED format
     *
     * @return bool
     */
    public function isBed(): bool
    {
        return $this->type === self::BED;
    }

    /**
     * Returns the path of the GFF3 file associated to this annotation
     * (A GFF3 file is used for viewing of GTF annotations in JBrowse2)
     * @return string
     */
    public function getGFF3Path(): string
    {
        return dirname($this->path) . '/' . basename($this->path, '.gtf') . '.gff3.gz';
    }

    /**
     * Returns the URI of the GFF3 file associated to this annotation
     *
     * @return string
     */
    public function getGFF3Uri(): string
    {
        return '/annotations/' . str_ireplace(env('ANNOTATIONS_PATH'), '', $this->getGFF3Path());
    }

    /**
     * Checks if this annotation has any GFF3 file associated
     *
     * @return bool
     */
    public function hasGFF3(): bool
    {
        if (!$this->isGtf()) {
            return false;
        }

        return file_exists($this->getGFF3Path());
    }

    /**
     * @inheritDoc
     */
    public function delete(): ?bool
    {
        if (file_exists($this->path)) {
            @unlink($this->path);
        }
        if ($this->map_path && file_exists($this->map_path)) {
            @unlink($this->map_path);
        }

        return parent::delete();
    }
}
