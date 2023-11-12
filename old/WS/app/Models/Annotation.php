<?php
/**
 * RNADetector Web Service
 *
 * @author S. Alaimo, Ph.D. <alaimos at gmail dot com>
 */

namespace App\Models;

use Illuminate\Database\Eloquent\Model;

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
        return '/annotations/' . str_ireplace(config('rnadetector.annotations_path'), '', $this->getGFF3Path());
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
