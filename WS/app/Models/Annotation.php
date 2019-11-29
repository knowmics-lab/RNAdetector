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
 * @property string                          $path
 * @property \Illuminate\Support\Carbon|null $created_at
 * @property \Illuminate\Support\Carbon|null $updated_at
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation newModelQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation newQuery()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation query()
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereCreatedAt($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereId($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereName($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation wherePath($value)
 * @method static \Illuminate\Database\Eloquent\Builder|\App\Models\Annotation whereUpdatedAt($value)
 * @mixin \Eloquent
 */
class Annotation extends Model
{

    /**
     * The attributes that are mass assignable.
     *
     * @var array
     */
    protected $fillable = [
        'name',
        'path',
    ];

    /**
     * @inheritDoc
     */
    public function delete()
    {
        if (file_exists($this->path)) {
            unlink($this->path);
        }

        return parent::delete();
    }
}
