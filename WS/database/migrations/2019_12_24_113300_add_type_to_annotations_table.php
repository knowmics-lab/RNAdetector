<?php

use Illuminate\Database\Migrations\Migration;
use Illuminate\Database\Schema\Blueprint;
use Illuminate\Support\Facades\Schema;

class AddTypeToAnnotationsTable extends Migration
{
    /**
     * Run the migrations.
     *
     * @return void
     */
    public function up()
    {
        Schema::table(
            'annotations',
            static function (Blueprint $table) {
                $table->enum('type', ['gtf', 'bed'])->after('name')->default('gtf');
            }
        );
    }

    /**
     * Reverse the migrations.
     *
     * @return void
     */
    public function down()
    {
        Schema::table(
            'annotations',
            static function (Blueprint $table) {
                $table->removeColumn('type');
            }
        );
    }
}
