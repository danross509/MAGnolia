#!/usr/bin/env nextflow

include { GTDBTK_CLASSIFYWF     } from '../../../modules/nf-core/gtdbtk/classifywf/main.nf'


workflow BIN_CLASSIFICATION {
    
    take:
        bins
        gtdb_db_dir
    
    main:

    if ( gtdb_db_dir.isDirectory() ) {
        gtdbtk_db = ["gtdbtk", gtdb_db_dir]
    } else {
        exit 1 ("gtdb_db must be supplied as a directory")
    }

    GTDBTK_CLASSIFYWF (
        bins,
        gtdbtk_db,
        params.gtdbtk_use_tmp_pplacer_dir
    )


    //emit:

}