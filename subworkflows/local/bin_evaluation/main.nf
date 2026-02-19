/*
 * BUSCO/CheckM/CheckM2/GUNC: Quantitative measures for the assessment of genome assembly
 */

include { BUSCO_DB_PREPARATION             } from '../../../modules/nf-core_mag/busco/db_preparation/main.nf'
include { BUSCO                            } from '../../../modules/nf-core_mag/busco/busco/main.nf'
include { BUSCO_SAVE_DOWNLOAD              } from '../../../modules/nf-core_mag/busco/save_download/main.nf'
include { BUSCO_SUMMARY                    } from '../../../modules/nf-core_mag/busco/summary/main.nf'


include { CHECKM_QA                        } from '../../../modules/nf-core/checkm/qa/main.nf'
include { CHECKM_LINEAGEWF                 } from '../../../modules/nf-core/checkm/lineagewf/main.nf'
include { CHECKM2_PREDICT                  } from '../../../modules/nf-core/checkm2/predict/main.nf'


include { COMBINE_TSV as COMBINE_BINQC_TSV } from '../../../modules/nf-core_mag/combine_tsv/main.nf'


workflow BIN_EVALUATION {
    take:
    post_refinement_bins // [meta, [list, of, bin, fastas]]
    checkm2_db_dir // [dbmeta, db]
    checkm_db_dir // [db]S

    main:
    ch_versions = channel.empty()

    // MODIFY TO ACCOUNT FOR POSSIBILITY OF BOTH CHECKM AND CHECKM2

    evaluation_input = post_refinement_bins.transpose()
        .map { _meta, bin ->
            def meta_new = [:]
            meta_new.id = 'post_refinement_bins'
            [ meta_new, bin ]
        }
        .groupTuple( by: 0 )

    if ( params.checkm_version == "checkm2" ) {
        CHECKM2_PREDICT(evaluation_input, checkm2_db_dir)

        //COMBINE_BINQC_TSV(CHECKM2_PREDICT.out.checkm2_tsv.collect { summary -> summary[1] })

        bin_summary = CHECKM2_PREDICT.out.checkm2_tsv
        ch_versions = ch_versions.mix(
            CHECKM2_PREDICT.out.versions.first(),
          //  COMBINE_BINQC_TSV.out.versions
        )
    } else if ( params.checkm_version == "checkm" ) {
        bins_for_checkmlineagewf = evaluation_input
            .filter { meta, _bins ->
                meta.domain != "eukarya"
            }
            .multiMap { meta, fa ->
                reads: [meta, fa]
                ext: fa.extension.unique().join("")
            }

        CHECKM_LINEAGEWF(
            bins_for_checkmlineagewf.reads, 
            bins_for_checkmlineagewf.ext, 
            checkm_db_dir
        )
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        /*ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }*/

        //CHECKM_QA(ch_checkmqa_input, [])

        //COMBINE_BINQC_TSV(CHECKM_QA.out.output.collect { summary -> summary[1] })

        bin_summary = CHECKM_LINEAGEWF.out.checkm_tsv
        ch_versions = ch_versions.mix(
            CHECKM_QA.out.versions.first(),
            //COMBINE_BINQC_TSV.out.versions
        )
    }

    emit:
    bins_output     = evaluation_input
    bin_summary     = bin_summary
    //multiqc_files = ch_multiqc_files
    versions        = ch_versions
}