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
    bins // [meta, [list, of, bin, fastas]]
    checkm2_db_dir // [dbmeta, db]
    checkm_db_dir // [db]S

    main:
    ch_versions = channel.empty()

    /*if (params.binqc_tool == "busco") {
        if (!ch_busco_db.isEmpty()) {
            if (ch_busco_db.extension in ['gz', 'tgz']) {
                // Expects to be tar.gz!
                BUSCO_DB_PREPARATION(ch_busco_db)
                ch_db_for_busco = BUSCO_DB_PREPARATION.out.db.map { meta, db ->
                    [[id: meta, lineage: 'Y'], db]
                }
            }
            else if (ch_busco_db.isDirectory()) {
                // Set meta to match expected channel cardinality for BUSCO
                ch_db_for_busco = Channel
                    .of(ch_busco_db)
                    .collect { db ->
                        def basename = db.getBaseName()
                        def lineage = basename.contains('odb10') ? 'Y' : 'N'
                        [[id: basename, lineage: lineage], db]
                    }
            }
        }
        else {
            // Set BUSCO database to empty to allow for --auto-lineage
            ch_db_for_busco = Channel
                .of([[lineage: ''], []])
                .collect()
        }

        if (params.save_busco_db) {
            // publish files downloaded by Busco
            ch_downloads = BUSCO.out.busco_downloads
                .groupTuple()
                .map { _lin, downloads -> downloads[0] }
                .toSortedList()
                .flatten()
            BUSCO_SAVE_DOWNLOAD(ch_downloads)

            ch_versions = ch_versions.mix(BUSCO_SAVE_DOWNLOAD.out.versions.first())
        }

        BUSCO(ch_input_bins_for_qc, ch_db_for_busco)

        BUSCO_SUMMARY(
            BUSCO.out.summary_domain.collect { _meta, summary -> summary }.ifEmpty([]),
            BUSCO.out.summary_specific.collect { _meta, summary -> summary }.ifEmpty([]),
            BUSCO.out.failed_bin.collect { _meta, summary -> summary }.ifEmpty([])
        )

        ch_multiqc_files = ch_multiqc_files.mix(
            BUSCO.out.summary_domain.mix(BUSCO.out.summary_specific).map { _meta, summary -> summary }
        )
        qc_summary = BUSCO_SUMMARY.out.summary
        ch_versions = ch_versions.mix(BUSCO.out.versions.first())
    }
    */

    if (params.checkm_version == "checkm2") {
        CHECKM2_PREDICT(bins, checkm2_db_dir)

        COMBINE_BINQC_TSV(CHECKM2_PREDICT.out.checkm2_tsv.collect { summary -> summary[1] })

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(
            CHECKM2_PREDICT.out.versions.first(),
            COMBINE_BINQC_TSV.out.versions
        )
    }
    else if (params.checkm_version == "checkm") {
        bins_for_checkmlineagewf = bins
            .filter { meta, _bins ->
                meta.domain != "eukarya"
            }
            .multiMap { meta, fa ->
                reads: [meta, fa]
                ext: fa.extension.unique().join("")
            }

        CHECKM_LINEAGEWF(bins_for_checkmlineagewf.reads, bins_for_checkmlineagewf.ext, checkm_db_dir)
        ch_versions = ch_versions.mix(CHECKM_LINEAGEWF.out.versions.first())

        ch_checkmqa_input = CHECKM_LINEAGEWF.out.checkm_output
            .join(CHECKM_LINEAGEWF.out.marker_file)
            .map { meta, dir, marker ->
                [meta, dir, marker, []]
            }

        CHECKM_QA(ch_checkmqa_input, [])

        COMBINE_BINQC_TSV(CHECKM_QA.out.output.collect { summary -> summary[1] })

        qc_summary = COMBINE_BINQC_TSV.out.combined
        ch_versions = ch_versions.mix(
            CHECKM_QA.out.versions.first(),
            COMBINE_BINQC_TSV.out.versions
        )
    }

    /*if (params.run_gunc) {
        // GUNC
        ch_input_bins_for_gunc = bins
            .filter { meta, _bins ->
                meta.domain != "eukarya"
            }

        if (params.gunc_db) {
            ch_db_for_gunc = ch_gunc_db
        }
        else {
            ch_db_for_gunc = GUNC_DOWNLOADDB(params.gunc_database_type).db
            ch_versions.mix(GUNC_DOWNLOADDB.out.versions)
        }

        GUNC_RUN(ch_input_bins_for_gunc, ch_db_for_gunc)
        ch_versions.mix(GUNC_RUN.out.versions)

        // Make sure to keep directory in sync with modules.conf
        GUNC_RUN.out.maxcss_level_tsv
            .map { _meta, gunc_summary -> gunc_summary }
            .collectFile(
                name: "gunc_summary.tsv",
                keepHeader: true,
                storeDir: "${params.outdir}/GenomeBinning/QC/"
            )

        if (params.binqc_tool == 'checkm') {
            ch_input_to_mergecheckm = GUNC_RUN.out.maxcss_level_tsv.combine(CHECKM_QA.out.output, by: 0)

            GUNC_MERGECHECKM(ch_input_to_mergecheckm)
            ch_versions.mix(GUNC_MERGECHECKM.out.versions)

            // Make sure to keep directory in sync with modules.conf
            GUNC_MERGECHECKM.out.tsv
                .map { _meta, gunc_checkm_summary -> gunc_checkm_summary }
                .collectFile(
                    name: "gunc_checkm_summary.tsv",
                    keepHeader: true,
                    storeDir: "${params.outdir}/GenomeBinning/QC/"
                )
        }
    }*/

    emit:
    qc_summary    = qc_summary
    //multiqc_files = ch_multiqc_files
    versions      = ch_versions
}