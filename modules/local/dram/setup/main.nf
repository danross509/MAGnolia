#!/usr/bin/env nextflow

process DRAM_SETUP {
    //tag "${meta.id}"
    label 'process_high'

    container ""
    conda "${moduleDir}/environment.yml"

    //publishDir "${params.resultsDir}/KRAKEN2/${file_type}/${meta.id}", mode: 'symlink'

    input:
        val config_loc
        val db_dir
        val kegg_loc
        val skip_uniref

    output:
        path "dram_setup_finished.txt"

    script:
    def args = task.ext.args ?: ""
    def kegg = kegg_loc ? "--kegg_loc ${kegg_loc}" : ""
    def uniref = skip_uniref ? "--skip_uniref" : ""
    def command = config_loc ? "import_config --config_loc  ${config_loc}" : "prepare_databases --output_dir ${db_dir} ${kegg} ${uniref} ${args} --threads $task.cpus"

    """
    #Replace DRAM-setup.py with the fixed version in bin
    if [[ -f \$CONDA_PREFIX/bin/DRAM-setup.py ]]; then
        rm \$CONDA_PREFIX/bin/DRAM-setup.py
    fi 

    cp ${moduleDir}/DRAM-setup_fixed.py \$CONDA_PREFIX/bin/DRAM-setup.py 

    python \$CONDA_PREFIX/bin/DRAM-setup.py \\
        $command 

    touch dram_setup_finished.txt
    """
}
