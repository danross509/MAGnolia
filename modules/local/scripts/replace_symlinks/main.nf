#!/usr/bin/env nextflow

process REPLACE_SYMLINKS {
    label 'process_single'

    container ""
    conda ""

    //publishDir "${params.resultsDir}", mode: 'move'

    input:
        val placeholder     // This is a placeholder to ensure symlinks are not moved early : i.e. the counted output of the final process that requires them in the original work directory
        path symlinks       // The folder containing symlinks to replace

    output:        
        path "symlinks_replaced.txt"            , emit: symlinks_replaced, optional: true
        path "error_no_directory.txt"           , emit: error_no_directory, optional: true
        path "error_directory_empty.txt"        , emit: error_directory_empty, optional: true
        path "error_no_symlinks.txt"            , emit: error_no_symlinks, optional: true

    script:

/*    """
    target=\$(readlink -e $symlinks)
    sed -i '' $symlinks
    rm \$target
    """
*/
    """
    if [[ ! -d ${symlinks} ]]; then
        touch error_no_directory.txt
    elif [[ -z \$(ls -A ${symlinks}) ]]; then
        touch error_directory_empty.txt
    elif [[ -z \$(find ${symlinks} -type l) ]]; then
        touch error_no_symlinks.txt
    else
        working_directory=\$PWD
        cd ${symlinks}
        for f in \$(find . -maxdepth 1 -type l); do
            target=\$(readlink -e \$f)
            sed -i '' \$f
            rm \$target
        done

        cd \$working_directory
        touch symlinks_replaced.txt
    fi
    
    """

}