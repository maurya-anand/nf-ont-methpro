process SUMMARY {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
        path(stats_files)
        path(log_files)
        path(methylation_bed)
        path(methylation_bedgraph)
        path(dmr)

    output:
        path("multiqc_report.html"), emit: multiqc
        path("multiqc_data"), emit: multiqc_data

    script:
    """
    mkdir -p multiqc_input
    cp ${dmr} multiqc_input/

    for file in ${stats_files}; do
        cp \$file multiqc_input/
    done
    
    for file in ${log_files}; do
        cp \$file multiqc_input/
    done
    
    if [ -n "${methylation_bed}" ] && [ "${methylation_bed}" != "null" ]; then
        for file in ${methylation_bed}; do
            [ -f \$file ] && cp \$file multiqc_input/
        done
    fi

    if [ -n "${methylation_bedgraph}" ] && [ "${methylation_bedgraph}" != "null" ]; then
        for file in ${methylation_bedgraph}; do
            [ -f \$file ] && cp \$file multiqc_input/
        done
    fi
    
    multiqc multiqc_input -o .
    """
}