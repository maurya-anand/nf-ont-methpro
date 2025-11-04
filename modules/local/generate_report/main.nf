process SUMMARY {
    publishDir "${params.outdir}/report", mode: 'copy'

    input:
        path(stats_files)
        path(log_files)
        path(methylation_bed), optional: true

    output:
        path("multiqc_report.html"), emit: multiqc
        path("multiqc_data"), emit: multiqc_data

    script:
    """
    mkdir -p multiqc_input
    cp \${stats_files} multiqc_input/
    cp \${log_files} multiqc_input/
    if [ -n "\${methylation_bed}" ]; then
        cp \${methylation_bed} multiqc_input/
    fi
    multiqc multiqc_input -o .
    """
}