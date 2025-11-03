process FILTER_READS {
        input:
        tuple val(meta), path(bam)

        output:
        tuple val(meta), path("${meta.sampleid}.${meta.flowcellid}.filtered.mod.pass.bam")

        script:
        def filtered_reads = "${meta.sampleid}.${meta.flowcellid}.filtered.mod.bam"

        """
        micromamba run -n base ont_tools.py filter_bam \\
        ${params.filter_reads}
        -i ${bam} \\
        -o ${filtered_reads}
        """
}
