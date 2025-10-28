process TRIM_READS {
        input:
        tuple(val(meta), path(reads))

        output:
        tuple(val(meta), path("${meta.sampleid}.${meta.flowcellid}.trim.mod.bam"))

        script:
        """
        total_threads=\$(nproc)
        threads=\$((total_threads - 2))
        dorado trim \\
        --threads \$threads
        ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam > \\
        ${meta.sampleid}.${meta.flowcellid}.trim.mod.bam
        """
}