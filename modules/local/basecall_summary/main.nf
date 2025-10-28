process SUMMARY {
        input:
        tuple val(meta), path(reads)

        output:
        tuple val(meta), path("summary.${meta.sampleid}.${meta.flowcellid}.txt.gz")

        script:
        """
        dorado summary \\
        ${reads} | gzip -  > summary.${meta.sampleid}.${meta.flowcellid}.txt.gz
        """
}
