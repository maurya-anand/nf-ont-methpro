process METHYLATION_CALL {
    publishDir "${params.outdir}/${sampleid}/methylation", mode: 'copy'

    input:
    tuple val(sampleid), path(bam), path(bai)

    output:
    tuple path("${sampleid}.methylation.calls.bed"), path("${sampleid}.modkit.pileup.log"), emit: modbed
    tuple path("${sampleid}.modkit.pileup.bedgraph.log"), path("methylation_calls.bedgraph/*.bedgraph"), emit: bedgraph

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    mkdir -p methylation_calls.bedgraph
    
    modkit pileup \
    ${bam} \
    methylation_calls.bedgraph \
    --bedgraph \
    --ref ${params.reference} \
    --cpg \
    --combine-strands \
    --ignore h \
    --threads \${threads} &> ${sampleid}.modkit.pileup.bedgraph.log
    
    modkit pileup \
    ${bam} \
    ${sampleid}.methylation.calls.bed \
    --ref ${params.reference} \
    --cpg \
    --combine-strands \
    --ignore h \
    --threads \${threads} &> ${sampleid}.modkit.pileup.log
    """
}
