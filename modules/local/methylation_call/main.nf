process METHYLATION_CALL {
    publishDir "${params.outdir}/${sampleid}/methylation_haplotypes", mode: 'copy'

    input:
    tuple val(sampleid), path(bam), path(bai), val(haplotype)
    path reference
    path reference_fai

    output:
    tuple val(sampleid), val(haplotype), path("${sampleid}.${haplotype}.methylation.calls.bed"), path("${sampleid}.${haplotype}.modkit.pileup.bedmethyl.log"), emit: modbed
    tuple val(sampleid), val(haplotype), path("methylation_calls.${haplotype}.bedgraph/*.bedgraph"), path("${sampleid}.${haplotype}.modkit.pileup.bedgraph.log"), emit: bedgraph

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    mkdir -p methylation_calls.${haplotype}.bedgraph
    
    modkit pileup \
    ${bam} \
    methylation_calls.${haplotype}.bedgraph \
    --bedgraph \
    --ref ${reference} \
    --cpg \
    --combine-strands \
    --ignore h \
    --threads \${threads} &> ${sampleid}.${haplotype}.modkit.pileup.bedgraph.log
    
    modkit pileup \
    ${bam} \
    ${sampleid}.${haplotype}.methylation.calls.bed \
    --ref ${reference} \
    --cpg \
    --combine-strands \
    --ignore h \
    --threads \${threads} &> ${sampleid}.${haplotype}.modkit.pileup.bed.log
    """
}