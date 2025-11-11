process METHYLATION_CALL {
    publishDir "${params.outdir}/${sampleid}/methylation_haplotypes", mode: 'copy'

    input:
    tuple val(sampleid), path(bam), path(bai), val(haplotype)
    path reference
    path reference_fai

    output:
    tuple val(sampleid), val(haplotype), path("${sampleid}.${haplotype}.methylation.calls.bed"), path("${sampleid}.${haplotype}.modkit.pileup.bed.log"), emit: modbed
    tuple val(sampleid), val(haplotype), path("${sampleid}.${haplotype}.methylation.calls.bedgraph/*.bedgraph"), path("${sampleid}.${haplotype}.modkit.pileup.bedgraph.log"), emit: bedgraph

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    mkdir -p ${sampleid}.${haplotype}.methylation.calls.bedgraph
    
    modkit pileup \
    ${bam} \
    ${sampleid}.${haplotype}.methylation.calls.bedgraph \
    --bedgraph \
    --ref ${reference} \
    --cpg \
    --combine-strands \
    --ignore h \
    --threads \${threads} &> ${sampleid}.${haplotype}.modkit.pileup.bedgraph.log
    for file in ${sampleid}.${haplotype}.methylation.calls.bedgraph/*.bedgraph; do
        basename=\$(basename \$file)
        mv \$file ${sampleid}.${haplotype}.methylation.calls.bedgraph/${sampleid}.${haplotype}.\$basename
    done
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