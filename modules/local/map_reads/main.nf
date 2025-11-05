process ALIGNMENT {
    publishDir "${params.outdir}/${sampleid}/alignment", mode: 'copy'

    input:
    tuple val(sampleid), path(reads)

    output:
    tuple val(sampleid), path("${sampleid}.aligned.sorted.bam"), path("${sampleid}.aligned.sorted.bam.bai"), emit: bam
    path("${sampleid}.alignment.stats.txt"), emit: stats
    path("${sampleid}.bam.mod.tags.txt"), emit: mod_tags
    path("${sampleid}.minimap2.log"), emit: logs

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    samtools view -@ \$threads ${reads}.basecalls.mod.bam | head -1 | grep -o "MM:Z:\\|ML:B:" > ${sampleid}.bam.mod.tags.txt
    samtools fastq -T MM,ML,MN ${reads} | \
    minimap2 -ax map-ont -y --secondary=no -t \$threads ${params.reference} - 2> ${sampleid}.minimap2.log | \
    samtools view -@ \$threads -b | \
    samtools sort -@ \$threads -o ${sampleid}.aligned.sorted.bam && \
    samtools index ${sampleid}.aligned.sorted.bam && \
    samtools flagstat ${sampleid}.aligned.sorted.bam > ${sampleid}.alignment.stats.txt
    """
}
