process ALIGNMENT {
    publishDir "${params.outdir}/${sampleid}/alignment", mode: 'copy'

    input:
    tuple val(sampleid), path(reads)
    path reference
    path reference_fai

    output:
    tuple val(sampleid), path("${sampleid}.aligned.sorted.bam"), path("${sampleid}.aligned.sorted.bam.bai"), emit: bam
    path("${sampleid}.alignment.stats.txt"), emit: stats
    path("${sampleid}.bam.mod.tags.txt"), emit: mod_tags
    path("${sampleid}.minimap2.log"), emit: logs

    script:
    """
    total_threads=\$(nproc)
    available=\$((total_threads - 2))
    minimap_threads=\$(( available * 7 / 10 ))
    sort_threads=\$(( available * 2 / 10 ))
    fastq_threads=\$(( available * 1 / 10 ))
    [ \$minimap_threads -lt 1 ] && minimap_threads=1
    [ \$sort_threads -lt 1 ] && sort_threads=1
    [ \$fastq_threads -lt 1 ] && fastq_threads=1
    samtools view -@ \$available ${reads} | head -1 | grep -o "MM:Z:\\|ML:B:" > ${sampleid}.bam.mod.tags.txt
    samtools fastq -@ \$fastq_threads -T MM,ML,MN ${reads} | \
    minimap2 -ax map-ont -y --secondary=no -t \$minimap_threads ${reference} - 2> ${sampleid}.minimap2.log | \
    samtools sort -@ \$sort_threads -o ${sampleid}.aligned.sorted.bam -
    samtools index -@ \$available ${sampleid}.aligned.sorted.bam
    samtools flagstat -@ \$available ${sampleid}.aligned.sorted.bam > ${sampleid}.alignment.stats.txt
    """
}
