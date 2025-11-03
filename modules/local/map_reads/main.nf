process ALIGNMENT {
    input:
    tuple val(meta), path(reads)

    output:
        tuple val(meta), 
            path("${meta.sampleid}.aligned.sorted.bam"), 
            path("${meta.sampleid}.aligned.sorted.bam.bai"), 
            path("${meta.sampleid}.alignment.stats.txt")

    script:
    """
        total_threads=\$(nproc)
        threads=\$((total_threads - 2))
        samtools fastq -T MM,ML,MN ${reads} | \
        minimap2 -ax map-ont -y --secondary=no -t \$threads ${params.reference} - 2> minimap2.log | \
        samtools view -b | \
        samtools sort -@ \$threads -o ${meta.sampleid}.aligned.sorted.bam && \
        samtools index ${meta.sampleid}.aligned.sorted.bam && \
        samtools flagstat ${meta.sampleid}.aligned.sorted.bam > ${meta.sampleid}.alignment.stats.txt
        """
}
