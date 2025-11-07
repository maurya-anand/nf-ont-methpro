process SPLIT_BAM {
    publishDir "${params.outdir}/${sampleid}/haplotypes", mode: 'copy'

    input:
    tuple val(sampleid), path(bam), path(bai)

    output:
    tuple val(sampleid), path("${sampleid}.HP1.bam"), path("${sampleid}.HP1.bam.bai"), emit: haplotype1
    tuple val(sampleid), path("${sampleid}.HP2.bam"), path("${sampleid}.HP2.bam.bai"), emit: haplotype2
    tuple val(sampleid), path("${sampleid}.untagged.bam"), path("${sampleid}.untagged.bam.bai"), emit: untagged

    script:
    def bed_option = params.regions_bed ? "-L ${params.regions_bed}" : ""
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    # Split haplotype 1 (HP:i:1)
    samtools view -@ \$threads -h -d HP:1 ${bed_option} ${bam} -o ${sampleid}.HP1.bam
    samtools index -@ \$threads ${sampleid}.HP1.bam

    # Split haplotype 2 (HP:i:2)
    samtools view -@ \$threads -h -d HP:2 ${bed_option} ${bam} -o ${sampleid}.HP2.bam
    samtools index -@ \$threads ${sampleid}.HP2.bam

    # Extract untagged reads (no HP tag)
    samtools view -@ \$threads -h -e '!([HP])' ${bed_option} ${bam} -o ${sampleid}.untagged.bam
    samtools index -@ \$threads ${sampleid}.untagged.bam
    """
}