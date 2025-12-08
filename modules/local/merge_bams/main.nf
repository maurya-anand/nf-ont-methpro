process MERGE_BAMS {
    publishDir "${params.outdir}/${sampleid}/basecall", mode: 'copy'

    input:
    tuple val(sampleid), path(bams)

    output:
    tuple val(sampleid), path("${sampleid}.raw.mod.bam"), emit: bam

    script:
    def bam_list = bams instanceof List ? bams.join(' ') : bams
    def bam_count = bams instanceof List ? bams.size() : 1
    """
    if [ ${bam_count} -eq 1 ]; then
        mv ${bam_list} ${sampleid}.raw.mod.bam
    else
        samtools merge -@ \$(nproc) ${sampleid}.raw.mod.bam ${bam_list}
    fi
    """
}