process DMR_CALL {
    publishDir "${params.outdir}/${sampleid}/differentially_methylated_regions", mode: 'copy'

    input:
    tuple val(sampleid), path(hp1_bed), path(hp2_bed)
    path reference
    path reference_fai

    output:
    tuple val(sampleid), path("${sampleid}.haplotype.differentially.methylated.regions.bed"), emit: dmr
    path("${sampleid}.dmr.log"), emit: log

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    
    # Compress bedmethyl files if not already compressed
    if [[ ! -f ${hp1_bed}.gz ]]; then
        bgzip -c ${hp1_bed} > ${hp1_bed}.gz
        tabix -p bed ${hp1_bed}.gz
    fi
    
    if [[ ! -f ${hp2_bed}.gz ]]; then
        bgzip -c ${hp2_bed} > ${hp2_bed}.gz
        tabix -p bed ${hp2_bed}.gz
    fi
    
    modkit dmr pair \
      -a ${hp1_bed}.gz \
      -b ${hp2_bed}.gz \
      -o ${sampleid}.haplotype.differentially.methylated.regions.bed \
      --ref ${reference} \
      --base C \
      --log-filepath ${sampleid}.dmr.log
    """
}