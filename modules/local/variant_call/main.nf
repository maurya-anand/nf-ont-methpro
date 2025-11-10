process VARIANT_CALL {
    publishDir "${params.outdir}/${sampleid}/variant_call", mode: 'copy'

    input:
    tuple val(sampleid), path(bam), path(bai)
    path reference
    path reference_fai

    output:
    tuple val(sampleid), path("${sampleid}.aligned.sorted.haplotagged.bam"), path("${sampleid}.aligned.sorted.haplotagged.bam.bai"), emit: bam
    tuple val(sampleid), path("${sampleid}.vcf.gz"), path("${sampleid}.vcf.gz.tbi"), emit: vcf
    path("${sampleid}.visual_report.html"), emit: report
    tuple path("${sampleid}.pepper.margin.deepvariant.log"), path("logs/*.log"), emit: logs

    script:
    """
    total_threads=\$(nproc)
    threads=\$((total_threads - 2))
    export TMPDIR=\$PWD/tmp
    export TEMP=\$TMPDIR
    export TMP=\$TMPDIR
    export PARALLEL_TMPDIR=\$TMPDIR
    mkdir -p \$TMPDIR
    run_pepper_margin_deepvariant call_variant \
    -b ${sampleid}.aligned.sorted.bam \
    -f ${reference} \
    -o . \
    --gvcf \
    --sample_name ${sampleid} \
    -p ${sampleid} \
    -t \${threads} \
    --dv_sort_by_haplotypes true \
    --keep_intermediate_bam_files true \
    --ont_r10_q20 >> ${sampleid}.pepper.margin.deepvariant.log 2>&1
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam \
    ${sampleid}.aligned.sorted.haplotagged.bam
    mv intermediate_files/PHASED.PEPPER_MARGIN.haplotagged.bam.bai \
    ${sampleid}.aligned.sorted.haplotagged.bam.bai
    """
}
