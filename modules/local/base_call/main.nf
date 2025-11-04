process ONT_BASECALL {
    publishDir "${params.outdir}/${meta.sampleid}/basecall", mode: 'copy'

    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta.sampleid), path("${meta.sampleid}.raw.mod.bam"), emit: bam

    script:
    """
    if nvidia-smi &> /dev/null && [ \$(nvidia-smi -L | wc -l) -gt 0 ]; then
        DEVICE="cuda:all"
    else
        DEVICE="cpu"
    fi
    dorado basecaller \\
    ${params.basecall_model} \\
    ${pod5_dir}/ \\
    --emit-moves \\
    --device "\$DEVICE" \\
    --modified-bases ${params.basecall_modifications} > \\
    ${meta.sampleid}.raw.mod.bam
    """
}
