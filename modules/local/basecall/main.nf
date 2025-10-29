process ONT_BASECALL {
    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta), path("${meta.sampleid}.${meta.flowcellid}.raw.mod.bam")

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
    --no-trim \\
    --modified-bases ${params.basecall_modifications} > \\
    ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam
    """
}
