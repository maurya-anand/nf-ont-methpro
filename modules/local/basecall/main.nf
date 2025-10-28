process ONT_BASECALL {
    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta), path("${meta.sampleid}.${meta.flowcellid}.raw.mod.bam")

    script:
    """
    dorado basecaller \\
    ${params.basecall_model} \\
    ${pod5_dir}/ \\
    --emit-moves \\
    --device "cuda:all" \\
    --no-trim \\
    --modified-bases ${params.basecall_modifications} > \\
    ${meta.sampleid}.${meta.flowcellid}.raw.mod.bam
    """
}
