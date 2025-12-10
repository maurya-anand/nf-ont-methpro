process ONT_BASECALL {
    publishDir "${params.outdir}/${meta.sampleid}/basecall", mode: 'copy'

    input:
    tuple val(meta), path(pod5_dir)

    output:
    tuple val(meta.sampleid), path("${meta.sampleid}.${meta.run_id}.trim.mod.bam"), emit: bam

    script:
    """
    if nvidia-smi &> /dev/null && [ \$(nvidia-smi -L | wc -l) -gt 0 ]; then
        DEVICE="cuda:all"
    else
        DEVICE="cpu"
    fi
    if ls ${pod5_dir}/*.fast5 1> /dev/null 2>&1; then
        echo "Found FAST5 files, converting to POD5..."
        mkdir -p pod5
        pod5 convert fast5 ${pod5_dir}/*.fast5 --output pod5/ --one-to-one ${pod5_dir}/
        INPUT_DIR="pod5"
    else
        INPUT_DIR="${pod5_dir}"
    fi
    total_threads=\$(nproc)
    available=\$((total_threads - 2))
    dorado basecaller \\
    ${params.basecall_model} \\
    \$INPUT_DIR \\
    --emit-moves \\
    --device "\$DEVICE" \\
    --modified-bases ${params.basecall_modifications} > \\
    ${meta.sampleid}.${meta.run_id}.raw.mod.bam

    dorado trim \\
    --threads \$available \\
    ${meta.sampleid}.${meta.run_id}.raw.mod.bam | \\
    samtools view \\
    -h -b -q 7 -e 'length(seq) >= 50' \\
    -@ \$available \\
    -o ${meta.sampleid}.${meta.run_id}.trim.mod.bam
    """
}
