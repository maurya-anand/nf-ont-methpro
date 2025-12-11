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
    total_threads=\$(nproc)
    available=\$((total_threads - 2))
    if ls ${pod5_dir}/*.fast5 1> /dev/null 2>&1; then
        echo "Found FAST5 files, converting to POD5..."
        mkdir -p pod5
        export POD5_DEBUG=1
        timeout 4400 pod5 convert fast5 ${pod5_dir}/*.fast5 --output pod5/ --one-to-one ${pod5_dir}/ --threads \$available --force-overwrite 2>&1 || {
            echo "Batch conversion timed out or failed, converting files individually..."
            
            # Convert individually with per-file timeout
            for f in ${pod5_dir}/*.fast5; do
                fname=\$(basename "\$f" .fast5)
                timeout 360 pod5 convert fast5 "\$f" --output pod5/\${fname}.pod5 --force-overwrite 2>&1 | grep -v "Converting" || {
                    echo "WARNING: Skipping problematic file: \$f" >&2
                }
            done
        }
        
        # Check if we got any pod5 files
        if [ ! -d pod5 ] || [ -z "\$(ls -A pod5/*.pod5 2>/dev/null)" ]; then
            echo "ERROR: No POD5 files were created" >&2
            exit 1
        fi
        INPUT_DIR="pod5"
    else
        INPUT_DIR="${pod5_dir}"
    fi
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
