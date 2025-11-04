#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONT_BASECALL as BASECALLING } from './modules/local/base_call'
include { ALIGNMENT as MAPPING } from './modules/local/map_reads'
include { VARIANT_CALL as VARIANT_CALLING } from './modules/local/variant_call'
include { METHYLATION_CALL as METHYLATION_CALLING } from './modules/local/methylation_call'
include { SUMMARY as REPORT } from './modules/local/generate_report'

workflow {
    ont_reads_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sampleid,
                pod5_dir: row.pod5_dir,
            ]
            def pod5_dir = file("${row.pod5_dir}/pod5")
            [meta, pod5_dir]
        }
    BASECALLING(ont_reads_ch)
    MAPPING(BASECALLING.out.bam)
    VARIANT_CALLING(MAPPING.out)
    METHYLATION_CALLING(VARIANT_CALLING.out.bam)
    summary_stats = MAPPING.out.stats
    variant_logs = VARIANT_CALLING.out.logs
    methylation_bed = METHYLATION_CALLING.out.modbed.map { i -> i[0] }
    REPORT(
        summary_stats.collect(),
        variant_logs.collect(),
        methylation_bed.collect(),
    )
}
