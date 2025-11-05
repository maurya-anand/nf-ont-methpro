#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONT_BASECALL as BASECALLING } from './modules/local/base_call'
include { ALIGNMENT as MAPPING } from './modules/local/map_reads'
include { VARIANT_CALL as VARIANT_CALLING } from './modules/local/variant_call'
include { METHYLATION_CALL as METHYLATION_CALLING } from './modules/local/methylation_call'
include { SPLIT_BAM as EXTRACT_READS_BY_HAPLOTYPE } from './modules/split_bam/'
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
    VARIANT_CALLING(MAPPING.out.bam)
    EXTRACT_READS_BY_HAPLOTYPE(MAPPING.out.bam)
    hp1_ch = EXTRACT_READS_BY_HAPLOTYPE.out.haplotype1.map { sampleid, bam, bai -> tuple(sampleid, bam, bai, "HP1") }
    hp2_ch = EXTRACT_READS_BY_HAPLOTYPE.out.haplotype2.map { sampleid, bam, bai -> tuple(sampleid, bam, bai, "HP2") }
    haplotype_bams_ch = hp1_ch.mix(hp2_ch)
    METHYLATION_CALLING(haplotype_bams_ch)
    REPORT(
        MAPPING.out.stats.collect(),
        VARIANT_CALLING.out.logs.collect(),
        METHYLATION_CALLING.out.modbed.map { bed, _log -> bed },
        METHYLATION_CALLING.out.bedgraph.map { bedgraphs, _log -> bedgraphs }
    )
}
