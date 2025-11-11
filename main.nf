#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONT_BASECALL as METHYLATED_BASECALLING } from './modules/local/base_call'
include { ALIGNMENT as MAPPING } from './modules/local/map_reads'
include { VARIANT_CALL as VARIANT_CALL_AND_PHASING } from './modules/local/variant_call'
include { SPLIT_BAM as EXTRACT_READS_BY_HAPLOTYPE } from './modules/split_bam/'
include { METHYLATION_CALL as METHYLATION_CALLING } from './modules/local/methylation_call'
include { DMR_CALL as DIFFERENTIAL_MODIFICATION } from './modules/local/dmr_call'
include { SUMMARY as REPORT } from './modules/local/generate_report'

workflow {
    ont_reads_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sampleid,
                data_dir: row.data_dir,
            ]
            def pod5_dir = file("${row.data_dir}/pod5")
            [meta, pod5_dir]
        }
    METHYLATED_BASECALLING(ont_reads_ch)
    MAPPING(METHYLATED_BASECALLING.out.bam, file(params.reference), file("${params.reference}.fai"))
    VARIANT_CALL_AND_PHASING(MAPPING.out.bam, file(params.reference), file("${params.reference}.fai"))
    variant_logs_ch = VARIANT_CALL_AND_PHASING.out.main_log.mix(VARIANT_CALL_AND_PHASING.out.logs.flatten())
    EXTRACT_READS_BY_HAPLOTYPE(VARIANT_CALL_AND_PHASING.out.bam, file(params.regions_bed))
    hp1_ch = EXTRACT_READS_BY_HAPLOTYPE.out.haplotype1.map { sampleid, bam, bai -> tuple(sampleid, bam, bai, "HP1") }
    hp2_ch = EXTRACT_READS_BY_HAPLOTYPE.out.haplotype2.map { sampleid, bam, bai -> tuple(sampleid, bam, bai, "HP2") }
    haplotype_bams_ch = hp1_ch.mix(hp2_ch)
    METHYLATION_CALLING(haplotype_bams_ch, file(params.reference), file("${params.reference}.fai"))
    hp1_modbed = METHYLATION_CALLING.out.modbed
        .filter { _sampleid, haplotype, _bed, _log -> haplotype == "HP1" }
        .map { sampleid, _haplotype, bed, _log -> tuple(sampleid, bed) }
    hp2_modbed = METHYLATION_CALLING.out.modbed
        .filter { _sampleid, haplotype, _bed, _log -> haplotype == "HP2" }
        .map { sampleid, _haplotype, bed, _log -> tuple(sampleid, bed) }
    dmr_input = hp1_modbed.join(hp2_modbed)
    DIFFERENTIAL_MODIFICATION(dmr_input, file(params.reference), file("${params.reference}.fai"))
    REPORT(
        MAPPING.out.stats.collect(),
        variant_logs_ch.collect(),
        METHYLATION_CALLING.out.modbed.map { _sampleid, _haplotype, bed, _log -> bed },
        METHYLATION_CALLING.out.bedgraph.map { sampleid, haplotype, bedgraphs, _log ->
            [bedgraphs]
                .flatten()
                .collect { bedgraph_file ->
                    file(bedgraph_file).copyTo("${sampleid}.${haplotype}.${file(bedgraph_file).name}")
                }
        }.flatten().collect(),
        DIFFERENTIAL_MODIFICATION.out.log,
    )
}
