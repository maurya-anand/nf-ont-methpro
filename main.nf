#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { ONT_BASECALL as BASECALL_READS } from './modules/local/basecall'
include { SUMMARY as BASECALL_SUMMARY } from './modules/local/basecall_summary'
include { FILTER_READS as FILTER_READS } from './modules/local/filter_reads'

workflow {
    ont_reads_ch = channel.fromPath(params.sample_sheet)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [
                sampleid: row.sampleid,
                input_dir: row.input_dir,
            ]
            def pod5_dir = file("${row.input_dir}/pod5")
            [meta, pod5_dir]
        }

    // Base modification calling of ONT reads. Raw fast5 files were converted to pod5 file format using ONT
    // pod5 convert fast5 v0.1.5. Then, pod5 files were used for simplex base and modification calling of 5mC and
    // 5hmC using the super-accurate basecalling model “dna_r9.4.1_e8_sup@v3.3” with ONT dorado basecaller
    // v0.7.2 (https://github.com/nanoporetech/dorado), and parameters “--emit-moves --device ‘cuda:all’ --no-trim
    // --modified-bases '5mCG_5hmCG'”. The resulting raw modbam files were input to dorado trim with default
    // parameters for adapter trimming. Then, reads smaller than 50 bases or with mean quality of less than 7 were
    // discarded with pysam v0.19.0 (https://github.com/pysam-developers/pysam).

    // Reads in modbam were transiently converted to fastq format and preserving all read tags with samtools fastq
    // v1.17, using the parameter “-T '*'”. The fastq reads were piped to Minimap2 v2.24-r1122 (Li, 2018) for
    // alignment to GRCh38 using the parameters “-a -x map-ont --rmq=yes --MD --cs -L -y”. After alignment,
    // unmapped reads were discarded using samtools view v1.17. Read group tags with the sequencing flow cell id
    // was added to modbam files using samtools addreplacerg. Modbams were, then, pooled by sample name
    // using samtools merge and sorted by coordinate and indexed with samtools sort and samtools index for
    // downstream analyses.
    // Reads overlapping the 294 gene regions and with a matching read id from LociTyper were retained and
    // assigned as haplotype 1 or haplotype2 in independent modbam files using samtools view. Modification ratios
    // were computed using ONT modkit pileup v0.2.2 (https://nanoporetech.github.io/modkit/) over each
    // haplotype with parameters “--sampling-frac 1 --seed 1234 “, and only restricting the pileup of modifications
    // to CpG dinucleotides. To annotate the modification ratios of each CpG, bedmethyl files were overlapped
    // with the 294 gene regions, using bedtools intersect v2.30.0
    // 38
    // .
    BASECALL_READS(ont_reads_ch)
    BASECALL_SUMMARY(TRIM_READS.out)
    FILTER_READS(TRIM_READS.out)
}
