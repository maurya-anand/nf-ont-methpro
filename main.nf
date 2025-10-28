#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ONT_BASECALL as BASECALL_READS } from './modules/local/basecall'

workflow {
        ont_reads_ch = channel.fromPath( params.sample_sheet  )
           .splitCsv( header:true,sep:"," )
           .map { 
            row -> 
            def meta = [
                sampleid:   row.sampleid,
                flowcellid: row.flowcellid,
                input_dir:  row.input_dir ]
            def pod5_dir = file("${row.input_dir}/pod5")
            [ meta, pod5_dir ] 
           }
    BASECALL_READS( ont_reads_ch )
}