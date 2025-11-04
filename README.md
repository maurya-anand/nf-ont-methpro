# nf-ont-methpro: Oxford Nanopore Methylation Profiling Pipeline

This pipeline processes Oxford Nanopore long-read sequencing data for methylation and variant calling using Nextflow DSL2. It is designed for sample-wise organization and modular analysis.

## Features

- Basecalling using `Dorado`
- Read alignment using `minimap2`
- Variant calling and haplotagging using `PEPPER-Margin-DeepVariant`
- Methylation calling using `modkit`
- Summary reporting using `MultiQC`

## Requirements

- Nextflow version 22.10.0 or higher
- Singularity or Docker (or a compatible container runtime)

## Quick Start

1. **Prepare your sample sheet** (CSV):

    ```csv
    sampleid,pod5_dir
    sample1,/path/to/sample1
    sample2,/path/to/sample2
    ```

2. **Configure your reference and parameters** in `nextflow.config`:

    ```groovy
    params {
        sample_sheet = "samplesheet.csv"
        basecall_model = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
        basecall_modifications = "5mCG_5hmCG"
        filter_reads = "-ml 50 -mq 7 -rd"
        sequencing_kit = "SQK-LSK110"
        reference = "GrCh38.fa"
        outdir = "results"
    }
    ```

3. **Run the pipeline:**

    ```bash
    nextflow run main.nf -profile docker
    ```

## Workflow Overview

The pipeline consists of the following main steps:

1. **Basecalling** (`ONT_BASECALL`): - Uses Dorado for basecalling and modified base detection. - Output: unaligned BAM with modifications.

2. **Alignment** (`ALIGNMENT`): - Aligns reads to the reference genome with minimap2. - Output: aligned BAM, BAM index, alignment stats, minimap2 log

3. **Variant Calling** (`VARIANT_CALL`): - Calls variants and haplotags reads using PEPPER-Margin-DeepVariant. - Output: haplotagged BAM, VCF, DeepVariant logs, visual report

4. **Methylation Calling** (`METHYLATION_CALL`): - Calls methylation using modkit, outputs both BED and bedGraph formats. - Output: methylation BED, bedGraph, modkit logs

5. **Summary Reporting** (`SUMMARY`): - Aggregates stats and logs using MultiQC. - Output: multiqc report

## Output Directory Structure

All outputs are organized **per sample**:

```text
results/
    sample1/
        basecall/
            sample1.raw.mod.bam
        alignment/
            sample1.aligned.sorted.bam
            sample1.aligned.sorted.bam.bai
            sample1.alignment.stats.txt
            sample1.minimap2.log
        variant_call/
            sample1.aligned.sorted.haplotagged.bam
            sample1.aligned.sorted.haplotagged.bam.bai
            sample1.vcf
            sample1.vcf.gz.tbi
            sample1.visual_report.htm
            sample1.pepper.margin.deepvariant.log
            logs/*.log
        methylation/
            sample1.methylation.calls.bed
            sample1.modkit.pileup.log
            sample1.modkit.pileup.bedgraph.log
            methylation_calls.bedgraph/*.bedgraph
    sample2/
        ...
    report/
        multiqc_report.html
        multiqc_data/
```

## Customization

- Adjust resource requirements and containers in `nextflow.config` and `conf/base.config`.
- Add or remove modules as needed for your workflow.

## Citation

todo
