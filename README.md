# nf-ont-methpro: Oxford Nanopore Methylation Profiling Pipeline

A Nextflow DSL2 pipeline for processing Oxford Nanopore long-read sequencing data to generate haplotype-resolved variant calls and DNA methylation profiles.

## Features

- Basecalling from POD5 files using `Dorado`
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
   sampleid,data_dir
   sample1,/path/to/sample1/pod5/
   sample2,/path/to/sample2/pod5/
   ```

   **Note:** The pipeline expects the full path of the directory containing the `*.pod5` or `*.fast5` files. Your directory structure should look like:

   ```text
   /path/to/sample1/
   └── pod5/
       ├── file1.pod5
       ├── file2.pod5
       └── ...
   ```

   ```text
   /path/to/sample1/
   └── fast5/
       ├── file1.fast5
       ├── file2.fast5
       └── ...
   ```

2. **Target regions** (BED)

   ```tsv
   chr1    1000000 2000000
   chr2    3000000 4000000
   ```

3. **Configure your reference and parameters** in `nextflow.config`:

   ```groovy
   params {
       sample_sheet = "samplesheet.csv"
       basecall_model = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
       basecall_modifications = "5mCG_5hmCG"
       reference = "GrCh38.fa"
       regions_bed = "regions.bed"
       outdir = "results"
   }
   ```

4. **Run the pipeline:**

   ```bash
   nextflow run main.nf -profile docker
   ```

> [!NOTE]
>
> - The reference genome file must have a corresponding index (`.fai`) file in the same directory. Generate it with: `samtools faidx GrCh38.fa`
> - For information on selecting basecall models and modification options, refer to the [Dorado documentation](https://software-docs.nanoporetech.com/dorado/latest/models/selection/#selecting-modified-base-models).
> - When `--regions_bed` is provided:
>   - Only reads overlapping the specified regions are extracted during haplotype splitting.
>   - Methylation calling and DMR detection are restricted to these regions.

## Workflow Overview

The pipeline consists of the following main steps:

1. **Basecalling** (`ONT_BASECALL`):

   - Performs basecalling and modified base detection using `Dorado`. The output BAM file contains MM (modification type) and ML (modification likelihood) tags.
   - Automatically detects and converts FAST5 files to POD5 format if needed before basecalling.
   - Output: Unaligned BAM with methylation tags.

2. **Alignment** (`ALIGNMENT`):

   - Aligns reads to the reference genome with `minimap2`.
   - Output: Aligned BAM, alignment statistics.

3. **Variant Calling and Phasing** (`VARIANT_CALL`):

   - Calls variants and haplotags reads using `PEPPER-Margin-DeepVariant`.
   - Output: Haplotagged BAM and VCF.

4. **Haplotype Splitting** (`SPLIT_BAM`):

   - Splits haplotagged BAM into haplotype-specific files (HP1, HP2, and untagged reads).
   - Optionally filters reads by genomic regions when `--regions_bed` is provided.
   - Output: Three BAM files per sample (HP1, HP2, untagged) with their indices.

5. **Haplotype-Resolved Methylation Calling** (`METHYLATION_CALL`):

   - Extracts and aggregates methylation calls separately for each haplotype using `modkit`.
   - Generates both BED and bedGraph formats.
   - Output: Methylation BED and bedGraph files for HP1 and HP2.

6. **Differentially Methylated Regions (DMR)** (`DMR_CALL`):

   - Identifies regions with significant methylation differences between haplotypes using `modkit dmr pair`.
   - Output: DMR BED file and analysis log.

7. **Summary** (`SUMMARY`):
   - Generates aggregated QC report using `MultiQC`.
   - Output: MultiQC report

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
            sample1.bam.mod.tags.txt
            sample1.minimap2.log
        variant_call/
            sample1.aligned.sorted.haplotagged.bam
            sample1.aligned.sorted.haplotagged.bam.bai
            sample1.vcf
            sample1.vcf.gz.tbi
            sample1.visual_report.html
            sample1.pepper.margin.deepvariant.log
            logs/*.log
        haplotypes/
            sample1.HP1.bam
            sample1.HP1.bam.bai
            sample1.HP2.bam
            sample1.HP2.bam.bai
            sample1.untagged.bam
            sample1.untagged.bam.bai
        methylation_haplotypes/
            sample1.HP1.methylation.calls.bed
            sample1.HP1.modkit.pileup.bed.log
            sample1.HP1.modkit.pileup.bedgraph.log
            methylation_calls.HP1.bedGraph/*.bedgraph
            sample1.HP2.methylation.calls.bed
            sample1.HP2.modkit.pileup.bed.log
            sample1.HP2.modkit.pileup.bedgraph.log
            methylation_calls.HP2.bedGraph/*.bedgraph
        differentially_methylated_regions/
            sample1.haplotype.differentially.methylated.regions.bed
            sample1.dmr.log
    sample2/
        ...
    report/
        multiqc_report.html
        multiqc_data/
```

> [!NOTE]
>
> - For `*.methylation.calls.bed` file format and column descriptions, see [modkit documentation](https://nanoporetech.github.io/modkit/intro_pileup.html#bedmethyl-column-descriptions).
> - For `*.differentially.methylated.regions.bed` file format and column descriptions, see [modkit documentation](https://nanoporetech.github.io/modkit/intro_dmr.html#differential-methylation-output-format).

## Customization

- Adjust resource requirements and containers in `nextflow.config` and `conf/base.config`.
- Add or remove modules as needed for your workflow.

## Components

| Component                                                            | Version        |
| -------------------------------------------------------------------- | -------------- |
| [Dorado](https://github.com/nanoporetech/dorado)                     | 1.2.0+f9443bb8 |
| [minimap2](https://github.com/lh3/minimap2)                          | 2.30-r1287     |
| [samtools](http://www.htslib.org/)                                   | 1.13           |
| [PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper) | r0.8           |
| [modkit](https://github.com/nanoporetech/modkit)                     | 0.5.0          |
| [MultiQC](https://multiqc.info/)                                     | 1.32           |
| [bedtools](https://bedtools.readthedocs.io/)                         | 2.30.0         |
| [bcftools](http://samtools.github.io/bcftools/)                      | 1.13           |

### Container Images

- PEPPER-Margin-DeepVariant: `kishwars/pepper_deepvariant:r0.8`

## Citation

todo
