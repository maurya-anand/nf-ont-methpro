# Changelog

All notable changes to this project will be documented in this file.

## [1.1.0] - 2025-11-05

### Added

- Haplotype-specific methylation analysis workflow
  - New `SPLIT_BAM` module to split haplotagged BAM files by haplotype (HP1, HP2, untagged)
  - Optional region-based filtering using BED file (`--regions_bed` parameter)
  - Haplotype-specific methylation calling with separate outputs for each haplotype
  - Support for targeted methylation analysis on specific genomic regions

### Changed

- Renamed `VARIANT_CALL` process to `VARIANT_CALL_AND_PHASING` to better reflect phasing functionality
- Updated `METHYLATION_CALL` process to accept haplotype-specific BAM inputs
- Modified workflow to perform methylation calling per haplotype instead of genome-wide
- Updated `SUMMARY` (REPORT) process to handle haplotype-specific methylation outputs (BED and bedGraph files)

### Enhanced

- Methylation calling now produces separate outputs for:
  - HP1 (haplotype 1) methylation calls
  - HP2 (haplotype 2) methylation calls
  - Both BED and bedGraph formats per haplotype
- Improved output organization with haplotype-specific subdirectories

## [1.0.0] - 2025-11-04

### Release Notes

- Initial release of the nf-ont-methpro pipeline for Oxford Nanopore methylation profiling.
- Implements modular Nextflow DSL2 workflow with sample-wise output organization.
- Major features:
  - Basecalling using Dorado
  - Read alignment with minimap2
  - Variant calling and haplotagging using PEPPER-Margin-DeepVariant
  - Methylation calling using modkit (BED and bedGraph outputs)
  - Summary reporting with MultiQC
- Supports execution with Docker, Singularity, Apptainer, and SLURM HPC environments.
