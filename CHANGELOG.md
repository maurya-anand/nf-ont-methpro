# Changelog

All notable changes to this project will be documented in this file.

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
