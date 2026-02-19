# Changelog

All notable changes to `nf-sarek-clinical` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.0.0] – 2024-XX-XX

### Initial release

#### Added
- Nextflow DSL2 pipeline wrapping post-sarek clinical interpretation workflow
- **Germline SNV** module: dbNSFP annotation → 90-gene cancer predisposition panel filter → TSV
- **Pharmacogenomics** module: DPYD/*2A/*13/HapB3 + UGT1A1/*28 rsID filter → TSV
- **Somatic SNV** module: dbNSFP annotation → somatic panel filter → CIViC annotation → TSV
- **Somatic CNV (ASCAT/CNVkit)** module: HOMDEL/HETDEL/AMP (38 somatic genes) + LOH 2nd-hit (90 germline genes)
- **TMB** module: PASS + HIGH/MODERATE + depth/AF/AD filter → mut/Mb
- **Mutational Signatures** module: SigProfilerMatrixGenerator + SigProfilerAssignment (COSMIC v3.4, SBS/DBS/ID)
- **Gene Coverage** module: mosdepth regions BED → per-gene mean/min/max coverage via R/data.table
- **Final Report** module: interactive HTML report with KPI dashboard (TMB, MSI, LOH, GI, WGD, ploidy, purity, coverage) + data tables
- Built-in gene panels:
  - `germline_90_cancer_predisposition.txt` (90 genes, ACMG SF v3.2 + ESMO)
  - `somatic_38_oncogenes.txt` (38 actionable somatic targets)
  - `somatic_snv_panel.txt` (60+ somatic SNV panel)
  - `pharmacogenomics_dpyd_variants.txt` (DPYD + UGT1A1 PGx variants)
- Docker and Singularity profile support
- Module-level skip flags for flexible partial runs
- Nextflow execution reports (timeline, trace, DAG)

#### Notes
- Requires nf-core/sarek ≥ 3.3.0 run with `haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff`
- CIViC annotation requires internet access during pipeline execution
- SigProfiler module requires a conda environment or dedicated container with SigProfilerAssignment
