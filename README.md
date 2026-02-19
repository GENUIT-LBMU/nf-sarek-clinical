# nf-sarek-clinical

**Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exome outputs.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Overview

`nf-sarek-clinical` is a Nextflow DSL2 pipeline that takes the output of [nf-core/sarek](https://nf-co.re/sarek) (tumour–normal WES) and performs:

| Module | Analysis | Input (Sarek) | Output |
|--------|----------|---------------|--------|
| **Germline SNV** | dbNSFP annotation → 90-gene panel + HIGH/MODERATE filter | HaplotypeCaller snpEff VCFs | `*.germline_snvs.tsv` |
| **Pharmacogenomics** | DPYD/UGT1A1 rsID filter | HaplotypeCaller snpEff VCFs | `*.pharmaco.tsv` |
| **Somatic SNV** | dbNSFP + somatic panel filter + CIViC annotation | Mutect2 snpEff VCFs | `*.somatic_snvs.tsv`, `*.civic.exploded.tsv` |
| **Somatic CNV (ASCAT/CNVkit)** | HOMDEL/HETDEL/AMP (38-gene panel) + LOH 2nd-hit (90 genes) | ASCAT segments+metrics, CNVkit .cns | `*_somatic_cnv_report.tsv/txt` |
| **TMB** | PASS + impact + depth filter → mut/Mb | Mutect2 snpEff VCFs | `*.tmb.tsv` |
| **Mutational Signatures** | SigProfilerMatrixGenerator + SigProfilerAssignment (COSMIC v3.4) | Mutect2 VCFs | SBS/DBS/ID profiles + assignment plots |
| **Gene Coverage** | Per-gene mean/min/max coverage from mosdepth | mosdepth regions BED | `*_gene_coverage.tsv` |
| **HTML Report** | Interactive clinical report with gauges and tables | All above outputs | `*_molecular_report.html` |

## Pipeline DAG

```
Sarek outputs
    │
    ├─ annotation/haplotypecaller/  ──► GERMLINE_POSTPROCESS ──► snv_tsv, pharmaco_tsv
    │                                     (dbNSFP → filter → extract)
    │
    ├─ annotation/mutect2/          ──► SOMATIC_POSTPROCESS
    │                                     ├── SNV (dbNSFP → CIViC → extract)  ──► somatic_snv_tsv, civic_tsv
    │                                     ├── CNV (ASCAT+CNVkit interpret)     ──► cnv_tsv, cnv_txt
    │                                     ├── TMB (PASS filter → count)        ──► tmb_tsv
    │                                     └── Signatures (SigProfiler)         ──► assignment plots
    │
    ├─ variant_calling/cnvkit/      ──┘
    ├─ variant_calling/ascat/       ──┘
    │
    ├─ reports/mosdepth/            ──► GENE_COVERAGE ──► gene_coverage_tsv
    │
    └─ (all outputs)                ──► FINAL_REPORT ──► molecular_report.html
```

## Requirements

- Nextflow ≥ 23.04.0
- Docker or Singularity
- nf-core/sarek ≥ 3.3.0 run completed with tools: `haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff`

### External databases (not included)

| File | Source | Notes |
|------|--------|-------|
| `dbNSFP4.x.txt.gz` | [Liu et al.](https://sites.google.com/site/jpopgen/dbNSFP) | ~30 GB; bgzip + tabix indexed |
| `SnpSift.jar` | [snpEff package](https://pcingola.github.io/SnpEff/) | Same version as used in Sarek |
| WES capture BED | Kit vendor (e.g. Agilent SureSelect) | For TMB panel Mb calculation |
| Gene coordinates BED | UCSC/Ensembl | For CNV gene intersection |

## Quick Start

```bash
# Clone
git clone https://github.com/your-org/nf-sarek-clinical.git
cd nf-sarek-clinical

# Run with Docker
nextflow run main.nf \
  -profile docker \
  --sarek_outdir /data/runs/MY_RUN/results \
  --run_name     MY_RUN \
  --dbnsfp       /data/dbNSFP4.1a.txt.gz \
  --snpsift_jar  /opt/snpEff/SnpSift.jar \
  --panel_bed    /data/hg38/S07604514_Regions.bed \
  --outdir       /results/MY_RUN_clinical
```

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--sarek_outdir` | Path to nf-core/sarek results directory |
| `--run_name` | Run identifier (used in output filenames) |
| `--dbnsfp` | Path to dbNSFP `.txt.gz` (tabix indexed) |
| `--snpsift_jar` | Path to `SnpSift.jar` |

### Optional (with defaults)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--panel_bed` | — | WES capture BED for TMB Mb calculation |
| `--gene_coords_bed` | — | BED with gene coordinates for CNV annotation |
| `--germline_panel` | `assets/gene_panels/germline_90_cancer_predisposition.txt` | 90-gene hereditary cancer panel |
| `--somatic_cnv_panel` | `assets/gene_panels/somatic_38_oncogenes.txt` | 38-gene somatic CNV panel |
| `--somatic_snv_panel` | `assets/gene_panels/somatic_snv_panel.txt` | Somatic SNV gene list |
| `--pharmaco_variants` | `assets/pharmacogenomics_dpyd_variants.txt` | PGx rsID list |
| `--genome_build` | `GRCh38` | `GRCh38` or `GRCh37` |
| `--amp_threshold` | `5` | Copy number threshold for amplification calls |
| `--tmb_min_dp` | `10` | Min tumor DP for TMB variants |
| `--tmb_min_af` | `0.05` | Min VAF for TMB variants |
| `--tmb_min_alt_reads` | `10` | Min ALT read count for TMB |
| `--germline_min_dp` | `5` | Min DP for germline SNV filter |
| `--outdir` | `./results` | Output directory |

### Module toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_germline_snv` | `false` | Skip germline SNV postprocessing |
| `--skip_somatic_snv` | `false` | Skip somatic SNV postprocessing |
| `--skip_somatic_cnv` | `false` | Skip somatic CNV interpretation |
| `--skip_tmb` | `false` | Skip TMB calculation |
| `--skip_signatures` | `false` | Skip mutational signature analysis |
| `--skip_pharmaco` | `false` | Skip pharmacogenomics |
| `--skip_coverage` | `false` | Skip gene coverage |
| `--skip_report` | `false` | Skip HTML report generation |

## Output Structure

```
results/
├── germline/
│   ├── snv/
│   │   ├── dbnsfp/    *.dbnsfp.vcf.gz
│   │   ├── filtered/  *.filtered.vcf.gz
│   │   └── extract/   *.germline_snvs.tsv
│   └── pharmaco/      *.pharmaco.tsv
├── somatic/
│   ├── snv/
│   │   ├── dbnsfp/         *.dbnsfp.vcf.gz
│   │   ├── filtered/       *.filtered.vcf.gz
│   │   ├── civic/          *.civic.vcf.gz
│   │   ├── extract/        *.somatic_snvs.tsv
│   │   └── civic_extract/  *.civic.exploded.tsv
│   ├── cnv/ascat_cnvkit/
│   │   ├── *_somatic_cnv_report.tsv
│   │   └── *_somatic_cnv_report.txt
│   ├── tmb/
│   │   └── *.tmb.tsv
│   └── signatures/
│       ├── filtered_vcfs/
│       ├── matrix_profiles/
│       └── assignment/
├── qc/
│   └── gene_coverage/  *_gene_coverage.tsv
├── reporting/
│   ├── *_molecular_report.html
│   └── *_molecular_report.txt
└── pipeline_info/
    ├── execution_timeline_*.html
    ├── execution_report_*.html
    └── execution_trace_*.txt
```

## Gene Panels

Built-in panels can be found in `assets/gene_panels/`:

- **`germline_90_cancer_predisposition.txt`** — 90 hereditary cancer predisposition genes (ACMG SF v3.2 compliant). Used for germline SNV filtering and ASCAT LOH 2nd-hit detection.
- **`somatic_38_oncogenes.txt`** — 38 actionable somatic genes. Used for CNVkit HOMDEL/HETDEL/AMP reporting.
- **`somatic_snv_panel.txt`** — Somatic SNV panel covering oncogenes, TSGs, and DNA repair genes.

To use a custom panel: `--germline_panel /path/to/my_panel.txt` (one gene symbol per line, `#` for comments).

## Citation

If you use `nf-sarek-clinical` in your research, please cite:

> [Your citation here]

Also cite the tools used:
- **nf-core/sarek**: Garcia M, et al. (2020). *Genome Medicine*.
- **dbNSFP**: Liu X, et al. (2020). *Genome Medicine*.
- **CIViC/civicpy**: Griffith M, et al. (2017). *Nature Genetics*.
- **SigProfilerAssignment**: Bergstrom EN, et al. (2022). *Cell Genomics*.
- **ASCAT**: Van Loo P, et al. (2010). *PNAS*.

## License

MIT — see [LICENSE](LICENSE).
