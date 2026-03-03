# nf-sarek-clinical

**Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exome outputs.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

## Overview

`nf-sarek-clinical` is a Nextflow DSL2 pipeline that takes the output of [nf-core/sarek](https://nf-co.re/sarek) (tumour–normal WES) and produces a clinical-grade molecular report (HTML + PDF).

| Module | Analysis | Output |
|--------|----------|--------|
| **Germline SNV** | dbNSFP annotation → 90-gene panel + HIGH/MODERATE filter + ClinVar benign exclusion + gnomAD rescue (Ensembl VEP) + VAF filter + multi-hit same-exon filter → ACMG P/LP/VUS classification with 9-tool in silico panel + individual criteria codes | `*.germline_snvs.tsv`, `*.germline_all_filtered.tsv`, `*.germline_snvs_acmg_annotated.tsv` |
| **Pharmacogenomics** | DPYD/UGT1A1 rsID filter | `*.pharmaco.tsv` |
| **Germline CNV (GATK)** | Cohort PoN-based germline copy number calling + MANE annotation + quality filter (NP ≥ 5 AND QA ≥ 10) | `*_segments_annotated.tsv` |
| **Somatic SNV** | dbNSFP + somatic panel filter + CIViC annotation + AMP/ASCO/CAP tier classification | `*.somatic_snvs.tsv`, `*.civic.exploded.tsv` |
| **Somatic CNV (ASCAT/CNVkit)** | HOMDEL/HETDEL/AMP (38-gene panel) + LOH 2nd-hit + somatic biallelic event detection | `*_somatic_cnv_report.tsv` |
| **Somatic CNV (GATK)** | Allele-specific somatic copy-number with denoising PoN | `*_annotated.tsv` |
| **TMB** | PASS + impact + depth filter → mut/Mb | `*.tmb.tsv` |
| **Mutational Signatures** | SigProfilerMatrixGenerator (observed spectra) + SigProfilerAssignment (COSMIC v3.4) SBS96/DBS/ID | `matrix_profiles/` (observed matrices), assignment activities, 4 PNG plots |
| **HRD** | HRProfiler: SBS + indel + CNV48 features → logistic model | `*_hrd_predictions.tsv` |
| **Gene Coverage** | Per-gene mean/min/max from mosdepth | `*_gene_coverage.tsv` |
| **Germline–Somatic Correlation** | Two-hit analysis per gene: germline P/LP SNV or CNV + somatic second hit detection + somatic biallelic event detection + SBS and ID signature correlation (MMR: SBS6/15/20/21/26/44 + ID1/ID2; HRR: SBS3/SBS8) | `correlation_section.json`, `correlation_section.html` |
| **HTML + PDF Report** | Interactive clinical report + printable PDF with Major Clinical Findings cards | `*_molecular_report.html`, `*_clinical_report.pdf` |

---

## Running the pipeline

There are two ways to run depending on your starting point:

| | Option A — Full pipeline | Option B — Clinical report only |
|---|---|---|
| **Starts from** | Raw FASTQ files | Completed nf-core/sarek results |
| **Script** | `./run_full.sh` | `nextflow run main.nf` |
| **When to use** | Standard end-to-end workflow | sarek already ran separately |
| **Duration** | ~10–14 h (sarek) + ~1 h (clinical) | ~1 h |

---

### Option A — Full pipeline from FASTQ (recommended)

`run_full.sh` is the orchestrator included in this repo. It runs nf-core/sarek (alignment + variant calling) followed by nf-sarek-clinical (clinical report) in sequence, using lab-fixed parameters stored in `conf/sarek.config` and `full_pipeline.params.yml`.

#### Step 1 — Create the samplesheet

Create `/home/user/runs/<RUN_NAME>/samplesheet.csv` following nf-core/sarek format:

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
1229,XY,0,1229_N07,2,/home/user/reads/.../1229_N07_L2_1.fq.gz,/home/user/reads/.../1229_N07_L2_2.fq.gz
1229,XY,1,1229_T08,1,/home/user/reads/.../1229_T08_L1_1.fq.gz,/home/user/reads/.../1229_T08_L1_2.fq.gz
1229,XY,1,1229_T08,2,/home/user/reads/.../1229_T08_L2_1.fq.gz,/home/user/reads/.../1229_T08_L2_2.fq.gz
```

| Column | Description |
|--------|-------------|
| `patient` | Patient / run ID — shared across all samples of the same patient |
| `sex` | `XX` (female) or `XY` (male) — required by ASCAT |
| `status` | `0` = normal, `1` = tumor |
| `sample` | Unique sample name. Convention: `<patient>_N##` (normal), `<patient>_T##` (tumor) |
| `lane` | Lane number from the FASTQ filename (e.g. `L2` → `2`). Add one row per lane |
| `fastq_1` / `fastq_2` | Absolute paths to R1 and R2 FASTQ files |

> The normal sample name must match `.*_N\d+` (e.g. `1229_N07`) and the tumor `.*_T\d+` (e.g. `1229_T08`). The clinical pipeline will create the pair `1229_T08_vs_1229_N07`.

#### Step 2 — Launch

```bash
cd /home/user/nf-sarek-clinical

./run_full.sh \
  --run        <RUN_NAME> \
  --samplesheet /home/user/runs/<RUN_NAME>/samplesheet.csv
```

**Real example (run 20260225_871):**

```bash
cd /home/user/nf-sarek-clinical

./run_full.sh \
  --run        20260225_871 \
  --samplesheet /home/user/runs/20260225_871/samplesheet.csv
```

#### Common flags

```bash
# Resume both phases from last checkpoint (after failure or interruption)
./run_full.sh \
  --run        20260225_871 \
  --samplesheet /home/user/runs/20260225_871/samplesheet.csv \
  --resume

# Clinical report only (sarek already completed)
./run_full.sh --run 20260225_871 --skip-sarek

# Resume clinical report only
./run_full.sh --run 20260225_871 --skip-sarek --resume
```

#### Output layout

```
/home/user/runs/<RUN_NAME>/
├── samplesheet.csv          ← input (you create this)
├── run.config               ← auto-generated by run_full.sh
├── work/                    ← sarek work dir (can be deleted after completion)
├── results/                 ← nf-core/sarek outputs (BAMs, VCFs, CNVs, QC)
├── clinical_work/           ← nf-sarek-clinical work dir
└── clinical_results/        ← HTML + PDF clinical reports  ← final output
```

#### Configuration files (edit when lab resources change)

| File | What it controls |
|------|-----------------|
| `conf/sarek.config` | Fixed sarek params: genome, tools, capture kit BED, ASCAT resources |
| `full_pipeline.params.yml` | Fixed sarek-clinical paths: dbNSFP, SnpSift, CNV models, reference FASTAs |

---

### Option B — Clinical report only (sarek already completed)

Use this when nf-core/sarek was run separately and you only need to generate the clinical report.

#### Step 1 — Install Nextflow and Docker

```bash
# Install Nextflow (requires Java 11+)
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/

# Verify
nextflow --version   # should be ≥ 23.04.0
docker --version     # should be ≥ 20.10
```

### Step 2 — Clone the pipeline

```bash
git clone https://github.com/YOUR_ORG/nf-sarek-clinical.git
cd nf-sarek-clinical
```

### Step 3 — Build and configure Docker images

Three custom images must be built and pushed to your registry:

```bash
REGISTRY="ghcr.io/YOUR_ORG/nf-sarek-clinical"

docker build -f docker/Dockerfile.civicpy    -t ${REGISTRY}/civicpy:3.1.0      docker/
docker build -f docker/Dockerfile.sigprofiler -t ${REGISTRY}/sigprofiler:0.1.7  docker/
docker build -f docker/Dockerfile.r_annotation -t ${REGISTRY}/r-annotation:1.0.0 docker/

docker push ${REGISTRY}/civicpy:3.1.0
docker push ${REGISTRY}/sigprofiler:0.1.7
docker push ${REGISTRY}/r-annotation:1.0.0
```

Then edit `conf/docker.config` and replace `ghcr.io/your-org` with your registry prefix:

```groovy
withLabel:python      { container = 'ghcr.io/YOUR_ORG/nf-sarek-clinical/civicpy:3.1.0'      }
withLabel:civicpy     { container = 'ghcr.io/YOUR_ORG/nf-sarek-clinical/civicpy:3.1.0'      }
withLabel:sigprofiler { container = 'ghcr.io/YOUR_ORG/nf-sarek-clinical/sigprofiler:0.1.7'  }
withLabel:r_base      { container = 'ghcr.io/YOUR_ORG/nf-sarek-clinical/r-annotation:1.0.0' }
```

The following images are pulled automatically from public registries:

| Label | Image |
|-------|-------|
| `snpsift` | `quay.io/biocontainers/snpeff:5.1--hdfd78af_2` |
| `bcftools` | `quay.io/biocontainers/bcftools:1.20--h8b25389_0` |
| `gatk` | `broadinstitute/gatk:4.6.2.0` |
| `hrprofiler` | `florpio/hrprofiler:1.0-full` |

### Step 4 — Prepare reference databases

#### dbNSFP (required)

```bash
# Download from https://sites.google.com/site/jpopgen/dbNSFP (~30 GB)
bgzip dbNSFP4.x.txt
tabix -s 1 -b 2 -e 2 dbNSFP4.x.txt.gz
```

#### SnpSift (required)

```bash
wget https://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
# SnpSift.jar is at snpEff/SnpSift.jar
```

#### GATK Germline gCNV assets (required when `--skip_germline_cnv false`)

These are built **once** from a cohort of normal samples. See the [GATK gCNV tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035531092).

| Asset | Created by |
|-------|-----------|
| `exomas.interval_list` | `FilterIntervals` on normals |
| `ploidy-model/` | `DetermineGermlineContigPloidy --run-mode COHORT` |
| `PoN_CNV/` | `GermlineCNVCaller --run-mode COHORT` |
| `hg38.fa` + `hg38.dict` | Full hg38 reference FASTA |
| `exons_mane.txt` | MANE transcript exon coordinates |
| `exon_annotation.R` | Annotation R script |

> **Important:** The interval list used for patient `CollectReadCounts` must exactly match the one used to build the ploidy model and PoN. All must use the same GATK version (`broadinstitute/gatk:4.6.2.0`).

#### GATK Somatic CNV assets (required when `--skip_somatic_cnv_gatk false`)

| Asset | Created by |
|-------|-----------|
| `wes-do-gc.pon.hdf5` | `CreateReadCountPanelOfNormals` |
| `preprocessed.interval_list` | `FilterIntervals` on normals |
| `common_snps_preprocessed.vcf.gz` | Common gnomAD SNPs (bgzipped + tabix) |
| `hg38.main_chroms.fa` + `.dict` | Main-chromosomes-only reference FASTA |
| `annotate_cnvs_complete_fixed.R` | CNV annotation R script |

### Step 5 — Create the run directory

Organise each sequencing run in its own directory:

```
/home/user/runs/MY_RUN/
├── results/          # nf-core/sarek output directory  ← --sarek_outdir
├── clinical_results/ # nf-sarek-clinical output        ← --outdir        (created by pipeline)
└── clinical_work/    # Nextflow work directory          ← -work-dir       (created by pipeline)
```

#### Required sarek output structure inside `results/`

The pipeline reads the following files from `--sarek_outdir`. All paths are relative to that directory. `{NORMAL}` and `{TUMOR}` are the sample names; `{PAIR}` = `{TUMOR}_vs_{NORMAL}`.

```
{sarek_outdir}/
│
├── annotation/
│   ├── haplotypecaller/
│   │   └── {NORMAL}/
│   │       ├── {NORMAL}.haplotypecaller.filtered_snpEff.ann.vcf.gz      ← germline SNVs
│   │       └── {NORMAL}.haplotypecaller.filtered_snpEff.ann.vcf.gz.tbi
│   └── mutect2/
│       └── {PAIR}/
│           ├── {PAIR}.mutect2.filtered_snpEff.ann.vcf.gz                ← somatic SNVs
│           └── {PAIR}.mutect2.filtered_snpEff.ann.vcf.gz.tbi
│
├── preprocessing/
│   └── recalibrated/
│       ├── {NORMAL}/
│       │   ├── {NORMAL}.recal.bam     ← used for GATK germline gCNV
│       │   └── {NORMAL}.recal.bam.bai
│       └── {TUMOR}/
│           ├── {TUMOR}.recal.bam      ← used for GATK somatic CNV
│           └── {TUMOR}.recal.bam.bai
│
├── variant_calling/
│   ├── ascat/
│   │   └── {PAIR}/
│   │       ├── {PAIR}.metrics.txt     ← purity, ploidy, LOH, WGD, GI  [required for reports]
│   │       └── {PAIR}*segment*.txt    ← ASCAT segments  [required for HRD]
│   ├── cnvkit/
│   │   └── {TUMOR}/
│   │       └── {TUMOR}.somatic.call.cns   ← CNVkit for somatic CNV interpretation
│   └── msisensorpro/
│       └── {PAIR}                     ← MSI value (file with no extension)
│
└── reports/
    └── mosdepth/
        ├── {NORMAL}/
        │   └── {NORMAL}.mosdepth.summary.txt   ← normal coverage
        └── {TUMOR}/
            └── {TUMOR}.mosdepth.summary.txt    ← tumor coverage
```

> **Important:** GATK modules (germline gCNV and GATK somatic CNV) are executed **from scratch** using the BAMs — they do not need to be pre-computed. All other files must come from a complete nf-core/sarek run with the following tools: `haplotypecaller`, `mutect2`, `ascat`, `cnvkit`, `msisensorpro`, `snpeff`.

### Step 6 — Run the clinical pipeline (Option B)

Full command with all modules active. Replace paths with your own installation:

```bash
cd /home/user/nf-sarek-clinical

nextflow run main.nf \
  -profile docker \
  -work-dir /home/user/runs/MY_RUN/clinical_work \
  \
  --sarek_outdir  /home/user/runs/MY_RUN/results \
  --run_name      MY_RUN \
  --outdir        /home/user/runs/MY_RUN/clinical_results \
  \
  --dbnsfp        /home/user/dbNSFP/dbNSFP4.1a.txt.gz \
  --snpsift_jar   /home/user/snpEff/SnpSift.jar \
  --panel_bed     /home/user/supporting_files/hg38/agilent/S07604514_hs_hg38/S07604514_Regions.bed \
  \
  --gcnv_ploidy_model   /home/user/supporting_files/cnv/ploidy-model \
  --gcnv_pon_dir        /home/user/supporting_files/cnv/PoN_CNV \
  --gcnv_interval_list  /home/user/supporting_files/cnv/exomas.interval_list \
  --gcnv_fasta          /home/user/supporting_files/hg38/hg38.fa \
  --gcnv_dict           /home/user/supporting_files/hg38/hg38.dict \
  --gcnv_mane_txt       /home/user/supporting_files/cnv/exons_mane.txt \
  --gcnv_exon_r         /home/user/supporting_files/cnv/exon_annotation.R \
  \
  --somatic_pon_hdf5       /home/user/supporting_files/cnv/wes-do-gc.pon.hdf5 \
  --somatic_interval_list  /home/user/supporting_files/cnv/preprocessed.interval_list \
  --common_snps_vcf        /home/user/supporting_files/cnv/common_snps_preprocessed.vcf.gz \
  --somatic_fasta          /home/user/supporting_files/cnv/main_chroms/hg38.main_chroms.fa \
  --somatic_dict           /home/user/supporting_files/cnv/main_chroms/hg38.main_chroms.dict \
  --somatic_annotate_r     /home/user/supporting_files/cnv/annotate_cnvs_complete_fixed.R
```

#### Modules that can be skipped when data is missing

| Flag | When to use |
|------|---------------|
| `--skip_somatic_cnv true` | No ASCAT/CNVkit output in `variant_calling/` |
| `--skip_hrd true` | No ASCAT segments available (required by HRProfiler) |
| `--skip_germline_cnv true` | No normal BAMs available |
| `--skip_somatic_cnv_gatk true` | No tumor BAMs available |
| `--skip_signatures true` | Too few somatic variants for SigProfiler |
| `--skip_coverage true` | No mosdepth results in `reports/` |
| `--skip_pdf_report true` | Skip clinical PDF (HTML only) |

### Step 7 — Resume after failure

All completed tasks are cached by Nextflow. After fixing an error, simply add `-resume`:

```bash
nextflow run main.nf -resume [all other params as above...]
```

If a task keeps failing despite fixes (stale work directory), remove its work dir:

```bash
rm -rf /home/user/runs/MY_RUN/clinical_work/XX/XXXXXX...
# Then resume
nextflow run main.nf -resume [all other params...]
```

---

## Minimal run (SNV only)

To run just germline and somatic SNV annotation (skip CNV, signatures, HRD):

```bash
nextflow run main.nf \
  -profile docker \
  --sarek_outdir  /home/user/runs/MY_RUN/results \
  --run_name      MY_RUN \
  --outdir        /home/user/runs/MY_RUN/clinical_results \
  --dbnsfp        /home/user/dbNSFP/dbNSFP4.1a.txt.gz \
  --snpsift_jar   /home/user/snpEff/SnpSift.jar \
  --skip_germline_cnv    true \
  --skip_somatic_cnv_gatk true \
  --skip_signatures      true \
  --skip_hrd             true
```

---

## Requirements

| Tool | Version |
|------|---------|
| **Nextflow** | ≥ 23.04.0 |
| **Docker** | ≥ 20.10 |
| RAM | ≥ 32 GB (GATK gCNV peaks at ~56 GB) |
| Disk | ≥ 200 GB (reference databases + work dir) |

**Sarek prerequisite:** nf-core/sarek ≥ 3.3.0 completed with tools: `haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff`

---

## Parameters

### Required for all runs

| Parameter | Description |
|-----------|-------------|
| `--sarek_outdir` | Path to the completed nf-core/sarek results directory |
| `--run_name` | Run identifier (used in filenames and report header) |
| `--outdir` | Where clinical results will be written |
| `--dbnsfp` | Path to `dbNSFP4.x.txt.gz` (bgzipped + tabix indexed) |
| `--snpsift_jar` | Path to `SnpSift.jar` |

### SNV filtering

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--germline_panel` | `assets/gene_panels/germline_90_cancer_predisposition.txt` | 90-gene hereditary cancer panel |
| `--somatic_snv_panel` | `assets/gene_panels/somatic_snv_panel.txt` | Somatic SNV gene list |
| `--germline_min_dp` | `5` | Minimum depth for germline variants |

### TMB

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--panel_bed` | — | WES capture BED for Mb calculation (required for TMB) |
| `--tmb_min_dp` | `10` | Min tumor depth |
| `--tmb_min_af` | `0.05` | Min allele fraction |
| `--tmb_min_alt_reads` | `10` | Min ALT reads |

### HRD

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hrd_organ` | `BREAST` | Organ model: `BREAST` or `OVARIAN` |
| `--hrd_threshold` | `0.45` | HRD probability threshold (≥ 0.45 → HRD+, < 0.45 → HRD−) |

### GATK Germline gCNV (required when `--skip_germline_cnv false`)

| Parameter | Description |
|-----------|-------------|
| `--gcnv_interval_list` | Preprocessed interval list (must match PoN build) |
| `--gcnv_ploidy_model` | Pre-built ploidy model directory |
| `--gcnv_pon_dir` | Pre-built gCNV PoN directory |
| `--gcnv_fasta` | Full hg38 reference FASTA |
| `--gcnv_dict` | hg38 sequence dictionary |
| `--gcnv_mane_txt` | MANE transcript exon coordinates |
| `--gcnv_exon_r` | Exon annotation R script |

### GATK Somatic CNV (required when `--skip_somatic_cnv_gatk false`)

| Parameter | Description |
|-----------|-------------|
| `--somatic_pon_hdf5` | Somatic PoN HDF5 |
| `--somatic_interval_list` | Preprocessed somatic interval list |
| `--common_snps_vcf` | Common SNPs VCF (bgzipped + tabix) |
| `--somatic_fasta` | Main-chromosomes reference FASTA |
| `--somatic_dict` | Sequence dictionary for main-chrom FASTA |
| `--somatic_annotate_r` | CNV annotation R script |

### Module toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_germline_snv` | `false` | Skip germline SNV postprocessing |
| `--skip_germline_cnv` | `false` | Skip GATK germline gCNV workflow |
| `--skip_somatic_snv` | `false` | Skip somatic SNV postprocessing |
| `--skip_somatic_cnv` | `false` | Skip ASCAT + CNVkit interpretation |
| `--skip_somatic_cnv_gatk` | `false` | Skip GATK somatic CNV workflow |
| `--skip_tmb` | `false` | Skip TMB calculation |
| `--skip_signatures` | `false` | Skip mutational signature analysis |
| `--skip_hrd` | `false` | Skip HRD (HRProfiler) |
| `--skip_pharmaco` | `false` | Skip pharmacogenomics |
| `--skip_coverage` | `false` | Skip gene coverage |
| `--skip_report` | `false` | Skip HTML report generation |
| `--skip_pdf_report` | `false` | Skip clinical PDF report |

---

## Output

```
clinical_results/
├── germline/
│   ├── snv/extract/         *.germline_snvs.tsv                    — final germline SNV table (all filters passed)
│   │                        *.germline_all_filtered.tsv             — ALL rejected variants with FILTER_REASON column
│   │                        *.germline_snvs_acmg_annotated.tsv      — all candidate variants with ACMG_class, ACMG_criteria, insil_dam/ben/total columns
│   ├── pharmaco/            *.pharmaco.tsv               — DPYD/UGT1A1 variants
│   └── cnv_gatk/annotated/  *_segments_annotated.tsv     — MANE-annotated germline CNV
├── somatic/
│   ├── snv/extract/         *.somatic_snvs.tsv           — final somatic SNV table
│   ├── snv/civic_extract/   *.civic.exploded.tsv         — CIViC evidence per variant
│   ├── cnv/                 *_somatic_cnv_report.tsv     — ASCAT/CNVkit CNV events
│   ├── cnv_gatk/annotated/  *_annotated.tsv              — GATK somatic CNV segments
│   ├── tmb/                 *.tmb.tsv                    — TMB calculation
│   ├── signatures/
│   │   ├── matrix_profiles/output/  SBS/, ID/, DBS/ observed count matrices
│   │   └── assignment/              COSMIC v3.4 assignment results + 4 PNG plots per context
│   └── hrd/                 *_hrd_predictions.tsv        — HRD probability + call
├── qc/gene_coverage/        *_gene_coverage.tsv          — per-gene coverage stats
└── reporting/
    ├── *_molecular_report.html    — interactive HTML clinical report
    ├── *_molecular_report.txt     — plain-text summary
    └── *_clinical_report.pdf      — printable clinical PDF report
```

---

## SNV Variant Filtering

### Germline Variants

**Input file** (nf-core/sarek output):
```
{sarek_outdir}/annotation/haplotypecaller/{NORMAL}/{NORMAL}.haplotypecaller.filtered_snpEff.ann.vcf.gz
```
> Already filtered by HaplotypeCaller and annotated by SnpEff with functional effects (`ANN[*]`).

**Step 1 — dbNSFP annotation**

SnpSift dbnsfp adds pathogenicity and population frequency fields to each variant:

| Category | Fields added |
|-----------|----------------|
| Pathogenicity | SIFT, SIFT4G, PolyPhen2 (HDIV/HVAR), LRT, MutationTaster, MutationAssessor, FATHMM, PROVEAN, REVEL, MutPred, PrimateAI, ClinPred, CADD, DANN, fathmm-MKL |
| Conservation | GERP++, phyloP100way, phastCons100way |
| Population | gnomAD exomes AF, gnomAD genomes AF |
| ClinVar | ID, clinical significance, disease, review, HGVS, MedGen/OMIM/Orphanet IDs |
| Function | Interpro domain |

**Step 1.5 — gnomAD AF rescue** (applied before Step 2)

For indels (and any variant missing gnomAD AF from dbNSFP), the pipeline queries the [Ensembl VEP REST API](https://rest.ensembl.org/documentation/info/vep_region_get) to retrieve the gnomAD population allele frequency. This is necessary because dbNSFP has limited indel coverage: many common insertions/deletions have no AF entry, which would otherwise cause them to spuriously pass the BA1/BS1 filters.

The rescue selects the highest AF among all colocated gnomAD variants within ±5 bp (to account for left-normalisation differences between the VCF representation and gnomAD). If the API returns no colocated variant, the AF is left as `0.0`.

**Step 2 — Panel and quality filter** (`SNPSIFT_FILTER_GERMLINE`)

Only variants meeting **all** of the following conditions are retained:

```
( ANN[*].GENE in germline_panel )
AND ( ANN[*].IMPACT = 'HIGH' OR ANN[*].IMPACT = 'MODERATE' )
AND ( DP >= germline_min_dp )
```

| Criterion | Default value | Parameter |
|----------|-------------------|-----------|
| Gene panel | 90 hereditary cancer genes | `--germline_panel` |
| Minimum impact | HIGH or MODERATE (SnpEff) | — |
| Minimum depth (DP) | 5 reads | `--germline_min_dp` |

**Step 3 — Post-filtering on extraction** (applied after `SnpSift extractFields`)

Three sequential filters are applied to the extracted TSV. All rejected variants are collected with a `FILTER_REASON` column in `*.germline_all_filtered.tsv` for manual review.

**3a. ClinVar Benign filter** (AWK):
Variants with `clinvar_clnsig` containing `benign` or `likely_benign` are removed, unless they also contain `pathogenic`.
Reason tag: `clinvar_benign (clnsig=<value>)`

**3b. Heterozygosity / VAF filter** (AWK):
VAF is computed from `GEN[0].AD` as `alt / (ref + alt)`. Thresholds depend on impact:

| Condition | Threshold | Reason tag |
|-----------|-----------|------------|
| VAF below minimum | HIGH impact: VAF < 0.25 · MODERATE: VAF < 0.35 | `vaf_too_low (VAF=X.XXX, impact=Y)` |
| VAF too high (likely homozygous / artefact) | VAF > 0.75 | `vaf_too_high (VAF=X.XXX)` |

Variants with missing AD pass unchanged. HIGH-impact variants use a lower minimum VAF (0.25) to preserve potential PVS1-class calls with lower read support.

**3c. Multi-hit same-exon filter** (Python):
Flags variants that co-occur in the same gene AND in the same exon (identical `ANN[*].RANK`) OR within 30 bp of another variant in the same gene. Within each such cluster, a variant is removed if:
- Its `ANN[*].EFFECT` is `disruptive_inframe_deletion` or `conservative_inframe_deletion` (regardless of VAF), **OR**
- Its VAF is < 0.40 (unbalanced allele fraction suggestive of sequencing noise)

**Exemptions** — the following are always kept and never filtered by this step:
- Variants with `ANN[*].IMPACT = HIGH` (PVS1-level evidence)
- Variants with ClinVar `clinvar_clnsig` containing `pathogenic` (P or LP)

Reason tags:
- `multi_hit_same_exon_inframe_deletion (effect=<effect>, gene=<G>, exon=<E>)`
- `multi_hit_same_exon_low_vaf (VAF=X.XXX, gene=<G>, exon=<E>)`

**Output files:**
- `*.germline_snvs.tsv` — clean variants passing all filters (used in report)
- `*.germline_all_filtered.tsv` — all rejected variants with `FILTER_REASON` column

---

### Somatic Variants

**Input file** (nf-core/sarek output):
```
{sarek_outdir}/annotation/mutect2/{PAIR}/{PAIR}.mutect2.filtered_snpEff.ann.vcf.gz
```
> Already filtered by Mutect2 (FilterMutectCalls) and annotated by SnpEff.

**Step 1 — dbNSFP annotation**

Same set of fields as germline (SIFT, PolyPhen, REVEL, gnomAD, ClinVar, etc.).

**Step 2 — Panel and quality filter** (`SNPSIFT_FILTER_SOMATIC`)

Only variants meeting **all** of the following conditions are retained:

```
( FILTER = 'PASS' OR FILTER = 'haplotype' )
AND ( ANN[*].GENE in somatic_panel )
AND ( ANN[*].IMPACT = 'HIGH' OR ANN[*].IMPACT = 'MODERATE' )
AND ( GEN[1].DP >= tmb_min_dp )
```

| Criterion | Default value | Parameter |
|----------|-------------------|-----------|
| FILTER status | PASS or haplotype | — |
| Gene panel | Somatic oncology panel | `--somatic_snv_panel` |
| Minimum impact | HIGH or MODERATE (SnpEff) | — |
| Minimum tumor depth (GEN[1].DP) | 10 reads | `--tmb_min_dp` |

> **Note**: `GEN[1]` corresponds to the tumor sample in Mutect2 paired VCFs (tumor-normal).

**Step 3 — CIViC annotation**

`civicpy annotate-vcf` adds the `INFO/CIVIC` field with therapeutic evidence (evidence level A–E, clinical significance, diseases, therapies, etc.). If annotation fails, the original VCF passes through unchanged.

**Step 4 — TSV extraction**

`SnpSift extractFields` generates the final TSV. **No additional filters are applied** in this step.

---

### Summary comparison

| Criterion | Germline | Somatic |
|----------|----------|----------|
| Input VCF | HaplotypeCaller (already filtered) | Mutect2 FilterMutectCalls (already filtered) |
| Gene panel | 90 hereditary cancer genes | Somatic oncology panel |
| Impact | HIGH / MODERATE | HIGH / MODERATE |
| VCF FILTER field | Not applied (HaplotypeCaller filters by default) | PASS or haplotype only |
| Minimum depth | DP ≥ 5 (normal) | GEN[1].DP ≥ 10 (tumor) |
| ClinVar filter | Removes Benign/Likely_benign | Not applied |
| VAF filter | HIGH: 0.25–0.75 · MODERATE: 0.35–0.75 | Not applied |
| Multi-hit same-exon filter | VAF < 0.40 or inframe deletion in cluster (P/LP + HIGH exempt) | Not applied |
| Additional annotation | dbNSFP | dbNSFP + CIViC |

---

## ACMG Germline Variant Classification

All germline variants that pass the SNV filters are classified following ACMG/AMP 2015 guidelines using a point-based evidence scoring system. The classification logic is implemented in three places that must stay in sync:

| File | Used by |
|------|---------|
| `bin/final_report_acmg.py` | HTML report (`modules/local/final_report/main.nf`) |
| `modules/local/snpsift_germline/main.nf` Step 4 (inline Python) | ACMG-annotated TSV (`*.germline_snvs_acmg_annotated.tsv`) |
| `bin/clinical_pdf_report.py` — `classify_acmg_with_criteria()` | Clinical PDF report |
| `bin/germline_somatic_correlate.py` — `acmg_class()` | Germline–somatic correlation analysis |

### Classification priority

1. **ClinVar annotation** (highest priority): if `dbNSFP_clinvar_clnsig` is present and unambiguous, the ClinVar classification is used directly and overrides all computed scores.

| ClinVar value (case-insensitive) | ACMG class | Criteria code |
|----------------------------------|------------|---------------|
| contains `pathogenic` (not `likely`, not `conflicting`) | P | `ClinVar:P` |
| contains `likely_pathogenic` | LP | `ClinVar:LP` |
| contains `uncertain` | VUS | `ClinVar:VUS` |
| contains `likely_benign` | LB | `ClinVar:LB` |
| contains `benign` (not `pathogenic`) | B | `ClinVar:B` |
| contains `conflicting` (e.g. `Conflicting_interpretations_of_pathogenicity`) | — | Falls through to evidence scoring |

> **Important:** `Conflicting_interpretations_of_pathogenicity` contains the substring `pathogenic`, but it is **not** classified as P. It falls through to score-based classification because `conflicting` is checked first.

2. **Population frequency fast-track**: if gnomAD exomes AF > 0.05, the variant is classified as **B** with criterion `BA1` (standalone benign — very high frequency), bypassing all further scoring. The AF used here is the gnomAD-rescued value from Step 1.5 when the dbNSFP field is missing.

3. **Evidence scoring** (applied when ClinVar is absent/ambiguous and gnomAD AF ≤ 0.05):

#### REVEL-based in silico criteria

| Criterion | Code | Score | Condition |
|-----------|------|-------|-----------|
| Very strong pathogenic | `PVS1` | +8 | `ANN[*].IMPACT = HIGH` (frameshift, nonsense, splice) |
| Moderate pathogenic | `PM2` | +2 | gnomAD exomes AF < 0.001 (absent/very rare in population) |
| Strong pathogenic in-silico | `PP3_St` | +4 | REVEL score ≥ 0.932 |
| Moderate pathogenic in-silico | `PP3` | +2 | REVEL score 0.644–0.931 |
| Supporting pathogenic in-silico | `PP3_Sup` | +1 | REVEL score 0.500–0.643 |
| Strong benign (common variant) | `BS1` | −4 | gnomAD exomes AF > 0.01 |
| Strong benign in-silico | `BP4_St` | −4 | REVEL score < 0.150 |
| Moderate benign in-silico | `BP4` | −2 | REVEL score 0.150–0.289 |

#### 9-tool in silico majority vote

In addition to REVEL, the pipeline scores an ensemble vote across **9 pathogenicity prediction tools** from dbNSFP. For each tool, the first prediction is used (dbNSFP stores multi-transcript values comma-separated; only the primary transcript is read):

| Tool | Damaging call | Benign call |
|------|--------------|------------|
| SIFT | D | T |
| SIFT4G | D | T |
| PolyPhen2-HDIV | D, P | B |
| PolyPhen2-HVAR | D, P | B |
| LRT | D | N |
| MutationTaster | A, D | N, P |
| MutationAssessor | H, M | L, N |
| FATHMM | D | T |
| PROVEAN | D | N |

The majority vote result is added on top of the REVEL score:

| Condition | Code | Score |
|-----------|------|-------|
| Damaging ≥ 75% of voting tools | `PP3(dam/total)` | +2 |
| Damaging majority < 75% | `PP3_Sup(dam/total)` | +1 |
| Benign ≥ 75% of voting tools | `BP4(ben/total)` | −2 |
| Benign majority < 75% | `BP4_Sup(ben/total)` | −1 |

Tools with missing/unknown values (`.`, `NA`) are abstentions and do not count. If no tools vote, this step contributes 0.

4. **Final score → class mapping**:

| Score | ACMG Class |
|-------|-----------|
| ≥ 10 | **P** (Pathogenic) |
| 6–9 | **LP** (Likely Pathogenic) |
| −5 to 5 | **VUS** (Variant of Uncertain Significance) |
| ≤ −6 | **LB** (Likely Benign) |

### Criteria display in reports

Each variant in the germline SNV table shows a **Criteria** column with colour-coded badges for every evidence code that fired:

| Badge colour | Codes | Meaning |
|-------------|-------|---------|
| Red | `PVS1`, `PM2`, `PP3_St`, `PP3`, `PP3_Sup`, `PP3(n/m)`, `PP3_Sup(n/m)` | Pathogenic evidence |
| Blue | `BA1`, `BS1`, `BP4_St`, `BP4`, `BP4(n/m)`, `BP4_Sup(n/m)` | Benign evidence |
| Purple | `ClinVar:P`, `ClinVar:LP`, `ClinVar:VUS`, `ClinVar:LB`, `ClinVar:B` | Direct ClinVar classification |

Example for a HIGH-impact variant absent from gnomAD with 7/9 tools damaging: **PVS1 · PM2 · PP3(7/9)** (score = 8+2+2 = 12 → P)

---

## Report Filtering Criteria by Section

This section documents **exactly which records appear in each section of the HTML and PDF reports** and which are excluded.

### 1. Germline SNVs (ACMG Classification)

**Source file:** `clinical_results/germline/snv/extract/*.germline_snvs.tsv`

**Displayed:** Only variants classified as **P, LP, or VUS** are shown. B and LB variants are excluded from the report (they passed all pipeline filters but are not clinically actionable).

| Filter | Kept | Excluded |
|--------|------|---------|
| ACMG class | P, LP, VUS | B (Benign), LB (Likely Benign) |

**Columns shown:** Gene · HGVS_c · HGVS_p · Impact · DP · VAF · REVEL · gnomAD · ClinVar · ACMG class · **Criteria codes**

**Criteria column:** Each row shows individual ACMG evidence codes (e.g. `PVS1 · PM2 · PP3`) colour-coded by direction (red = pathogenic, blue = benign, purple = ClinVar). If only a ClinVar code fires, it appears alone (e.g. `ClinVar:LP`).

**Interactive filter buttons** (HTML only): All variants · P + LP + VUS · P + LP only

---

### 2. Germline CNVs (GATK gCNV)

**Source file:** `clinical_results/germline/cnv_gatk/annotated/*_segments_annotated.tsv`

**Two-stage filter** — both conditions must be met simultaneously:

| Filter | Threshold | Clinical rationale |
|--------|-----------|-------------------|
| **NP** (number of exon capture intervals) | ≥ 5 | Segment must span ≥ 5 exon probes. Rules out single-exon artefacts and single-probe noise. |
| **QA** (ALT allele Phred quality) | ≥ 10 | ≥ 90% Phred-confidence in the DEL/DUP genotype call (10^(−10/10) = 10% maximum error rate). |

Additionally:
- Only **CN < 2** (deletions) are reported as first germline hits. Duplications (CN > 2) are excluded — they are not classic tumour suppressor gene inactivation events.
- CN = 0 is reported as **HOMDEL** (P); CN = 1 as **HETDEL** (LP).

**What is filtered out:**
- Segments with NP < 5 (single-probe or very small segments — high false-positive rate)
- Segments with QA < 10 (low-confidence genotype calls)
- CN ≥ 2 segments (diploid or copy gains)

**Example false-positive removed by these filters:**
```
Gene: CHEK2 | CN=0 | NP=2 | QA=2 | size=133 bp
→ Filtered: NP=2 < 5 AND QA=2 < 10
```

**Columns shown:** Locus · Gene(s) · Type (HETDEL/HOMDEL) · CN · NP · QA · QS · Quality

If no germline CNV segment passes both quality filters, a styled callout box is shown:
> **Germline CNV Analysis — No Reportable Events**: The GATK germline copy-number variant analysis did not identify any CN ≠ 2 events meeting quality thresholds (NP ≥ 5 intervals and QA ≥ 10) in the germline cancer predisposition panel.

---

### 3. Somatic SNVs (AMP/ASCO/CAP Tier Classification)

**Source file:** `clinical_results/somatic/snv/extract/*.somatic_snvs.tsv`

All variants in the somatic SNV file are displayed — no further post-extraction filtering is applied at the report level. The filtering was done upstream (Step 2 of somatic variant filtering: PASS/haplotype, somatic panel, HIGH/MODERATE impact, tumor DP ≥ 10).

**AMP/ASCO/CAP tiers** (shown as colour-coded badges):

| Tier | Colour | Criteria |
|------|--------|----------|
| IA | Red | FDA-approved therapy in the same cancer type (CIViC evidence level A/B) |
| IB | Red | FDA-approved therapy in a different cancer type or clinical trial |
| IIA | Orange | Well-powered clinical study in the same cancer type |
| IIB | Orange | Well-powered clinical study in a different cancer type |
| III | Yellow | Clinical evidence: case reports / small series |
| IV | Grey | Preclinical or biological evidence only |

Tier is assigned by cross-referencing with CIViC evidence. Variants without any CIViC match are displayed without a tier badge.

---

### 4. Somatic CNVs — ASCAT/CNVkit

**Source file:** `clinical_results/somatic/cnv/ascat_cnvkit/*_somatic_cnv_report.tsv`

**Displayed:** Only events in the **38-gene somatic oncogene/TSG panel** (`assets/gene_panels/somatic_38_oncogenes.txt`) meeting the following event-type criteria:

| Event type | CN | Interpretation |
|-----------|-----|---------------|
| HOMDEL | 0 | Complete gene loss (both alleles deleted) |
| HETDEL | 1 | Heterozygous deletion (one allele lost) |
| AMP | ≥ 5–6 (context-dependent) | High-level amplification |
| 2nd_hit_LOH | 2 (minor=0) | Copy-neutral LOH — somatic loss of wild-type allele |

Only the events listed above are shown. Diploid segments (CN = 2, minor ≥ 1) are excluded.

---

### 5. Somatic CNVs — GATK

**Source file:** `clinical_results/somatic/cnv_gatk/annotated/*_annotated.tsv`

All annotated segments are displayed — no event-type filter is applied. This section shows the raw allele-specific copy-number landscape for manual review.

---

### 6. Germline–Somatic Correlation Analysis

**Source:** `germline_somatic_correlate.py` — runs at report generation time using the upstream filtered files as input.

This section integrates germline variants with somatic second hits and pathway biomarkers to assess biallelic inactivation. It is split into two sub-tables.

#### 6a. Germline–Somatic Correlation table

**Input variants (first germline hit):**

| Source | Filter applied |
|--------|---------------|
| Germline SNVs | ACMG class = **P or LP only** (from `*.germline_snvs.tsv`). VUS variants are excluded from the correlation table — they appear only in the Germline SNV section. |
| Germline CNVs | CN < 2, **NP ≥ 5**, **QA ≥ 10** (same thresholds as the Germline CNV report section) |

**Second hit detection (somatic events):**
- Somatic SNV in the same gene (HIGH or MODERATE impact, from `*.somatic_snvs.tsv`)
- Somatic CNV/LOH in the same gene (HETDEL, HOMDEL, 2nd_hit_LOH from `*_somatic_cnv_report.tsv`)

**Pathway biomarker evaluation:**

| Pathway | Biomarker | Threshold | Active if |
|---------|-----------|-----------|-----------|
| HRR (BRCA1/2, PALB2, RAD51C/D, ATM, CHEK2, …) | HRD probability | ≥ 0.45 | HRD+ |
| HRR | SBS3 / SBS8 signatures | detected | any activity > 0 |
| MMR (MLH1, MSH2, MSH6, PMS2, EPCAM, …) | TMB | ≥ 10 mut/Mb | HIGH |
| MMR | MSI | ≥ 3.5% | HIGH |
| MMR | SBS6/15/20/21/26/44 signatures (SBS) | detected | any activity > 0 |
| MMR | ID1 / ID2 signatures (indel) | detected | any activity > 0 |

> The MMR biomarker evaluation includes both **SBS** (SBS6/15/20/21/26/44) and **ID** (ID1/ID2) canonical MMR signatures. These are stored in `detected_sigs` and `detected_id_sigs` respectively in the JSON output. A patient may have active ID MMR signatures without any SBS MMR signatures (or vice versa).

**Display filter (critical):** A germline variant is only included in this table if it has **at least one** of:
- A detected somatic second hit (`has_second_hit = True`), **OR**
- A positive pathway biomarker (`biomarker_match = True`)

Germline-only variants (no second hit, no biomarker concordance) are excluded. These are still reported in the Germline SNV section above but have no somatic actionability to show in a correlation context.

**Interpretation labels:**

| Label | Criteria | Badge colour |
|-------|----------|-------------|
| CONCORDANT | Biallelic (2nd hit) + pathway biomarkers active | Red |
| Biallelic — pathway biomarkers absent | 2nd hit detected, pathway inactive | Orange |
| Pathway active — 2nd hit undetected | No 2nd hit, but biomarkers elevated | Orange |
| Biallelic — Other TSG | 2nd hit in non-MMR/HRR gene | Orange |
| VUS — Evidence supports pathogenicity | VUS + 2nd hit + biomarker active | Orange |
| VUS — Insufficient somatic evidence | VUS without somatic confirmation | Yellow |
| Germline only | No somatic evidence | Grey — **not shown** |

#### 6b. Somatic Biallelic Events (two somatic hits, no germline anchor)

**Filter:** Genes that have **both**:
1. A somatic HIGH or MODERATE impact SNV (from `*.somatic_snvs.tsv`)
2. A somatic loss-of-function CNV event: **HETDEL, HOMDEL, or 2nd_hit_LOH** (from `*_somatic_cnv_report.tsv`)

**Exclusion:** Genes already reported in the Germline–Somatic Correlation table (6a) are excluded to avoid double-reporting.

**Excluded somatic CNV events:** Amplifications (AMP) do not qualify as a second inactivating hit for this analysis.

**Mechanism interpretation:** The report generates a clinical narrative for each event. For HETDEL events, the somatic SNV allele fraction (VAF) is used to infer allele-specific behaviour:
- VAF > 0.50: "*The elevated VAF is consistent with preferential loss of the wild-type allele*" — classic biallelic somatic inactivation pattern
- VAF ≤ 0.50: "*May represent biallelic inactivation depending on clonal architecture*"

**Example (run 20260225_871):**
```
TP53:  Somatic MODERATE p.Met237Ile (VAF=0.708) + Somatic HETDEL (CN=1)
       → Biallelic somatic inactivation: elevated VAF (>0.71) consistent with
         preferential loss of the wild-type allele.
```

---

## HTML & PDF Report Sections

### Major Clinical Findings (PDF only)

The PDF report opens with a **Major Clinical Findings** section immediately after the KPI dashboard. It renders colour-coded cards summarising the most actionable findings. Each card is only shown when its trigger conditions are met:

#### Finding 1 — HRD Pathway

**Trigger:** HRD probability ≥ 0.45 **OR** at least one germline P/LP variant in an HRR gene (BRCA1/2, PALB2, RAD51C/D, ATM, CHEK2, …).

Content shown:
- HRD Status (HRD+ / HRD−) + probability
- Each germline HRR hit: gene, ACMG class, HGVS_p/c, somatic 2nd hit(s)
- HRD signatures (SBS3/SBS8) detected or not detected

#### Finding 2 — MMR Deficiency

**Trigger:** MSI ≥ 3.5% **OR** germline P/LP MMR variant **OR** SBS MMR signature detected (SBS6/15/20/21/26/44) **OR** ID MMR signature detected (ID1/ID2) **OR** somatic biallelic MMR event (see 6b).

Content shown:
- MSI (%) + MSI status + TMB (mut/Mb)
- Each germline MMR hit: gene, ACMG class, HGVS_p/c, somatic 2nd hit(s)
- **Biallelic somatic** MMR events (e.g. somatic splice variant + cnLOH), shown with explicit `Biallelic somatic:` label
- SBS MMR signatures: listed if detected, "not detected" if absent
- ID MMR signatures: listed if detected (ID1/ID2), "not detected" if absent

> SBS and ID signatures are reported on **separate lines** so that a patient with ID2 (indel MMR signature) but no SBS MMR signatures still shows `ID MMR signatures: ID2` clearly, rather than "MMR signatures not detected".

#### Finding 3 — Double Somatic Hits

**Trigger:** Only genes that appear in the germline–somatic correlation results (with `biomarker_match = True`) **OR** in the somatic biallelic events (6b) AND have ≥ 2 somatic SNVs. Genes with two somatic SNVs but **no pathway biomarker support** are excluded to avoid noise.

#### Finding 4 — High-Confidence CNV Events

**Trigger:** At least one HIGH-confidence HOMDEL (CN=0) or high-level AMP in the harmonised somatic CNV data.

---

### KPI Dashboard

TMB (mut/Mb) · MSI · LOH · Genomic Instability · WGD · Ploidy · Tumor Purity · Coverage Normal/Tumor · HRD Probability · HRD Status

**MSI display logic:**

| MSIsensorPro value | Displayed label |
|--------------------|----------------|
| < 5% | `MSI not detected (MSS) — X.XX%` |
| ≥ 5% | `MSI detected — X.XX%` |

**HRD Status display logic:**

| HRProfiler probability | Status shown |
|------------------------|-------------|
| ≥ 0.45 | `HRD+` |
| < 0.45 | `HRD−` |

### Report Section Order

| # | Section | PDF page | Key filter applied |
|---|---------|----------|--------------------|
| 0 | **Major Clinical Findings** | Page 2 | HRD/MMR/CNV triggers (see Major Clinical Findings above) |
| 1 | **Germline–Somatic Correlation** | Page 3 | P/LP only; has_second_hit OR biomarker_match |
| 1b | **Somatic Biallelic Events** | Page 3 (sub-table) | somatic SNV + somatic LOH/DEL, not already in 1 |
| 2 | **Germline SNVs** | Page 4 | ACMG P/LP/VUS only |
| 3 | **Somatic SNVs** | Page 5 | All from filtered TSV (no additional filter) |
| 4 | **Somatic CNVs (ASCAT/CNVkit)** | Page 6 | 38-gene panel, HOMDEL/HETDEL/AMP/LOH |
| 5 | **TMB** | Page 7 | N/A |
| 6 | **Pharmacogenomics** | Page 8 | DPYD/UGT1A1 |
| 7 | **Gene Coverage** | Page 9 | All genes in panels |
| 8 | **HRD Predictions** | Page 10 | N/A |
| 9 | **Mutational Signatures** | Page 11 | SBS96 + ID83 spectrum + etiology plots |
| 10 | **Germline CNVs (GATK)** | Page 12 | CN < 2, NP ≥ 5, QA ≥ 10 |
| 11 | **Somatic CNVs (GATK)** | Page 13 | All annotated segments |

---

## Mutational Signatures

### Overview

The signatures module runs in two steps:

1. **SigProfilerMatrixGenerator** — generates per-sample observed mutation count matrices (SBS96, ID83, DBS78) from the filtered VCFs. Output is saved to `somatic/signatures/matrix_profiles/output/`.
2. **SigProfilerAssignment** — decomposes the observed spectrum against COSMIC v3.4 reference signatures (SBS-96, DBS, ID contexts). Output is in `somatic/signatures/assignment/`.

After SigProfilerAssignment completes, `bin/generate_signature_plots.py` is automatically called to produce four publication-quality PNG plots per run.

### Plots generated

| Plot | File | Content |
|------|------|---------|
| SBS96 matrix profile | `*_SBS96/…/SBS96_matrix_profile.png` | 96-channel trinucleotide substitution spectrum using **observed** counts (y-axis: "Mutations (observed)") |
| SBS96 etiology activity | `*_SBS96/…/etiology_activity_plot.png` | Horizontal bar chart per active signature, colour-coded with COSMIC v3.4 aetiology (Clock-like, APOBEC, HRD, MMR, POLE/POLD1, Tobacco, UV, etc.) |
| ID83 matrix profile | `*_ID/…/ID83_matrix_profile.png` | 83-channel indel spectrum using observed counts |
| ID etiology activity | `*_ID/…/etiology_activity_plot.png` | Etiology-annotated bar chart for indel signatures (ID1–ID23) |

### COSMIC Aetiology column

The HTML report's signature assignment tables include a **COSMIC Aetiology** column showing a colour-coded badge for each active signature. The full aetiology map covers all COSMIC v3.4 SBS (SBS1–SBS99) and ID (ID1–ID23) signatures.

### Minimum variant requirement

SigProfiler requires a minimum number of somatic variants to produce reliable decompositions. If the filtered VCF has very few PASS variants (< ~30 SNVs), the assignment may not detect any active signatures or may fail silently. In such cases, use `--skip_signatures true`.

---

## Germline–Somatic Integration Details

The pipeline performs automated two-hit analysis for germline variants in cancer predisposition genes. See [Report Filtering Criteria — Germline–Somatic Correlation](#6-germlinendashsomatic-correlation-analysis) for the full filtering specification.

**Pathway gene assignments:**

| Pathway | Genes |
|---------|-------|
| **HRR** | BRCA1, BRCA2, PALB2, RAD51C, RAD51D, ATM, CHEK2, BARD1, BRIP1, NBN, RAD50, MRE11, FANCA, FANCC, FANCD2, FANCE, FANCF, FANCG, FANCL, FANCM, CDK12, ABRAXAS1, ABRAXAS2 |
| **MMR** | MLH1, MSH2, MSH6, PMS2, EPCAM, MSH3, MLH3, PMS1 |
| **Other TSG** | All remaining genes in the germline panel (no specific pathway biomarker applies) |

**Germline CNV deletions** affecting known panel genes are treated as first germline hits (HETDEL → LP; HOMDEL → P) and searched for second hits identically to SNVs, using the same NP ≥ 5 and QA ≥ 10 quality filter.

---

## Gene Panels

Built-in panels in `assets/gene_panels/`:

- **`germline_90_cancer_predisposition.txt`** — 90 hereditary cancer predisposition genes (ACMG SF v3.2). Used for germline SNV filtering, germline CNV annotation, LOH 2nd-hit detection, and germline–somatic correlation.
- **`somatic_38_oncogenes.txt`** — 38 actionable somatic genes. Used for ASCAT/CNVkit HOMDEL/HETDEL/AMP reporting.
- **`somatic_snv_panel.txt`** — Somatic SNV panel covering oncogenes, TSGs, and DNA repair genes.

Custom panel: `--germline_panel /path/to/my_panel.txt` (one gene symbol per line).

---

## GATK Germline gCNV — cohort model

The germline CNV analysis runs in **CASE mode**: it processes each normal sample independently against a pre-built cohort model.

### Operational modes

| Mode | When to use | How |
|------|-------------|-----|
| **CASE (default)** | Routine — use existing model | Set `gcnv_ploidy_model` + `gcnv_pon_dir` as usual |
| **Build new cohort model** | New capture kit, cohort change, GATK version mismatch | Follow cohort protocol below, then re-run as CASE |
| **Skip gCNV** | No normal BAM available, or skipping CNV analysis | Add `--skip_germline_cnv true` to the run command |

> **Minimum cohort size**: ≥ 20 unrelated normal samples from the same sequencing platform and capture kit. More samples → better noise model. The current lab model (`ploidy-model/`, `PoN_CNV/`) was built with `broadinstitute/gatk:4.6.2.0` (gcnvkernel 0.9) and must not be mixed with other versions.

---

### Building a new cohort model (COHORT mode)

All steps below use the GATK Docker container directly. Run from a machine with Docker and at least 64 GB RAM.

#### Step 1 — PreprocessIntervals

```bash
docker run --rm \
  -v /home/user/supporting_files:/ref \
  broadinstitute/gatk:4.6.2.0 \
  gatk PreprocessIntervals \
    -R /ref/hg38/hg38.fa \
    -L /ref/hg38/agilent/S07604514_hs_hg38/S07604514_Regions_Annotated.bed \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O /ref/cnv/new_cohort/preprocessed_raw.interval_list
```

#### Step 2 — AnnotateIntervals (GC content)

```bash
docker run --rm \
  -v /home/user/supporting_files:/ref \
  broadinstitute/gatk:4.6.2.0 \
  gatk AnnotateIntervals \
    -R /ref/hg38/hg38.fa \
    -L /ref/cnv/new_cohort/preprocessed_raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    -O /ref/cnv/new_cohort/annotated_intervals.tsv
```

#### Step 3 — CollectReadCounts on all normal samples

```bash
mkdir -p /home/user/supporting_files/cnv/new_cohort/counts

for BAM in ${NORMALS_DIR}/*/results/preprocessing/recalibrated/*_N*/*.recal.bam; do
  SAMPLE=$(basename "$(dirname "$BAM")")
  docker run --rm \
    -v /home/user:/home/user \
    broadinstitute/gatk:4.6.2.0 \
    gatk CollectReadCounts \
      -I "$BAM" \
      -L "$INTERVAL_LIST" \
      --interval-merging-rule OVERLAPPING_ONLY \
      --format HDF5 \
      -O "${COUNTS_DIR}/${SAMPLE}.counts.hdf5"
done
```

#### Step 4 — FilterIntervals

```bash
HDF5_ARGS=$(ls /home/user/supporting_files/cnv/new_cohort/counts/*.hdf5 | \
            awk '{print "-I "$1}' | tr '\n' ' ')

docker run --rm \
  -v /home/user/supporting_files:/ref \
  broadinstitute/gatk:4.6.2.0 bash -c \
  "gatk FilterIntervals \
    -L /ref/cnv/new_cohort/preprocessed_raw.interval_list \
    --annotated-intervals /ref/cnv/new_cohort/annotated_intervals.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    ${HDF5_ARGS} \
    -O /ref/cnv/new_cohort/exomas.interval_list"
```

> This filtered interval list **must** be used for all subsequent steps. Save it as `gcnv_interval_list` in `full_pipeline.params.yml`.

#### Step 5 — DetermineGermlineContigPloidy (COHORT)

```bash
HDF5_ARGS=$(ls /home/user/supporting_files/cnv/new_cohort/counts/*.hdf5 | \
            awk '{print "--input "$1}' | tr '\n' ' ')

docker run --rm \
  -v /home/user/supporting_files:/ref \
  broadinstitute/gatk:4.6.2.0 bash -c \
  "gatk DetermineGermlineContigPloidy \
    -L /ref/cnv/new_cohort/exomas.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    ${HDF5_ARGS} \
    --contig-ploidy-priors /ref/cnv/new_cohort/contig_ploidy_priors.tsv \
    --output /ref/cnv/new_cohort/ \
    --output-prefix ploidy \
    --run-mode COHORT"
```

#### Step 6 — GermlineCNVCaller (COHORT)

```bash
HDF5_ARGS=$(ls /home/user/supporting_files/cnv/new_cohort/counts/*.hdf5 | \
            awk '{print "-I "$1}' | tr '\n' ' ')

docker run --rm \
  -v /home/user/supporting_files:/ref \
  broadinstitute/gatk:4.6.2.0 bash -c \
  "gatk GermlineCNVCaller \
    --run-mode COHORT \
    -L /ref/cnv/new_cohort/exomas.interval_list \
    ${HDF5_ARGS} \
    --contig-ploidy-calls /ref/cnv/new_cohort/ploidy-calls \
    --annotated-intervals /ref/cnv/new_cohort/annotated_intervals.tsv \
    --interval-merging-rule OVERLAPPING_ONLY \
    --output /ref/cnv/new_cohort/ \
    --output-prefix cnv_pon \
    --run-mode COHORT"
```

#### Step 7 — Register the new model

Update `full_pipeline.params.yml`:

```yaml
gcnv_interval_list: /home/user/supporting_files/cnv/new_cohort/exomas.interval_list
gcnv_ploidy_model:  /home/user/supporting_files/cnv/new_cohort/ploidy-model
gcnv_pon_dir:       /home/user/supporting_files/cnv/new_cohort/cnv_pon-model
```

---

## Troubleshooting

### Session lock error
```bash
rm -f /home/user/nf-sarek-clinical/.nextflow/cache/*/db/LOCK
```

### Stale work directory (task fails despite code fix)
```bash
# Find the work dir hash from the Nextflow error output, then delete it:
rm -rf /home/user/runs/MY_RUN/clinical_work/XX/XXXXXXXXXXXXXX
nextflow run main.nf -resume [params...]
```

### Out of memory
Edit `conf/base.config` to lower the `process_high` memory limit:
```groovy
withLabel:process_high {
    memory = { check_max( 40.GB * task.attempt, 'memory' ) }
}
```

### gCNV model version mismatch
The ploidy model and PoN must be built with the same GATK version as the pipeline (`broadinstitute/gatk:4.6.2.0`). Rebuild if created with a different version.

### CIViC annotation fails
CIViC requires internet access on first run. On failure the step passes through unannotated and continues.

### ACMG criteria not visible in PDF germline table
The ACMG classification logic exists in **three separate places** that must stay in sync:
1. `bin/final_report_acmg.py` — imported by the HTML report module
2. `modules/local/snpsift_germline/main.nf` Step 4 (inline Python heredoc) — generates `*.germline_snvs_acmg_annotated.tsv`
3. `bin/clinical_pdf_report.py` `classify_acmg_with_criteria()` — used by the PDF report directly

If ACMG scores differ between the HTML and PDF reports, check that all three are consistent (same ClinVar conflicting logic, same 9-tool panel, same comma-split normalization for multi-transcript dbNSFP fields).

---

## Citation

If you use `nf-sarek-clinical` in your research, please cite the following tools:

- **nf-core/sarek**: Garcia M, et al. (2020). *Genome Medicine*.
- **dbNSFP**: Liu X, et al. (2020). *Genome Medicine*.
- **CIViC/civicpy**: Griffith M, et al. (2017). *Nature Genetics*.
- **SigProfilerAssignment**: Bergstrom EN, et al. (2022). *Cell Genomics*.
- **HRProfiler**: Chopra N, et al. (2024). [github.com/FlorPio/HRProfiler-workflow](https://github.com/FlorPio/HRProfiler-workflow).
- **ASCAT**: Van Loo P, et al. (2010). *PNAS*.

## License

MIT — see [LICENSE](LICENSE).
