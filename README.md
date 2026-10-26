# nf-sarek-clinical

**Clinical interpretation pipeline for tumor–normal whole-exome sequencing, on top of nf-core/sarek.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

![Pipeline overview](metro_map_sarek_clinical.png)

---

## What this pipeline does

`nf-sarek-clinical` takes the output of [nf-core/sarek](https://nf-co.re/sarek) for
tumor–normal pairs (WES) and produces, for every pair, a laboratory HTML report and a
clinical PDF report built from one report model. With `run_full.sh` it runs both
phases: sarek (FASTQ → alignment, variant calling, copy number, MSI, QC) and this
pipeline (clinical interpretation → reports).

**Analyses**

| Area | What is done |
|---|---|
| Germline SNVs / indels | 105-gene hereditary cancer panel; MANE Select nomenclature (snpEff MANE database, MANE Plus Clinical shown as secondary); local gnomAD v4.1 (AF, FAF95, grpmax), ClinVar (with review stars), REVEL v1.3, SpliceAI; ACMG/AMP 2015 with ClinGen SVI points, gene-specific ClinGen VCEP rules for 17 genes, PS1 / PM5 from ClinVar protein coordinates |
| Germline CNVs | GATK gCNV in CASE mode against a laboratory cohort model; MANE exon annotation; low-mappability segments listed apart |
| Pharmacogenomics | DPYD and UGT1A1 genotypes → phenotype from both alleles (DPYD activity score, CPIC 2017; UGT1A1 decreased-function alleles, DPWG for irinotecan) → dose recommendation; variants outside the capture reported as not assessable |
| Somatic SNVs / indels | Mutect2 (report tables) and Mutect2 ∩ Strelka2 consensus (TMB, signatures, HRD); three gene tracks (core 60 genes, second-hit germline genes, TSO500 pan-tumor); AMP/ASCO/CAP tiers IA–IV from a frozen CIViC snapshot and the patient's tumor type; Tier IV for common gnomAD polymorphisms; TERT promoter assessability |
| Somatic CNVs | ASCAT + CNVkit (from sarek) and GATK somatic CNV (with the ASCAT purity), harmonised: HOMDEL, AMP (CN ≥ 5), second-hit LOH in the germline panel |
| Biomarkers | TMB (non-synonymous, capture ∩ MANE CDS), MSI (MSIsensor-pro), mutational signatures (SigProfilerAssignment, COSMIC v3.4 SBS / DBS / ID), HRD (HRProfiler), ASCAT purity / ploidy / WGD / LOH / GI with fit-quality alerts |
| Integration | Germline–somatic correlation (two hits in MMR / HRR genes, biallelic somatic events), major clinical findings, recruiting clinical trials (ClinicalTrials.gov) |
| QC | Per-gene coverage of normal and tumor; biomarkers "Indeterminate" when tumor coverage < 50X |

**Deliverables per tumor/normal pair** (`<outdir>/reporting/`)

- `<pair>_clinical_report.pdf` — clinical report: patient and specimen identification, every biomarker and the major findings on page 1, then germline P/LP/VUS, germline CNVs, somatic Tier I–III, signatures, somatic CNV / LOH (with the events left out stated), correlation, pharmacogenomics, trials, coverage, methodology and versions. English or Spanish (`--report_language`).
- `<pair>_molecular_report.html` — laboratory report: the whole model (every germline class with filters, every tier, Tier IV and TSO500 tables, all LOH events, DBS signatures, raw TMB and HRD tables, gene context per variant).
- `<pair>_molecular_report.txt` — plain-text summary.
- `<pair>.report_model.json` — the report model: every value and interpretation both reports show, with thresholds, versions and resources; the record of what the reports said.
- Structured TSV / VCF files of every step (see [Output](#output)).

The two reports never compute anything themselves: `REPORT_MODEL` (`bin/report_model.py`)
builds the model once and `bin/render_html.py` / `bin/clinical_pdf_report.py` render it,
so the HTML and the PDF cannot disagree.

---

## Quick Start

**1. Install Nextflow (≥ 24.04), Java 17 and Docker**, and Python 3 with `reportlab`
(and `pymupdf`) on the host — see [Requirements](#requirements).

**2. Get the code and build the custom images once**

```bash
git clone https://github.com/GENUIT-LBMU/nf-sarek-clinical.git && cd nf-sarek-clinical
bash docker/build_local.sh          # civicpy, r-annotation, sigprofiler images
```

**3. Build / download the reference resources** — see [Reference files](#reference-files).

**4. Create `lab_params.yml`** (in the repository directory, where `run_full.sh` looks for it)

```bash
cp full_pipeline.params.yml lab_params.yml   # replace every /path/to/ with your paths
```

**5. Write the samplesheet** (sarek format; the extra columns are optional and only
identify patient and specimens in the reports — see [Samplesheet](#samplesheet))

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
P001,XX,0,P001_N01,L001,/data/P001_N01_R1.fastq.gz,/data/P001_N01_R2.fastq.gz
P001,XX,1,P001_T01,L001,/data/P001_T01_R1.fastq.gz,/data/P001_T01_R2.fastq.gz
```

**6. Run**

```bash
# FASTQ → reports (sarek 3.5.1 + this pipeline)
./run_full.sh --run MY_RUN --samplesheet /path/to/samplesheet.csv

# sarek already run, results in $RUNS_BASE/MY_RUN/results
./run_full.sh --run MY_RUN --skip-sarek [--samplesheet /path/to/samplesheet.csv]

# resume either: add --resume
```

**7. Results** in `$RUNS_BASE/MY_RUN/clinical_results/reporting/` (`RUNS_BASE` defaults to `$HOME/runs`).

---

## Requirements

### Software

| Tool | Version | Notes |
|---|---|---|
| Nextflow | ≥ 24.04.0 | With Nextflow ≥ 26 set `NXF_SYNTAX_PARSER=v1` (the pipeline uses the classic script syntax); `run_full.sh` exports it |
| Java | 17 | For Nextflow |
| Docker | ≥ 20.10 | Or Singularity (`-profile singularity`) |
| Python 3 + `reportlab` | host | `REPORT_MODEL`, `FINAL_REPORT` and `CLINICAL_PDF_REPORT` run on the host Python (no container label); `pymupdf` (`fitz`) embeds the SigProfiler PDF plots |
| `bgzip` / `tabix` | host | Only for the `-stub` test |

Every bioinformatics tool runs in a container (`conf/containers.config`):

| Label / step | Image |
|---|---|
| SnpSift filter / extract, TMB | `quay.io/biocontainers/snpeff:5.1` |
| snpEff MANE re-annotation | `quay.io/biocontainers/snpeff:5.2` |
| SpliceAI | `quay.io/biocontainers/spliceai:1.3.1` |
| bcftools (consensus, gnomAD / ClinVar / REVEL annotation) | `quay.io/biocontainers/bcftools:1.20` |
| GATK (germline and somatic CNV) | `broadinstitute/gatk:4.6.2.0` (gcnvkernel 0.9) |
| R annotation (CNV, coverage) | `${docker_namespace}/r-annotation:1.0.0` |
| civicpy (CIViC annotation, somatic CNV interpretation) | `${docker_namespace}/civicpy:5.4.0` |
| SigProfiler | `${docker_namespace}/sigprofiler:0.1.8` |
| HRProfiler | `florpio/hrprofiler:1.0-full` |

The three `${docker_namespace}` images (default `nf-sarek-clinical`) are built locally
with `bash docker/build_local.sh`, or pulled from a registry with `--docker_namespace`
(`-profile mirror` uses digest-pinned `ghcr.io/genuit-lbmu/*` images; see
[`docker/README.md`](docker/README.md)).

**Network:** the reports query ClinicalTrials.gov (HTTPS) when they are built; without
access they print "ClinicalTrials.gov API unavailable" and the query terms. Nothing else
is downloaded at run time (gnomAD, ClinVar, REVEL and CIViC are local, versioned files).

### nf-core/sarek

`run_full.sh` runs nf-core/sarek **3.5.1** (`SAREK_REVISION`, overridable) with
`--tools haplotypecaller,mutect2,strelka,manta,cnvkit,ascat,msisensorpro,snpeff`,
from the run directory (so this repository's `nextflow.config` does not reach sarek),
and records the release in `results/pipeline_info/sarek_release.txt`. For your own sarek
runs, use the same tools: Strelka2 + Manta feed the consensus (without them each pair
falls back to Mutect2 PASS), ASCAT feeds purity / ploidy / LOH / HRD, MSIsensor-pro the
MSI score.

### Hardware

| Resource | Recommended | Notes |
|---|---|---|
| RAM | 64 GB | `process_high` (gCNV ploidy and caller) asks for 56 GB. On a smaller host set `--max_memory` (e.g. `32.GB`): every task is capped at `--max_cpus` / `--max_memory` / `--max_time` (`process.resourceLimits`) |
| CPU | 8–16 cores | `--max_cpus` default 16 |
| Disk | ≥ 200 GB | References + work directories (sarek work for one WES pair: ~100 GB) |

---

## First-time setup

### Step 1 — Clone the repository

```bash
git clone https://github.com/GENUIT-LBMU/nf-sarek-clinical.git /path/to/nf-sarek-clinical
cd /path/to/nf-sarek-clinical
```

### Step 2 — Params file

`run_full.sh` uses `<repository>/lab_params.yml` when it exists, otherwise
`full_pipeline.params.yml` (the template). Copy the template and set every path:

```yaml
# ── Annotation resources ─────────────────────────────────────────────────────
snpsift_jar:      /data/refs/snpEff/snpEff/SnpSift.jar
snpeff_mane_data: /data/refs/snpeff/data                         # snpEff MANE database
gnomad_tsv:       /data/refs/gnomad/gnomad_v4.1_panel.tsv.gz     # germline: required
clinvar_tsv:      /data/refs/clinvar/clinvar_20260822_panel.tsv.gz        # germline: required
revel_tsv:        /data/refs/revel/revel_v1.3_panel.tsv.gz       # germline: required
clinvar_protein_index: /data/refs/clinvar/clinvar_20260822_protein_index.tsv.gz  # PS1 / PM5
civic_cache:      /data/refs/civic/civicpy_cache_20261024.pkl    # somatic: required
tumor_type:       "Colorectal Cancer|Colon Cancer|Rectum Cancer|Colorectal Adenocarcinoma|Colon Adenocarcinoma|Rectum Adenocarcinoma"

# ── Capture and MANE ─────────────────────────────────────────────────────────
panel_bed:        /data/refs/hg38/agilent/S07604514_Regions.bed            # TMB territory, TERT note
gene_coords_bed:  /data/refs/hg38/agilent/S07604514_Padded.withGenes.bed   # somatic CNV interpretation
mane_cds_bed:     /data/refs/mane/cds_mane_v1.2.bed                        # TMB territory

# ── Germline GATK CNV (case mode against a cohort model) ────────────────────
gcnv_ploidy_model:  /data/refs/cnv/ploidy-model
gcnv_pon_dir:       /data/refs/cnv/PoN_CNV
gcnv_interval_list: /data/refs/cnv/exomas.interval_list
gcnv_fasta:         /data/refs/hg38/hg38.fa        # also: consensus (bcftools norm), SpliceAI
gcnv_dict:          /data/refs/hg38/hg38.dict
gcnv_mane_txt:      /data/refs/mane/exons_mane_v1.2.txt
gcnv_lowmap_bed:    /data/refs/cnv/lowmap/GRCh38_lowmappabilityall_panel.bed.gz

# ── Somatic GATK CNV ─────────────────────────────────────────────────────────
somatic_pon_hdf5:        /data/refs/cnv/wes-do-gc.pon.hdf5
somatic_interval_list:   /data/refs/cnv/preprocessed.interval_list
common_snps_vcf:         /data/refs/cnv/common_snps_preprocessed.vcf.gz
somatic_fasta:           /data/refs/cnv/main_chroms/hg38.main_chroms.fa
somatic_dict:            /data/refs/cnv/main_chroms/hg38.main_chroms.dict
```

The R annotation scripts (`exon_annotation.R`, `annotate_cnvs_complete_fixed.R`) are in `bin/`.

### Step 3 — sarek keys (FASTQ runs only)

The same params file is given to sarek. `conf/sarek.config` holds the non-path settings
(`genome = 'GATK.GRCh38'`, `wes = true`, `ascat_genome`, `save_output_as_bam = true`);
add the paths:

```yaml
intervals:        /data/refs/hg38/agilent/S07604514_Regions_Annotated.bed
ascat_alleles:    /data/refs/hg38/ascat/ascat_alleles.zip
ascat_loci:       /data/refs/hg38/ascat/ascat_loci.zip
ascat_loci_gc:    /data/refs/hg38/ascat/ascat_loci_gc.zip
ascat_loci_rt:    /data/refs/hg38/ascat/ascat_loci_rt.zip
```

### Step 4 — Runs directory

`run_full.sh` writes to `$RUNS_BASE/<run>/` (default `$HOME/runs`):

```bash
export RUNS_BASE=/your/runs/directory
```

### Step 5 — Custom images

```bash
bash docker/build_local.sh     # tags nf-sarek-clinical/{civicpy:5.4.0,r-annotation:1.0.0,sigprofiler:0.1.8}
```

---

## Inputs

### Sample naming

Normal samples end in `_N<digits>` and tumor samples in `_T<digits>`
(`<patient>_N01`, `<patient>_T01`); sarek names the pairs `<tumor>_vs_<normal>`. The
patient ID shown in the reports is the samplesheet `patient` (else the tumor name
before `_T`). One patient may have several tumors; its normal is called once.

### Samplesheet

The nf-core/sarek CSV (`patient, sex, status, sample, lane, fastq_1, fastq_2`; also
`bam`/`cram` inputs, see the sarek documentation). `sex` (`XX` / `XY`) is used by ASCAT
and gCNV; `status` 0 = normal, 1 = tumor; lanes of a sample are merged.

Optional **identification columns** (ignored by sarek, which only warns "unidentified
headers"; read by this pipeline with `--samplesheet`, which `run_full.sh` passes):

| Column | Level | Shown as |
|---|---|---|
| `patient_name`, `birth_date`, `medical_record`, `ordering_physician`, `institution`, `diagnosis` | patient (repeat on every row of the patient) | PDF cover / HTML header; `sex` from the sarek column |
| `accession_id`, `specimen_type`, `specimen_site`, `collection_date`, `received_date` | sample | tumor and normal specimen lines |
| `tumor_content` | tumor sample | pathology estimate (%), next to the ASCAT purity |

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2,patient_name,birth_date,medical_record,ordering_physician,institution,diagnosis,accession_id,specimen_type,specimen_site,collection_date,received_date,tumor_content
P001,XX,0,P001_N01,L001,/data/P001_N01_R1.fastq.gz,/data/P001_N01_R2.fastq.gz,Jane Doe,1965-03-02,MR-1234,Dr. Smith,Hospital X,Colon adenocarcinoma,A26-0002,Peripheral blood,,2026-09-01,2026-09-03,
P001,XX,1,P001_T01,L001,/data/P001_T01_R1.fastq.gz,/data/P001_T01_R2.fastq.gz,Jane Doe,1965-03-02,MR-1234,Dr. Smith,Hospital X,Colon adenocarcinoma,A26-0001,FFPE,Sigmoid colon,2026-08-20,2026-09-03,60
```

Every column is optional. A value that differs between rows of the same patient or
sample is reported as a samplesheet warning; a samplesheet without these columns keeps
the compact cover. `diagnosis` is free text; the AMP tiers use `--tumor_type` /
`--tumor_types` (CIViC disease names).

### What is read from the sarek output

Under `--sarek_outdir` (= `$RUNS_BASE/<run>/results` with `run_full.sh`):

| Path | Used for |
|---|---|
| `annotation/haplotypecaller/<normal>/*.vcf.gz` | germline SNVs, pharmacogenomics (one VCF per folder) |
| `annotation/mutect2/<pair>/*mutect2.filtered_snpEff.ann.vcf.gz` | somatic SNV tables, consensus |
| `variant_calling/strelka/<pair>/*.somatic_{snvs,indels}.vcf.gz` | consensus (optional) |
| `variant_calling/cnvkit/<pair>/*.somatic.call.cns` (else `*.call.cns`) | somatic CNV (HOMDEL / AMP) |
| `variant_calling/ascat/<pair>/*.metrics.txt`, `*segment*.txt` | purity, ploidy, WGD, LOH, GI, goodness of fit; LOH events; HRD; GATK purity |
| `variant_calling/msisensorpro/<pair>/<pair>` | MSI score |
| `reports/mosdepth/<sample>/*.recal.regions.bed.gz`, `*.mosdepth.summary.txt` | per-gene and mean coverage |
| `preprocessing/recalibrated/<sample>/*.recal.{bam,cram}` (+ index) | germline and somatic GATK CNV (BAM preferred; a CRAM needs `--cram_fasta`) |
| `pipeline_info/nf_core_sarek_software_mqc_versions.yml`, `sarek_release.txt` | tool versions and sarek release in the reports |

Pairs are those with a Mutect2 VCF or ASCAT metrics. A missing file is not an error: the
corresponding report section says "no result for this sample" (or "not analysed in this
run" when switched off with a `--skip_*` flag).

---

## Running the pipeline

### Mode A — FASTQ → reports

```bash
cd /path/to/nf-sarek-clinical
./run_full.sh --run 20260109_1229 --samplesheet /path/to/samplesheet.csv
```

```
$RUNS_BASE/20260109_1229/
├── run.config          ← work directory of sarek (written by run_full.sh)
├── work/               ← sarek work
├── results/            ← sarek output (+ pipeline_info/sarek_release.txt)
├── clinical_work/      ← nf-sarek-clinical work
└── clinical_results/   ← nf-sarek-clinical output; reports in reporting/
```

`run_full.sh` options: `--run` (required), `--samplesheet` (required unless
`--skip-sarek`; with `--skip-sarek` it still supplies the identification columns),
`--skip-sarek`, `--resume`. Other parameters go in the params file. Environment:
`RUNS_BASE`, `SAREK_REVISION` (3.5.1), `NXF_SYNTAX_PARSER` (v1).

### Mode B — sarek already run

```bash
# results in $RUNS_BASE/<run>/results
./run_full.sh --run 20260109_1229 --skip-sarek --resume

# or directly with Nextflow (any layout)
nextflow run /path/to/nf-sarek-clinical/main.nf -profile docker \
    -params-file lab_params.yml \
    -work-dir     /path/to/clinical_work \
    --sarek_outdir /path/to/sarek_results \
    --samplesheet  /path/to/samplesheet.csv \
    --run_name     20260109_1229 \
    --outdir       /path/to/clinical_results \
    [--skip_germline_cnv true ...] [-resume]
```

`nextflow run GENUIT-LBMU/nf-sarek-clinical -r main ...` pulls the repository from
GitHub instead (configure `~/.nextflow/scm` for a private repository).
`run_github.sh` does both phases like `run_full.sh` (same sarek tools and pinned
release, `--samplesheet` passed to phase 2, launched from the run directory), with
phase 2 pulled from GitHub (`CLINICAL_REPO`, `CLINICAL_REVISION` = `main`, `-latest`);
it needs `lab_params.yml` and `conf/sarek.config` next to the script.

### Several patients in one run

A samplesheet may hold any number of patients. Every pair gets its own model and
reports; germline results are keyed by the normal, somatic results by the pair, and
germline CNVs are called per normal in case mode, so a patient's results do not depend
on the other patients of the run.

### Resuming and re-running

Every report input is a task input (nothing is searched in `--outdir`), and a hash of
`bin/` and `assets/` is part of the report inputs: with `-resume` (or `--resume`), a
report is rebuilt when any of its inputs or the code changes, and taken from the cache
otherwise. Do not edit the repository (in particular `run_full.sh` and `bin/`) while a
run is going: tasks read `bin/` when they start, and bash reads `run_full.sh` as it runs.

---

## Workflow

```
sarek output ─┬─ GERMLINE_POSTPROCESS ── germline SNVs (ACMG) + pharmacogenomics
              ├─ CNV_GERMLINE_GATK ───── germline CNVs (gCNV case mode)
              ├─ SOMATIC_POSTPROCESS ─── somatic SNV tracks + CIViC / AMP, consensus,
              │                          ASCAT/CNVkit CNV, TMB, signatures, HRD
              ├─ CNV_SOMATIC_GATK ────── somatic CNVs (GATK, ASCAT purity)
              ├─ GENE_COVERAGE ───────── per-gene coverage (normal, tumor)
              └─ REPORT_MODEL (per pair: correlation + model) ─┬─ FINAL_REPORT (HTML)
                                                               └─ CLINICAL_PDF_REPORT (PDF)
```

### Germline SNVs and pharmacogenomics (`subworkflows/local/germline`)

1. `SNPEFF_MANE` — snpEff 5.2 with the MANE v1.2 RefSeq database (one NM per gene);
   MANE Plus Clinical transcripts moved after MANE Select and copied to
   `INFO/MANE_PLUS_CLINICAL` (only with `--snpeff_mane_data`).
2. `GNOMAD_ANNOTATE`, `CLINVAR_ANNOTATE`, `REVEL_ANNOTATE` — local, versioned resources
   (bcftools annotate); a variant inside the gnomAD resource regions and absent from it
   is "absent", not "unknown".
3. `SNPSIFT_FILTER_GERMLINE` — germline panel genes, HIGH / MODERATE impact (plus
   splice-region / synonymous / intronic candidates with `noncoding_splice`), DP ≥
   `germline_snv.min_dp`. Published as `germline/snv/filtered/<sample>.filtered.vcf.gz`.
4. `SPLICEAI` — SpliceAI 1.3.1 on these candidates (distance 50, masked).
5. `SNPSIFT_EXTRACT_GERMLINE` — table extraction and classification:
   field normalisation; canonical gnomAD (AF, FAF95, grpmax, source), ClinVar (class,
   stars), REVEL and SpliceAI columns; zygosity (`ZYGOSITY`, high VAF is never a
   rejection); ClinVar PS1 / PM5 evidence; filters with a `FILTER_REASON` (BA1 by the
   ACMG rule set — ClinGen VCEP gene threshold or gnomAD FAF95 — except ClinVar P/LP with
   ≥ `clinvar_min_stars` stars; ClinVar B/LB with enough stars; low VAF when DP ≥ 20;
   multi-allelic sites; non-coding candidates with SpliceAI < 0.2; nearby multi-hit
   alignment artefacts); then the ACMG classification of every candidate
   (`bin/acmg_shared.py`, see [ACMG points](#acmg-points)).
6. `PHARMACO_FILTER` — DPYD / UGT1A1 rsIDs of `--pharmaco_variants` from the
   HaplotypeCaller VCF. The report model counts the allele copies of each gene and
   derives the phenotype: DPYD activity score = 2 − 1 per no-function copy (*2A, *13)
   − 0.5 per decreased-function copy (c.2846A>T, HapB3), two variants of unknown phase
   taken as on different alleles (CPIC 2017: 2 normal; 1–1.5 intermediate, start at 50%
   and titrate; 0.5 / 0 poor, avoid); UGT1A1 decreased-function copies (*6, *28 / *37
   repeat insertions; DPWG irinotecan: ≥ 2 poor, start at 70%; 1 intermediate, no
   change). A genotyped variant outside `--panel_bed` is listed as not assessable
   (Agilent S07604514 ends ~50 bp after the UGT1A1 *28 promoter repeat). Only these six
   variants are genotyped; a recommendation must be confirmed by a clinical
   pharmacologist.

### Germline CNVs

GATK gCNV in CASE mode per normal — see
[Germline CNV: cohort mode and case mode](#germline-cnv-cohort-mode-and-case-mode).

### Somatic SNVs (`subworkflows/local/somatic`)

1. `MUTECT2_SAMPLE_ORDER` — Mutect2 VCF rewritten as normal, tumor (from its
   `##normal_sample` / `##tumor_sample` lines; the run stops if they are missing).
2. `SNPEFF_MANE` and `GNOMAD_ANNOTATE` as for germline.
3. Per gene track — **core** (60 genes, report tables), **second hit** (germline genes
   outside the core; correlation only), **TSO500** (pan-tumor genes outside the core):
   `SNPSIFT_FILTER_SOMATIC` (panel, FILTER exactly `PASS` or `haplotype`, HIGH /
   MODERATE, tumor DP ≥ `somatic_snv.min_dp`), `CIVIC_ANNOTATE` (frozen civicpy cache;
   the task fails rather than pass an unannotated VCF), `SNPSIFT_EXTRACT_SOMATIC`
   (`POP_COMMON` = gnomAD FAF95 / AF ≥ `max_pop_af`), `ANNOTATE_CONSENSUS` (`Callers`
   column: `Mutect2+Strelka2` or `Mutect2`, shown M+S / M), `CIVIC_EXTRACT`.
4. AMP tier per variant in the report model (`bin/amp_tier.py`, see the CIViC section).

### Consensus, TMB, signatures, HRD

- `SOMATIC_CONSENSUS` — Mutect2 ∩ Strelka2 (see
  [Somatic consensus calling](#somatic-consensus-calling-mutect2--strelka2)).
- `TMB_FILTER` + `TMB_CALCULATE` — consensus variants with FILTER `PASS` / `haplotype`,
  tumor DP ≥ 10, VAF ≥ 0.05, ≥ 3 ALT reads, HIGH / MODERATE impact (synonymous not
  counted), inside capture ∩ MANE Select CDS (`--panel_bed`, `--mane_cds_bed`);
  TMB = variants / territory Mb (1128 WES: 32.89 Mb). TMB-H ≥ 10 mut/Mb.
- `MUTSIG_FILTER` + `MUTSIG_SIGPROFILER` — consensus variants with ≥ 10 ALT reads;
  SigProfilerMatrixGenerator + SigProfilerAssignment (COSMIC v3.4 SBS96, DBS78, ID83),
  one run per pair. A signature counts as MMR / HRD evidence only with ≥ 10% of the
  assigned mutations, ≥ 10 mutations and ≥ 100 SBS / 50 ID in the sample (fewer:
  "not assessable").
- `HRD_PROFILER` — HRProfiler (`--hrd_organ` BREAST or OVARIAN model, the only validated
  tumor types) on the consensus VCF and the ASCAT segments; HRD+ when the probability ≥
  `biomarkers.hrd_high` (0.45, also given to HRProfiler); the `status` column records
  why a probability is missing (e.g. fewer than 5 somatic variants). Pairs without
  ASCAT segments are left out.
- MSI — the MSIsensor-pro score of sarek; MSI-H ≥ 3.5%, borderline ≥ 1.5%.

### Somatic CNVs

- `SOMATIC_CNV_INTERPRET` (`bin/parse_somatic_cnv.py`) — CNVkit calls (HOMDEL, HETDEL,
  AMP) in the 38-gene somatic CNV panel and ASCAT segments (copy-neutral LOH, HETDEL +
  LOH, HOMDEL) in the 105 germline genes as second hits; chrX / chrY LOH is not reported
  in male tumors (ASCAT sex). Also the ASCAT metrics: purity, ploidy, goodness of fit,
  WGD, LOH fraction, genomic instability.
- `CNV_SOMATIC_GATK` — CollectReadCounts (tumor, normal), DenoiseReadCounts (somatic
  PoN), CollectAllelicCounts (common SNPs), ModelSegments (paired), CallCopyRatioSegments,
  plots, and `annotate_cnvs_complete_fixed.R`: MANE exons of `--cnv_gatk_genes`,
  absolute copy number corrected by the ASCAT purity of the tumor (`PURITY_SOURCE`;
  the ASCAT ploidy is not used). A failure of the R script fails the task; a tumor
  without events gets a header-only table.
- The report model harmonises both (`bin/cnv_shared.py`): HOMDEL and AMP (CN ≥ 5) in
  the somatic CNV panel, LOH in the germline panel, with the tools that support each
  event and a confidence (HIGH: ASCAT / CNVkit and GATK both report it).

### Coverage

`GENE_COVERAGE` — mean / min / max depth per gene from the mosdepth regions of each
sample (germline panel genes). A panel gene the capture BED names otherwise (RAD51D
appears as RAD51L3-RFFL) is found by its MANE exon coordinates (`--gcnv_mane_txt`).

### Report model and reports

`REPORT_MODEL` (one task per pair) runs `germline_somatic_correlate.py` (two-hit events
in MMR / HRR and other tumor-suppressor genes: germline P/LP + somatic SNV / LOH /
deletion; biallelic somatic events; pathway biomarkers — MMR: MSI-H, or TMB-H with an
MMR signature; HRR: HRD+ or SBS3 / SBS8) and `bin/report_model.py`, which reads every
input once and writes `<pair>.report_model.json`:

- identification (samplesheet), IDs, date, versions and resources, thresholds in force;
- metrics and statuses (TMB-H, MSI status, HRD status, "Indeterminate" when tumor
  coverage < 50X), ASCAT alerts (purity < 20%, goodness of fit < 90%), ploidy / WGD note,
  TERT promoter note;
- germline variants with ACMG class and criteria; germline CNV segments (reportable,
  low-mappability);
- somatic variants with tier, tier-defining CIViC record, evidence summary, gene context
  and consequence; Tier IV; TSO500 Tier III; second-hit variants; actionable genes;
- somatic CNV (HOMDEL, AMP, LOH; the LOH the PDF lists and the ones it leaves out);
- signatures (assigned, percentages, COSMIC aetiology, MMR / HRD calls);
- correlation, CPIC matches, per-gene coverage, raw TMB / HRD tables;
- clinical trials (queried once; date, queries, trials or the connection error);
- analyses not run / without result.

`FINAL_REPORT` (`bin/render_html.py`) and `CLINICAL_PDF_REPORT`
(`bin/clinical_pdf_report.py`, ReportLab) render it.

**Clinical trials.** Recruiting studies on ClinicalTrials.gov (API v2), up to 3 per
query: genes with an AMP Tier I/II variant, genes of the germline–somatic findings, and
the biomarkers (MSI-H, TMB-H, HRD+), each with the patient's tumor type as condition;
MSI-H and TMB-H also in any solid tumor. Trials restricted to the opposite MSI status
are dropped. The report prints the query date.

---

## Germline CNV: cohort mode and case mode

GATK gCNV (`DetermineGermlineContigPloidy` + `GermlineCNVCaller`) works in two modes:

| | **COHORT mode** | **CASE mode** |
|---|---|---|
| What it does | Learns, from a batch of samples, the contig-ploidy model and the read-depth denoising / calling model (per-interval bias, GC and capture effects, copy-number priors), and calls CNVs in those same samples | Calls CNVs in new samples against a **frozen** cohort model |
| Input | ≥ 30 normals (more is better) of the same capture kit and protocol | One sample (or a few) + the cohort model + the exact interval list the model was trained on |
| Output | `ploidy-model/`, `<prefix>-model/` (the "PoN"), and calls for the batch | Calls for the new sample(s) |
| Result depends on | The whole batch (adding or removing samples changes every call) | The sample and the versioned model only |
| Used by this pipeline | No — run once, offline, to build the model | **Yes**, for every normal sample |

### What the pipeline runs (CASE mode, per normal sample)

`CNV_GERMLINE_GATK` (`subworkflows/local/cnv_germline_gatk`) on the recalibrated
normal BAM (or CRAM, decoded with `--cram_fasta`) from sarek:

| Step | Process | GATK tool / script | Key settings |
|---|---|---|---|
| 1 | `GATK_COLLECT_READ_COUNTS_GERMLINE` | `CollectReadCounts` | `-L --gcnv_interval_list`, `--interval-merging-rule OVERLAPPING_ONLY`, HDF5 |
| 2 | `GATK_DETERMINE_GERMLINE_CONTIG_PLOIDY` | `DetermineGermlineContigPloidy` (case) | `--model --gcnv_ploidy_model` |
| 3 | `GATK_GERMLINE_CNV_CALLER` | `GermlineCNVCaller --run-mode CASE` | `--model --gcnv_pon_dir`, contig-ploidy calls of step 2 |
| 4 | `GATK_POSTPROCESS_GERMLINE_CNV_CALLS` | `PostprocessGermlineCNVCalls` | allosomal contigs chrX / chrY; denoised copy ratios, genotyped intervals and segments VCF |
| 5 | `MANE_ANNOTATE_GERMLINE_CNV` | `bin/exon_annotation.R` | MANE exons (`--gcnv_mane_txt`), germline panel genes minus `filters.germline_cnv.exclude_genes`, GIAB low-mappability fraction per segment (`--gcnv_lowmap_bed`) |

The report model then keeps the deletions / duplications (gCNV ALT `<DEL>` / `<DUP>`,
against the contig ploidy, so one chrX copy in a male is reference) in the germline
panel genes with QA ≥ `filters.germline_cnv.min_qa` (4), and lists apart, as not
reportable without orthogonal confirmation (MLPA), the segments with ≥
`filters.germline_cnv.lowmap_min_frac` (50 %) of their exon bases in low-mappability
regions (CHEK2 ex11–15, PMS2 / PMS2CL, SDHA, FANCD2, NF1 …).

Consequences of case mode:

- **Patient-independent results.** A sample is called against the same model whatever
  else is in the run; re-running a patient alone or in a batch gives the same calls.
- **Any batch size.** One patient or many per run (each tumor/normal pair gets its own
  reports; each normal is called once even if it pairs with several tumors).
- **Traceable.** The reports list the ploidy model, the cohort model and the interval
  list used (Versions and Resources table).
- **The model must match the data.** Same capture kit, same library protocol, same
  reference and aligner, same GATK / gcnvkernel version
  (`broadinstitute/gatk:4.6.2.0`, gcnvkernel 0.9; a 0.8 model cannot be read), and
  `--gcnv_interval_list` must be the filtered interval list the model was trained on.

Required parameters: `--gcnv_interval_list`, `--gcnv_ploidy_model`, `--gcnv_pon_dir`,
`--gcnv_fasta` (+ `.fai`), `--gcnv_dict`, `--gcnv_mane_txt`, `--gcnv_lowmap_bed`
(and `--cram_fasta` for CRAM input). `--skip_germline_cnv true` switches it off; the
reports then say "not analysed in this run".

### Building the cohort model (COHORT mode, offline)

Run once per capture kit / protocol / GATK version, with the normals of your own
laboratory (≥ 30; 50–100 or more gives a more stable model; both sexes; same kit and
protocol as the patients; drop samples with failed QC, contamination or known large
CNVs). The standard GATK sequence is below (recommended settings; the workflow used
to build the current model is in
[FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)):

```bash
GATK="docker run --rm -v $PWD:$PWD -w $PWD broadinstitute/gatk:4.6.2.0 gatk"

# 1. Capture intervals, no binning, 250 bp padding
$GATK PreprocessIntervals -R hg38.fa -L capture.bed --bin-length 0 --padding 250 \
     --interval-merging-rule OVERLAPPING_ONLY -O preprocessed.interval_list
# 2. GC content per interval
$GATK AnnotateIntervals -R hg38.fa -L preprocessed.interval_list \
     --interval-merging-rule OVERLAPPING_ONLY -O annotated_intervals.tsv
# 3. Read counts of every normal (recalibrated BAM / CRAM from sarek)
for s in normals/*.bam; do
  $GATK CollectReadCounts -I $s -L preprocessed.interval_list -R hg38.fa \
       --interval-merging-rule OVERLAPPING_ONLY --format HDF5 -O counts/$(basename $s .bam).counts.hdf5
done
# 4. Drop intervals with extreme GC, low mappability or low / variable coverage in the cohort
$GATK FilterIntervals -L preprocessed.interval_list --annotated-intervals annotated_intervals.tsv \
     $(printf -- '-I %s ' counts/*.hdf5) --interval-merging-rule OVERLAPPING_ONLY \
     -O exomas.interval_list                                   # → --gcnv_interval_list
# 5. Contig-ploidy model (COHORT)
$GATK DetermineGermlineContigPloidy -L exomas.interval_list $(printf -- '-I %s ' counts/*.hdf5) \
     --contig-ploidy-priors contig_ploidy_prior.tsv --interval-merging-rule OVERLAPPING_ONLY \
     --output . --output-prefix ploidy                          # ploidy-model/ → --gcnv_ploidy_model
# 6. Denoising + calling model (COHORT), one shard
$GATK GermlineCNVCaller --run-mode COHORT -L exomas.interval_list $(printf -- '-I %s ' counts/*.hdf5) \
     --contig-ploidy-calls ploidy-calls --annotated-intervals annotated_intervals.tsv \
     --interval-merging-rule OVERLAPPING_ONLY --output . --output-prefix PoN_CNV
                                                                # PoN_CNV-model/ → --gcnv_pon_dir
```

`contig_ploidy_prior.tsv` gives, per contig, the prior of ploidy 0–4 (autosomes 2 with
0.996; chrX 1 or 2 and chrY 0 or 1 for a mixed-sex cohort). The pipeline expects a
single model directory (one shard), so the cohort model is built without scattering.

The model in use (`gcnvkernel_version.json` 0.9) was built on ~221,600 exome intervals
with: max copy number 5, `p_alt` 5 × 10⁻⁴, `p_active` 0.1, CNV / class coherence length
10 kb, up to 6 bias factors with explicit GC-bias modelling (20 GC bins); ploidy model
`mapping_error_rate` 0.3, `psi_j_scale` 10⁻³. These values are stored in
`calling_config.json`, `denoising_config.json` and `ploidy_config.json` inside the
model directories.

Rebuild the model when the capture kit, the library protocol, the reference, the
aligner or the GATK version changes, and re-validate it on samples with known CNVs
(MLPA-confirmed deletions / duplications) before clinical use.

### Calling a batch directly in COHORT mode

Steps 1–6 above also return calls for every sample of the batch (`PoN_CNV-calls/`),
and `PostprocessGermlineCNVCalls` turns them into segment VCFs. This is how a research
batch of ≥ 30 samples can be called in cohort mode. The pipeline does not do it: for
clinical reports it calls each patient in CASE mode against the fixed, versioned
model, so that a patient's result does not depend on the other samples of the batch
and can be reproduced later. To include a patient's own normal in a new model, rebuild
the model (above) and point `--gcnv_ploidy_model` / `--gcnv_pon_dir` /
`--gcnv_interval_list` to the new directories.

### Somatic CNV (GATK): the same principle

`CNV_SOMATIC_GATK` calls each tumor (with its normal for allelic counts) against a
read-count Panel of Normals built offline with `CreateReadCountPanelOfNormals` from the
same normals (`--somatic_pon_hdf5`, `--somatic_interval_list`, `--common_snps_vcf`).
ASCAT and CNVkit (run by sarek) need no panel of normals: they compare each tumor with
its matched normal.

---

## Reports

Both reports render `<pair>.report_model.json`; a section of an analysis switched off
with `--skip_*` says "not analysed in this run (--skip_x true)", and an analysis without
a result for the pair says "no result for this sample", never an empty negative result.

| Clinical PDF (`<pair>_clinical_report.pdf`) | Laboratory HTML (`<pair>_molecular_report.html`) |
|---|---|
| **Page 1:** identification (patient and specimens from the samplesheet, run, date, pipeline and sarek versions), global biomarkers (TMB, MSI, HRD probability and status), secondary biomarkers (WGD, GI, LOH, ploidy with "(WGD)"), quality (purity, normal / tumor coverage), ASCAT fit and low-coverage alerts, **major clinical findings** (MMR deficiency, HRD, CIViC actionable variants with their evidence, high-confidence CNVs, germline P/LP), disclaimer | Header with identification and versions; dashboard (TMB, MSI, LOH, GI, WGD, ploidy, purity, coverage N / T, HRD probability and status, "Indeterminate" below 50X tumor coverage); ASCAT alerts |
| Germline variants P / LP / VUS (ACMG class, criteria, REVEL, SpliceAI, gnomAD AF / FAF, ClinVar with stars, zygosity, QC flags) | Germline SNVs: every candidate (P, LP, VUS, LB, B) with filter buttons |
| Germline CNVs (reportable segments; low-mappability segments apart) | Somatic SNVs: every tier with filters, expandable gene context and CIViC evidence per variant; Tier IV and TSO500 Tier III tables; TERT note |
| Somatic variants: actionable (Tier I–II), other (Tier III), TSO500 Tier III, Tier IV population polymorphisms; gene annotation notes | Raw TMB and HRD tables |
| Mutational signatures (SBS and ID, counts per type, COSMIC aetiology, spectra) | Pharmacogenomics (CPIC) |
| Somatic CNVs: HOMDEL, AMP, LOH limited to genes with a germline P/LP or somatic variant, with the number and names of the LOH left out; ploidy / WGD note | Per-gene coverage (normal and tumor, mean and minimum) |
| Germline–somatic correlation and somatic biallelic events | Signatures SBS96, DBS, ID with plots |
| Pharmacogenomics (CPIC) | Germline–somatic correlation |
| Clinical trials | Germline CNVs; somatic CNVs with every LOH (the ones the PDF lists marked) |
| Per-gene coverage grids (germline panel in the normal, somatic panel in the tumor) | Clinical trials |
| Methodology (tools with the versions recorded by sarek, thresholds, references) and Versions and Resources table | Methodology and Versions and Resources table |

The PDF is in English or Spanish (`--report_language en|es`); the HTML is in English.
The versions block lists the pipeline version and git commit (marked "+ local changes"
when the working tree differs), the sarek release, Nextflow, the tool versions sarek
recorded (BWA-MEM, fastp, GATK, Strelka2, ASCAT, CNVkit, MSIsensor-pro, mosdepth, snpEff)
and the resource files used (gnomAD, ClinVar, ClinVar protein index, REVEL, CIViC
cache, snpEff MANE database, MANE exons and CDS, capture BED, gCNV models and intervals,
somatic GATK PoN, COSMIC v3.4).

---

## Reference files

All Docker images contain the bioinformatics tools. The files below must be downloaded or
built once and stored locally. Configure their paths in your params file.

| Resource | Parameter | Needed for |
|---|---|---|
| SnpSift JAR | `snpsift_jar` | always |
| gnomAD v4.1 panel | `gnomad_tsv` | germline SNVs (required); somatic Tier IV |
| ClinVar panel | `clinvar_tsv` | germline SNVs (required) |
| REVEL v1.3 panel | `revel_tsv` | germline SNVs (required) |
| ClinVar protein index | `clinvar_protein_index` | PS1 / PM5 (recommended) |
| snpEff MANE database | `snpeff_mane_data` | MANE Select nomenclature (recommended) |
| CIViC cache | `civic_cache` | somatic SNVs (required) |
| MANE exons / CDS | `gcnv_mane_txt`, `mane_cds_bed` | germline CNV, GATK somatic CNV, coverage fallback; TMB territory |
| Capture BED | `panel_bed`, `gene_coords_bed` | TMB, TERT note; somatic CNV interpretation |
| hg38 FASTA + dict | `gcnv_fasta`, `gcnv_dict` (`somatic_fasta`, `somatic_dict`) | consensus, SpliceAI, GATK CNV |
| Low-mappability BED | `gcnv_lowmap_bed` | germline CNV (required) |
| gCNV cohort model | `gcnv_ploidy_model`, `gcnv_pon_dir`, `gcnv_interval_list` | germline CNV |
| Somatic GATK PoN | `somatic_pon_hdf5`, `somatic_interval_list`, `common_snps_vcf` | GATK somatic CNV |
| ASCAT resources | `ascat_*` (sarek) | FASTQ runs |

### SnpSift — VCF filtering and field extraction

Part of the **snpEff** package:

```bash
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip -d /data/refs/snpEff
# Jar → /data/refs/snpEff/snpEff/SnpSift.jar
```

### gnomAD v4.1 — population frequencies (required for germline SNVs)

A small local subset of gnomAD v4.1 *joint* (exomes + genomes, ~807k individuals)
restricted to the MANE exons (±100 bp) of the panel genes, and to the whole gene
body of the germline genes (`--gene-body`: intronic candidates kept for SpliceAI
need a frequency for BA1). No patient coordinates leave the host, every variant gets the same
gnomAD version, and a variant inside the resource regions that is not in
gnomAD is reported as truly absent (`gnomAD_AF_source = gnomAD_v4.1_absent`).
Sites failing gnomAD QC in both exomes and genomes (`*_qc_fail`, e.g. paralogous
regions such as FANCD2/FANCD2P2) keep their AF for BA1/BS1 but never support PM2.

```bash
# needs bcftools built with libcurl — the pipeline's bcftools container has it
docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -v /data/refs:/data/refs -w $PWD \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash scripts/build_gnomad_panel.sh \
        --genes <(awk -F'\t' 'NR>1 && $1!="RYBP" && $1!="TERC" {print $1}' assets/gene_panels/genes.tsv) \
        --gene-body assets/gene_panels/germline_cancer_predisposition.txt \
        --mane  /data/refs/mane/exons_mane_v1.2.txt \
        --outdir /data/refs/gnomad
# → gnomad_v4.1_panel.tsv.gz(.tbi), gnomad_v4.1_panel_regions.bed.gz(.tbi), manifest
```

All 541 genes of `genes.tsv` with a MANE v1.2 transcript (germline, somatic core,
TSO500, somatic CNV, MMR / HRR), plus the gene body of the 105 germline genes (~8 Mb):
~180 MB, about 6 h (remote reads of
the gnomAD bucket). The exon table must include the MANE Plus Clinical exons
(`scripts/build_mane_exons.py` does), or variants on BRAF NM_001374258.1, GNAS XLαs,
CDKN2A p14ARF, PTCH1, MITF, … get no frequency. Then set
`gnomad_tsv: /data/refs/gnomad/gnomad_v4.1_panel.tsv.gz`. Rebuild it whenever
`genes.tsv` changes. Germline TSVs gain `gnomAD_AF`, `gnomAD_AF_source`
and the `GNOMAD_*` columns (FAF95 grpmax, ancestry group, admixed-American AF,
homozygote counts, gnomAD FILTER, and grpmax AF / AC / AN for the ClinGen VCEP
rules; resources built before 2026-10-23 lack these and use the general rules).

### CIViC cache — somatic clinical evidence (required for somatic postprocessing)

civicpy reads CIViC from a pickled cache. The pipeline never downloads it at run
time: a dated snapshot is frozen once and passed with `--civic_cache`, so the same
cache always gives the same annotation (civicpy would otherwise replace a cache
older than 7 days with the current nightly). The pickle is tied to the civicpy
version, so build it with the pipeline's civicpy image:

```bash
docker run --rm -u $(id -u):$(id -g) -v $PWD:$PWD -v /data/refs:/data/refs -w $PWD \
    nf-sarek-clinical/civicpy:5.4.0 \
    bash scripts/build_civic_cache.sh --outdir /data/refs/civic
# → civicpy_cache_<YYYYMMDD>.pkl + manifest (snapshot time, civicpy version, counts, md5)
```

Then set `civic_cache: /data/refs/civic/civicpy_cache_<YYYYMMDD>.pkl`. Evidence with
status `accepted` and `submitted` (not yet reviewed by CIViC curators) is annotated;
`rejected` is not. CIViC data are CC0.

**AMP tier** (`bin/amp_tier.py`, the same for HTML, PDF and the findings card).
Only accepted evidence whose direction is `SUPPORTS` sets the tier; submitted
evidence and evidence against are listed but never raise it.

| CIViC level | Patient's tumor type (or pan-cancer) | Other tumor type |
|---|---|---|
| A | IA | IIC |
| B | IB | IIC |
| C | IIC | — |
| D, E | IID | — |

The tumor type is given as CIViC disease names joined by `|` (a name also covers
its subtypes: `Breast Cancer` matches `HER2 Positive Breast Cancer`), for the whole run
(`--tumor_type "Colorectal Cancer|Colon Cancer|Rectum Cancer|Colorectal Adenocarcinoma"`)
or per patient with `--tumor_types`, a TSV `id<TAB>terms` where `id` is the
patient, the tumor sample or the pair. Without either, every disease counts as
"other tumor type" and the highest tier is IIC (the run prints a warning).

**Tier IV** (AMP/ASCO/CAP 2017: benign / likely benign, MAF ≥ 1% in the general
population). With `--gnomad_tsv`, the somatic Mutect2 VCF is annotated with the local
gnomAD v4.1 resource (`GNOMAD_ANNOTATE`; the resource covers all 541 genes of `genes.tsv`, so the core, second-hit and TSO500 tracks); a
variant with gnomAD FAF95 (else AF) ≥ `filters.somatic_snv.max_pop_af` (0.01) is Tier IV
whatever its CIViC evidence. Tier IV variants are listed in a separate table as probable
germline polymorphisms and are not counted as somatic findings, double hits or second
hits (e.g. 1128 BRCA2 p.Thr1915Met, gnomAD AF 2.7%, called at VAF 0.125). Without
`--gnomad_tsv` there is no population frequency and no Tier IV (the run warns).

### Low-mappability BED — germline CNV flags (required for germline CNV)

GIAB genome stratifications v3.3 (`GRCh38_lowmappabilityall`) restricted to the
MANE exons (±1 kb) of the germline panel. Read counts in pseudogene / paralog
regions (CHEK2 exons 11–15, PMS2 / PMS2CL, SDHA, FANCD2, NF1, …) are unreliable:
reads map with low MAPQ and gCNV calls recurrent artefacts (CHEK2 exon 14 HOMDEL
in most normals). Segments with ≥ `filters.germline_cnv.lowmap_min_frac` (0.5) of
their panel exon bases in these regions are listed in a separate table that asks
for orthogonal confirmation (MLPA), and are not used by the germline–somatic
correlation. They are not hidden: real PMS2 exon 12–15 deletions fall there.

```bash
bash scripts/build_lowmap_panel.sh \
    --genes assets/gene_panels/germline_cancer_predisposition.txt \
    --mane  /data/refs/mane/exons_mane_v1.2.txt \
    --outdir /data/refs/cnv/lowmap
# → GRCh38_lowmappabilityall_panel.bed.gz, manifest (needs curl + bedtools)
```

Then set `gcnv_lowmap_bed: /data/refs/cnv/lowmap/GRCh38_lowmappabilityall_panel.bed.gz`.
The annotated germline CNV TSV gains `LOWMAP` (fraction of each exon) and
`LOWMAP_SEG` (fraction of the segment's panel exon bases).

### snpEff MANE database — germline and somatic nomenclature (recommended)

Sarek runs snpEff with `-canon`, which keeps the transcript with the longest CDS.
For ~1 in 10 panel variants that is not the MANE Select transcript, so the
reported HGVS differs from ClinVar and the literature (e.g. germline BRCA1
p.Ser1634Gly instead of p.Ser1613Gly; somatic BRAF p.Lys641Glu instead of
p.Lys601Glu, CHEK2 p.Arg389Cys instead of p.Arg346Cys). With `snpeff_mane_data`
set, the germline VCF and the somatic Mutect2 VCF are re-annotated (`SNPEFF_MANE`,
snpEff 5.2) with the MANE database, one RefSeq NM transcript per gene. The somatic
consensus, TMB and signatures keep the sarek annotation.

The database also holds the 64 MANE Plus Clinical transcripts, which snpEff may
list first (BRAF, FGFR2, KRAS, MUTYH and SMARCA4 among the panel genes).
`SNPEFF_MANE` moves them after the MANE Select transcript
(`--mane_plus_clinical`, default `assets/mane/mane_v1.2_plus_clinical.txt`, from
the NCBI MANE v1.2 summary) and copies their HGVS to `INFO/MANE_PLUS_CLINICAL`.
The reports show the transcript under HGVS_c and, when it differs, the MANE Plus
Clinical HGVS under HGVS_p (MUTYH p.Gly368Asp with `NM_001128425.2: p.Gly396Asp`).

```bash
mkdir -p /data/refs/snpeff && cd /data/refs/snpeff
curl -O https://snpeff-public.s3.amazonaws.com/databases/v5_2/snpEff_v5_2_GRCh38.mane.1.2.refseq.zip
unzip snpEff_v5_2_GRCh38.mane.1.2.refseq.zip     # → data/GRCh38.mane.1.2.refseq/ (~370 MB)
```

Then set `snpeff_mane_data: /data/refs/snpeff/data` (`snpeff_mane_db` defaults to
`GRCh38.mane.1.2.refseq`). MANE 1.2 transcripts are identical (same NM version) to
MANE 1.5 for all 122 genes of the three panels.

### ClinVar — germline classifications (required for germline SNVs)

A local subset of the official NCBI ClinVar GRCh38 VCF over the same panel
regions. Matching is by allele, so indels are covered, and the release date is recorded
in every output. When `clinvar_tsv` is set it is the **only** ClinVar source: a
variant absent from that release is reported as not in ClinVar.

```bash
# needs the gnomAD regions BED above and an hg38 FASTA (+ .fai)
docker run --rm -u $(id -u):$(id -g) -v /data:/data -w $PWD \
    quay.io/biocontainers/bcftools:1.20--h8b25389_0 \
    bash scripts/build_clinvar_panel.sh \
        --clinvar /data/refs/clinvar/clinvar_20260822.vcf.gz \
        --regions /data/refs/gnomad/gnomad_v4.1_panel_regions.bed.gz \
        --fasta   /data/refs/hg38/hg38.fa \
        --outdir  /data/refs/clinvar
# → clinvar_20260822_panel.tsv.gz(.tbi), manifest
```

Then set `clinvar_tsv: /data/refs/clinvar/clinvar_20260822_panel.tsv.gz`.
ClinVar is released monthly; rebuild to update. Germline TSVs gain
`ClinVar_sig`, `ClinVar_stars` (0–4 from the review status), `ClinVar_revstat`,
`ClinVar_id`, `ClinVar_conf`, `ClinVar_dn` and `ClinVar_source`; the reports show
the classification with its stars, e.g. "Benign/Likely benign (2*)".

The ClinVar class is final only with at least `filters.acmg.clinvar_min_stars`
review stars (default 2: multiple submitters without conflicts, expert panel or
practice guideline):

| ClinVar | >= min stars | below min stars |
|---|---|---|
| P / LP | final | ACMG score, never below VUS (except BA1); criteria carry `ClinVar:LP(1*)` for manual review |
| VUS | VUS, but BS1 / BA1 can bring it down to LB / B | ACMG score |
| B / LB | final; removed from the candidate table (`clinvar_benign`) | kept; ACMG score |
| Conflicting / other | ACMG score | ACMG score |

ClinVar P/LP still exempts a variant from the VAF floor and the multi-hit filter
regardless of stars.

### ClinVar protein index — ACMG PS1 / PM5 (recommended)

ClinVar P/LP missense variants in MANE Select protein coordinates. For a germline
missense variant, `ClinVar_PS1` lists P/LP missense with the same amino-acid change
(other nucleotide) and `ClinVar_PM5` those with another change at the same residue;
references need >= `clinvar_min_stars` stars and are not at exon edges
(`splice_region`). See the header of `scripts/build_clinvar_protein_index.py` for
the bcftools + snpEff MANE commands, then:

```bash
python3 scripts/build_clinvar_protein_index.py --vcf clinvar_panel.mane.vcf \
    --mane-plus-clinical assets/mane/mane_v1.2_plus_clinical.txt --outdir /data/refs/clinvar
# → clinvar_20260822_protein_index.tsv.gz, manifest
```

Set `clinvar_protein_index: /data/refs/clinvar/clinvar_20260822_protein_index.tsv.gz`
and rebuild it with each ClinVar release.

### REVEL v1.3 — missense scores per MANE transcript (required for germline SNVs)

A panel subset of REVEL v1.3 (Zenodo 7072866, ODbL). REVEL gives one score per
amino-acid change; the build keeps, for every allele, the score of the MANE
Select transcript (`REVEL_source = REVEL_v1.3_MANE`), or the highest score over
the other transcripts when MANE has none (`REVEL_v1.3_other_tx`).

```bash
# download the 25 zips of https://zenodo.org/records/7072866 into /data/refs/revel/source
python3 scripts/build_revel_panel.py \
    --revel-dir /data/refs/revel/source \
    --regions   /data/refs/gnomad/gnomad_v4.1_panel_regions.bed.gz \
    --mane      /data/refs/mane/MANE.GRCh38.v1.5.summary.txt.gz \
    --outdir    /data/refs/revel
# → revel_v1.3_panel.tsv.gz(.tbi), manifest   (needs python3, bgzip, tabix)
```

Then set `revel_tsv: /data/refs/revel/revel_v1.3_panel.tsv.gz`.

### SpliceAI — splicing predictions

With `spliceai: true` (default) and a GRCh38 FASTA (`spliceai_fasta`, defaults to
`gcnv_fasta`), SpliceAI 1.3.1 runs on the germline candidates that passed the
panel + impact filter (one to two hundred variants per sample, minutes on CPU;
distance 50, masked scores). `SpliceAI_DS` (max delta) gives PP3_Supporting when
>= 0.2 and blocks a REVEL BP4 when > 0.1 (ClinGen, Walker et al. 2023).

**Non-coding splice candidates** (`filters.germline_snv.noncoding_splice = true`,
needs SpliceAI): besides HIGH / MODERATE, the panel filter lets through LOW /
MODIFIER variants of panel genes whose effect is `splice_region_variant`,
`splice_donor_5th_base_variant`, `splice_donor_region_variant`,
`splice_polypyrimidine_tract_variant`, `synonymous_variant` or `intron_variant`
(UTR, upstream and downstream stay out). Step 2 keeps them only when
`SpliceAI_DS >= acmg.spliceai_pp3` (0.2) or ClinVar is P/LP (any stars); the rest
go to the filtered table as `spliceai_below_threshold (non-coding, DS=...)`. Kept
variants are scored like any other (typically PP3_Sup(splice) + PM2_Sup = VUS;
ClinVar >= 2 stars is final). On ClinVar >= 2 stars non-coding variants of the
panel genes, SpliceAI >= 0.2 recovers 471 / 503 P/LP and 18 / 1,478 B/LB (random
sample); P/LP below 0.2 (e.g. PALB2 c.3350+5G>A, RB1 c.264+5G>A) are kept through
the ClinVar exemption.

### hg38 reference genome

From the **GATK GCP bucket**:

```bash
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     /data/refs/hg38/hg38.fa
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai /data/refs/hg38/hg38.fa.fai
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict      /data/refs/hg38/hg38.dict

# Create main-chromosomes FASTA (required for somatic GATK CNV):
samtools faidx /data/refs/hg38/hg38.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 \
    chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
    > /data/refs/cnv/main_chroms/hg38.main_chroms.fa
samtools faidx /data/refs/cnv/main_chroms/hg38.main_chroms.fa
gatk CreateSequenceDictionary -R /data/refs/cnv/main_chroms/hg38.main_chroms.fa
```

### WES capture BED

Download from your kit vendor (`--panel_bed`). Used for the TMB territory (capture ∩ MANE Select CDS) and to tell whether the TERT promoter hotspots are covered.

- **Agilent SureSelect**: [earray.chem.agilent.com/suredesign](https://earray.chem.agilent.com/suredesign/) (free registration)
- **Twist / KAPA / Roche**: request from the manufacturer

### MANE exon coordinates and CDS

From the NCBI MANE v1.2 RefSeq GTF (the release of the snpEff MANE database):

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.refseq_genomic.gtf.gz
python3 scripts/build_mane_exons.py --gtf MANE.GRCh38.v1.2.refseq_genomic.gtf.gz \
    --out /data/refs/mane/exons_mane_v1.2.txt --cds-bed /data/refs/mane/cds_mane_v1.2.bed
```

`exons_mane_v1.2.txt` (chr, start, end, gene_symbol, exon, transcript_id; MANE Select
and MANE Plus Clinical exons) → `gcnv_mane_txt`; `cds_mane_v1.2.bed` (merged MANE Select
CDS) → `mane_cds_bed`.

### GATK CNV support files and Panel of Normals

The complete workflow (interval preprocessing, common SNPs VCF, PoN construction) is provided at:
**[github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)**

```bash
git clone https://github.com/FlorPio/gatk-cnv-somatic
# Follow that repo's README to generate:
#   exomas.interval_list          ← germline FilterIntervals output
#   preprocessed.interval_list   ← somatic PreprocessIntervals output
#   common_snps_preprocessed.vcf.gz + .tbi
#   ploidy-model/                 ← DetermineGermlineContigPloidy (COHORT mode)
#   PoN_CNV/                      ← GermlineCNVCaller (COHORT mode, ≥ 30 normals)
#   wes-do-gc.pon.hdf5            ← CreateReadCountPanelOfNormals
```

> PoN files are **institution-specific**: built from your own cohort of ≥ 30 normal samples with the same kit and the same GATK version (`broadinstitute/gatk:4.6.2.0`). Rebuild if GATK version changes. See [Germline CNV: cohort mode and case mode](#germline-cnv-cohort-mode-and-case-mode) for the germline model.

> The R annotation scripts (`exon_annotation.R`, `annotate_cnvs_complete_fixed.R`) from that repo are already **bundled in `bin/`** — no manual copy needed.

### ASCAT resources (Mode A only)

Required for ASCAT allele-specific copy-number in nf-core/sarek. Download from the **ASCAT GitHub**:
[github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES)

Download the four zip files for hg38 (`alleles`, `loci`, `loci_gc`, `loci_rt`) and set their paths in your params file (`ascat_alleles`, `ascat_loci`, `ascat_loci_gc`, `ascat_loci_rt`; see First-time setup, Step 3).

### Suggested directory layout

```
/data/refs/
├── snpEff/
│   └── snpEff/SnpSift.jar
├── hg38/
│   ├── hg38.fa  (+ .fai, .dict)
│   ├── ascat/   (alleles.zip, loci.zip, loci_gc.zip, loci_rt.zip)
│   └── <kit>/
│       ├── <kit>_Regions.bed
│       ├── <kit>_Regions_Annotated.bed   ← for sarek intervals
│       └── <kit>_Padded.withGenes.bed
└── cnv/
    ├── (exons_mane_v1.2.txt: see mane/, scripts/build_mane_exons.py)
    ├── exomas.interval_list
    ├── preprocessed.interval_list
    ├── common_snps_preprocessed.vcf.gz  (+ .tbi)
    ├── ploidy-model/
    ├── PoN_CNV/
    ├── wes-do-gc.pon.hdf5
    └── main_chroms/
        ├── hg38.main_chroms.fa  (+ .fai)
        └── hg38.main_chroms.dict
```

---

## Parameters reference

All thresholds are in [`conf/filters.config`](#filtering-thresholds-conffiltersconfig);
the parameters below select inputs, resources and modules.

### Run

| Parameter | Default | Description |
|---|---|---|
| `--sarek_outdir` | — | nf-core/sarek results directory (required) |
| `--outdir` | `./results` | Output directory |
| `--run_name` | — | Run identifier, shown in the reports |
| `--samplesheet` | — | sarek samplesheet; its optional identification columns go to the reports |
| `--sarek_version` | `<sarek_outdir>/pipeline_info/sarek_release.txt` | sarek release shown in the reports |
| `--genome_build` | `GRCh38` | Genome build (CIViC, SigProfiler, HRProfiler) |
| `--report_language` | `en` | PDF language: `en` or `es` |
| `--publish_dir_mode` | `copy` | Nextflow publishDir mode |
| `--max_cpus` / `--max_memory` / `--max_time` | `16` / `128.GB` / `240.h` | Cap of every task (`process.resourceLimits`) |
| `--docker_namespace` | `nf-sarek-clinical` | Prefix of the custom images (civicpy, r-annotation, sigprofiler) |

### Annotation resources

| Parameter | Default | Description |
|---|---|---|
| `--snpsift_jar` | — | `SnpSift.jar` (required) |
| `--gnomad_tsv` | — | Local gnomAD v4.1 panel (`scripts/build_gnomad_panel.sh`); path or https URL. **Required** for germline SNVs; somatic Tier IV |
| `--gnomad_regions` | `<gnomad_tsv>_regions.bed.gz` | Regions BED of that resource (only if renamed) |
| `--clinvar_tsv` | — | Local ClinVar panel (`scripts/build_clinvar_panel.sh`). **Required** for germline SNVs |
| `--revel_tsv` | — | Local REVEL v1.3 panel (`scripts/build_revel_panel.py`). **Required** for germline SNVs |
| `--clinvar_protein_index` | — | ClinVar P/LP missense in MANE protein coordinates, PS1 / PM5 (`scripts/build_clinvar_protein_index.py`) |
| `--snpeff_mane_data` | — | snpEff data dir with the MANE database (germline and somatic MANE Select nomenclature) |
| `--snpeff_mane_db` | `GRCh38.mane.1.2.refseq` | snpEff MANE database |
| `--mane_plus_clinical` | `assets/mane/mane_v1.2_plus_clinical.txt` | MANE Plus Clinical transcripts of that release |
| `--spliceai` | `true` | SpliceAI on the germline candidates |
| `--spliceai_fasta` | `--gcnv_fasta` | GRCh38 FASTA (+ `.fai`) for SpliceAI |
| `--civic_cache` | — | Frozen civicpy cache (`scripts/build_civic_cache.sh`). **Required** for somatic SNVs |
| `--tumor_type` | — | Tumor type for AMP tiers: CIViC disease names joined by `\|` (whole run) |
| `--tumor_types` | — | TSV `id<TAB>terms` (patient, tumor or pair); overrides `--tumor_type` |

### Capture, MANE, panels

| Parameter | Default | Description |
|---|---|---|
| `--panel_bed` | — | WES capture BED: TMB territory (with `--mane_cds_bed`) and TERT promoter note |
| `--mane_cds_bed` | — | MANE Select CDS BED; without it the TMB territory is the whole capture (warning) |
| `--gene_coords_bed` | — | Gene coordinates of the capture, for the ASCAT / CNVkit interpretation |
| `--germline_panel` | `assets/gene_panels/germline_cancer_predisposition.txt` | Germline SNV filter (105 genes) |
| `--somatic_snv_panel` | `assets/gene_panels/somatic_snv_panel.txt` | Core somatic SNV filter (60 genes) |
| `--somatic_cnv_panel` | `assets/gene_panels/somatic_38_oncogenes.txt` | Somatic CNV panel (38 genes) |
| `--second_hit_panel` | `assets/gene_panels/somatic_track_second_hit.txt` | Second-hit somatic track (49 genes; `null` = off) |
| `--tso500_panel` | `assets/gene_panels/somatic_track_tso500.txt` | TSO500 pan-tumor track (466 genes; `null` = off) |
| `--cnv_gatk_genes` | `assets/gene_panels/somatic_cnv_gatk_genes.txt` | Genes annotated by the GATK somatic CNV step (113) |
| `--pharmaco_variants` | `assets/pharmacogenomics_dpyd_variants.txt` | rsIDs extracted for pharmacogenomics |

### Germline CNV (GATK gCNV)

| Parameter | Description |
|---|---|
| `--gcnv_interval_list` | Interval list the cohort model was trained on (FilterIntervals output) |
| `--gcnv_ploidy_model` | Contig-ploidy model directory (DetermineGermlineContigPloidy COHORT) |
| `--gcnv_pon_dir` | gCNV cohort model directory (GermlineCNVCaller COHORT) |
| `--gcnv_fasta`, `--gcnv_dict` | hg38 FASTA (+ `.fai`) and dictionary; the FASTA also serves the consensus and SpliceAI |
| `--gcnv_mane_txt` | MANE exons (`exons_mane_v1.2.txt`); also GATK somatic CNV and the coverage fallback |
| `--gcnv_lowmap_bed` | Low-mappability BED (`scripts/build_lowmap_panel.sh`) |
| `--cram_fasta` | FASTA (+ `.fai`, `.dict`) the CRAMs were written with; only when a sample has a CRAM and no BAM |

### Somatic CNV (GATK)

| Parameter | Default | Description |
|---|---|---|
| `--somatic_pon_hdf5` | — | Read-count PoN (CreateReadCountPanelOfNormals) |
| `--somatic_interval_list` | — | Preprocessed intervals of the PoN |
| `--common_snps_vcf` | — | Common SNPs (bgzip + tabix) for allelic counts |
| `--somatic_fasta`, `--somatic_dict` | — | hg38 main chromosomes FASTA and dictionary |
| `--somatic_min_contig_length` | `46709983` | Smallest contig in the GATK plots |

### HRD

| Parameter | Default | Description |
|---|---|---|
| `--hrd_organ` | `BREAST` | HRProfiler model: `BREAST` or `OVARIAN` (the only validated tumor types; applied to any tumor type, as the methodology states) |

Thresholds are not command-line parameters: set them in `conf/filters.config`.

### Profiles

`docker` (runs as the host user), `singularity`, `test` (stub test), `mirror`
(digest-pinned `ghcr.io/genuit-lbmu/*` images, `conf/containers.mirror.config`),
`debug`.

---

## Skip flags

Set them in the params file (`run_full.sh` takes no other options) or on the
`nextflow run` command line. A switched-off analysis is reported as "not analysed in
this run".

| Flag | Effect |
|---|---|
| `--skip_germline_snv true` | No germline SNV filtering / SpliceAI / classification; gnomAD / ClinVar / REVEL no longer required |
| `--skip_pharmaco true` | No pharmacogenomics extraction (with `--skip_germline_snv` too, the germline subworkflow does not run) |
| `--skip_germline_cnv true` | No GATK gCNV (no normal BAM / CRAM, or no cohort model) |
| `--skip_somatic_cnv_gatk true` | No GATK somatic CNV (no tumor BAM / CRAM, or no PoN) |
| `--skip_tmb true` | No TMB |
| `--skip_signatures true` | No SigProfiler (below 100 SBS / 50 ID the signature calls are already "not assessable") |
| `--skip_hrd true` | No HRProfiler (pairs without ASCAT segments are left out automatically) |
| `--skip_strelka_consensus true` | Mutect2 PASS only for TMB, signatures and HRD (no Strelka2 intersection) |
| `--skip_coverage true` | No per-gene coverage |
| `--skip_report true` / `--skip_pdf_report true` | No HTML / no PDF (the report model is built while either is on) |
| `--skip_somatic_snv true` | No somatic SNV tracks (core, second hit, TSO500), no CIViC / AMP tiers; `--civic_cache` no longer required. The consensus, TMB, signatures and HRD still run |
| `--skip_somatic_cnv true` | No ASCAT / CNVkit interpretation (somatic CNV section and LOH second hits) |

---

## Filtering thresholds (`conf/filters.config`)

Every filter, ACMG weight and biomarker threshold is in this file; `nextflow.config`
exports them to the scripts (`FILT_*`, `QC_*`, `ACMG_*`, `BIO_*`, `SIG_*`) and the report
model records the values in force.

| Section | Keys (default) | Use |
|---|---|---|
| `germline_snv` | `min_dp=5`, `vaf_min_high=0.20`, `vaf_min_moderate=0.20`, `vaf_floor_min_dp=20`, `vaf_max=0.75`, `multi_hit_distance=30`, `multi_hit_vaf=0.40`, `multi_hit_min_gq=20`, `noncoding_splice=true` | Depth; low-VAF floor (only when DP ≥ 20, else QC flag "DP<20"); zygosity (`vaf_max`, not a filter); nearby multi-hit artefacts (low VAF, FILTER ≠ PASS or low GQ); non-coding splice candidates kept with SpliceAI ≥ 0.2 or ClinVar P/LP |
| `germline_cnv` | `min_qa=4`, `lowmap_min_frac=0.5`, `exclude_genes=AR,BRAF,CCND1,CCND2,CCND3,CCNE1,CDK6,EGFR,ERBB2,FGFR1,FGFR2,FGFR3,KRAS,MDM2,MDM4,MYC` | gCNV quality; low-mappability fraction listed apart; panel genes left out of germline CNV (somatic oncogenes) |
| `somatic_snv` | `min_dp=10`, `min_af=0.05`, `min_alt_reads=10`, `max_pop_af=0.01` | Tumor DP (all tracks, TMB); TMB VAF; signatures ALT reads; Tier IV population frequency (gnomAD FAF95, else AF) |
| `tmb` | `min_alt_reads=3` | TMB ALT reads |
| `ascat` | `min_purity=0.20`, `min_goodness_of_fit=90`, `ploidy_note=2.7` | Report alerts: low purity, poor ASCAT fit; WGD / high-ploidy note on the absolute AMP threshold |
| `somatic_cnv` | `amp_threshold=5` | AMP when CN ≥ 5 |
| `acmg` | `clinvar_min_stars=2`, `frequency_metric='faf95'`, `ba1_af=0.05`, `pm2_af=0.001`, `bs1_af=0.01`, `vcep_specs=true`, `pvs1_score=8`, `pvs1_strong_score=4`, `pvs1_moderate_score=2`, `pvs1_last_exon_min_lost=0.10`, `pm4_score=2`, `pm2_score=1`, `ps1_score=4`, `pm5_score=2`, `pm1_score=2`, `pm1_sup_score=1`, `bs1_score=-4`, `bs1_sup_score=-1`, `pvs1_not_applicable=ALK,AR,BRAF,…,RET`, REVEL PP3 0.644 / 0.773 / 0.932 and BP4 0.290 / 0.183 / 0.016, `spliceai_pp3=0.2`, `spliceai_bp4=0.1`, `score_p=10`, `score_lp=6`, `score_lb=-1`, `score_b=-7` | ACMG classification (see [ACMG points](#acmg-points)) |
| `biomarkers` | `tmb_high=10.0`, `tmb_borderline=6.0`, `msi_high=3.5`, `msi_borderline=1.5`, `hrd_high=0.45`, `hrd_borderline=0.30` | Biomarker statuses (correlation, findings, both reports); `hrd_high` is also HRProfiler's threshold |
| `qc` | `low_cov_threshold=50` | Tumor coverage below which TMB, MSI, HRD, signatures and CNV are "Indeterminate" |
| `signatures` | `min_fraction=0.10`, `min_count=10`, `min_total_sbs=100`, `min_total_id=50`, `mmr_sbs=SBS6,SBS15,SBS20,SBS21,SBS26,SBS44`, `mmr_id=ID1,ID2,ID7`, `hrd_sbs=SBS3,SBS8` | When a signature counts as MMR / HRD evidence |
| `consensus` | `min_callers=2` | 2 = Mutect2 ∩ Strelka2; 1 = every Mutect2 PASS (any other value stops the run) |

### ACMG points

Germline variants without a final ClinVar class are scored with the ClinGen SVI
point system (Tavtigian et al. 2020): P >= 10, LP 6..9, VUS 0..5, LB -1..-6, B <= -7.

| Criterion | Points | Rule |
|---|---|---|
| PVS1 | +8 | nonsense / frameshift outside the last exon, splice donor/acceptor +-1/2, other HIGH |
| PVS1_St | +4 | nonsense / frameshift in the last exon, >= 10% of the protein lost |
| PVS1_Mod | +2 | last exon with < 10% lost, or start_lost |
| PVS1_NA | 0 | gene in `pvs1_not_applicable` (loss of function is not the cancer mechanism: oncogenes, POLD1/POLE, HOXB13, ...) |
| PM4 | +2 | stop_lost |
| PM2_Sup | +1 | gnomAD AF < `pm2_af` and FAF95 < `pm2_af` (rare in every ancestry group; not for sites failing gnomAD QC) |
| PP3 / BP4 | +-1 / +-2 / +-4 | REVEL >= 0.644 / 0.773 / 0.932 ; <= 0.290 / 0.183 / 0.016 (Pejaver et al. 2022); only when the MANE effect is missense |
| PP3_Sup(splice) | +1 | SpliceAI max delta >= 0.2 (not added to a REVEL PP3; the stronger one counts) |
| BP4 (missense) | | needs SpliceAI <= 0.1 when SpliceAI ran (Walker et al. 2023) |
| PS1 | +4 | missense with the same amino-acid change as a ClinVar P/LP missense (>= 2 stars, not at an exon edge); not if SpliceAI >= 0.2 |
| PM5 | +2 | missense at a residue where another ClinVar P/LP missense exists (when PS1 does not apply) |
| PM1 / PM1_Sup | +2 / +1 | VCEP genes only: hotspot / critical domain (see below); PP3 + PM1 at most +4 (Pejaver et al. 2022) |
| BS1 | -4 | gnomAD FAF95 > `bs1_af`, or the VCEP gene threshold |
| BS1_Sup | -1 | VCEP gene threshold (BRCA1, BRCA2, PTEN, ENG) |
| BA1 | stand-alone B | gnomAD FAF95 > `ba1_af`, or the VCEP gene threshold (also removes the variant in the germline filter, except ClinVar P/LP with >= `clinvar_min_stars` stars) |

PP3 / BP4 are not applied to null variants (PVS1 already weighs them).

FAF95 is the maximum filtering allele frequency (95% CI) across gnomAD ancestry
groups (`fafmax_faf95_max_joint`): a variant common in one group, e.g. admixed
American, is benign even when its overall AF is low (BRCA2 p.Ile2490Thr: AF 0.5%,
FAF95 11%). Without FAF95 (legacy mode, outside the resource regions) the overall
AF is used; `frequency_metric = 'af'` restores AF-only behaviour. Reports show
"AF / FAF 0.113".

#### ClinGen VCEP gene specifications

For 17 genes the frequency rules (and PM1 for 7 of them) follow the Released
specifications of the ClinGen Variant Curation Expert Panels, copied from the
[CSpec registry](https://cspec.genome.network) on 2026-10-23 into
`bin/acmg_shared.py` (`VCEP_SPECS`, `VCEP_PM1`). `vcep_specs = false` turns them
off. The `ACMG_spec` column of the annotated TSV names the rule set used.
FAF = gnomAD grpmax FAF95 (0 when gnomAD computes none, e.g. AC = 1); grpmax
AF / AC / AN come from the gnomAD resource built with the grpmax columns (older
resources: general rules). gnomAD v4 grpmax and FAF exclude the bottlenecked
groups (ASJ, FIN, AMI, remaining).

| Gene | Spec | BA1 | BS1 | BS1_Sup | PM2_Sup |
|---|---|---|---|---|---|
| APC | GN089 v2.1 | FAF >= 0.1% | FAF >= 0.001% | | AF <= 0.0003% if AC > 1; AF < 0.001% if AC <= 1 |
| ATM | GN020 v1.6 | FAF > 0.5% | FAF > 0.05% | | grpmax AF <= 0.001%, or AC = 1 |
| BRCA1, BRCA2 | GN092 / GN097 v1.2 | FAF > 0.1% | FAF > 0.01% | FAF > 0.002% | absent (AC = 0) |
| CDH1 | GN007 v3.1 | FAF >= 0.2% | FAF >= 0.1% | | AF <= 1/100,000; grpmax AF <= 1/50,000 if grpmax AC >= 2 |
| DICER1 | GN024 v1.4 | grpmax AF > 0.3% (AN > 2,000, AC >= 5) | grpmax AF > 0.03% (idem) | | AF < 0.0005% and grpmax AC <= 1 |
| ENG | GN136 v1.1 | FAF >= 1% | FAF > 0.2%, or BS1_Sup with >= 2 homozygotes | FAF > 0.08% | AC < 6, or grpmax AF < 0.004% |
| KRAS, BRAF | GN044 / GN049 v2.3 | FAF >= 0.05% | FAF >= 0.025% | | absent |
| MLH1, MSH2 | GN115 / GN137 v2.0 | FAF >= 0.1% | FAF >= 0.01% | | AF < 0.002% |
| MSH6 | GN138 v2.0 | FAF >= 0.22% | FAF >= 0.022% | | AF < 0.002% |
| PMS2 | GN139 v2.0 | FAF >= 0.28% | FAF >= 0.028% | | AF < 0.002% |
| PALB2 | GN077 v1.3 | FAF > 0.1% | FAF > 0.01% | | grpmax AF <= 1/300,000, or AC = 1 |
| PTEN | GN003 v3.2 | FAF > 0.056% | FAF >= 0.0043% | FAF >= 0.00043% | AF < 0.001%; grpmax AF < 0.002% if grpmax AC >= 2 |
| TP53 | GN009 v2.4 | FAF >= 0.1% (grpmax AN >= 2,000, AC >= 2) | FAF >= 0.03% (idem) | | AF < 0.003%; grpmax AF < 0.004% if grpmax AC >= 2 |
| VHL | GN078 v1.1 | FAF >= 0.0156% | FAF >= 0.00156% | | FAF <= 0.000156% (or none) |

PM1 (missense, MANE protein numbering): TP53 codons 175, 245, 248, 249, 273,
282 or cancerhotspots.org >= 10 for the same change (Moderate), 2-9
(Supporting); PTEN 90-94, 123-130, 166-168; DICER1 metal-ion residues S1344,
E1705, D1709, D1713, G1809, D1810, E1813 (Moderate) and the RNase IIIb domain
1682-1846 (Supporting); ENG S278, F282 and cysteines C207, C350, C363, C382,
C394, C412, C516, C549, C582; VHL cancerhotspots.org >= 10 (Moderate), < 10
(Supporting); KRAS P-loop 10-17, SW1 25-40, SW2 57-64, SAK 145-156; BRAF exons
6 and 11, P-loop 459-474, CR3 594-627. Hotspot counts: cancerhotspots.org v2
(Chang et al. 2018; `scripts/build_vcep_hotspots.py`), check its terms of use.

Not automated: the non-cancer gnomAD v2.1 / v3.1 subsets and read depth of the
BRCA1/2 PM2 rule (AC = 0 in gnomAD v4.1 is used), the founder-variant exclusion
of the MMR BA1 / BS1 rules (ClinVar P/LP with >= 2 stars keeps its class and is
never removed by BA1), the TP53 CHIP recalculation by allele balance, and the
VHL germline hotspot table (not in the CSpec API). Not included: CHEK2, RAD51C,
MUTYH, SDHB, SMAD4 (specifications not released), NF1 (no rules in the API),
PIK3CA (brain-malformation / mosaic specification) and LZTR1 (RASopathy, not
schwannomatosis). Other VCEP rules (gene PVS1 trees, other PP3 tools, PS3 / PS4)
are not applied.

Rare missense variants cannot reach LP from in silico and frequency evidence
alone (maximum PP3_St + PM2_Sup = 5); a new null variant in a loss-of-function
gene is LP (PVS1 + PM2_Sup = 9). The penultimate-exon NMD rule (last 50 nt) is
not applied.

---

## Somatic consensus calling (Mutect2 + Strelka2)

The pipeline runs two independent somatic variant callers through nf-core/sarek and intersects their PASS calls to produce a high-confidence consensus VCF:

```
Mutect2 VCF (snpEff-annotated) ──┐
                                  ├─→ bcftools norm → bcftools isec -n+2 -w1 → consensus.vcf.gz
Strelka2 SNVs + Indels (PASS) ──┘
```

**How it works:**

0. The Mutect2 VCF is rewritten with its samples as normal,tumor (`MUTECT2_SAMPLE_ORDER`, from the `##normal_sample` / `##tumor_sample` header lines; the run stops if they are missing), so every somatic step reads the tumor as `GEN[1]`
1. Both VCFs are filtered for PASS variants (Mutect2: FILTER exactly `PASS` or exactly `haplotype`, i.e. phased with a filtered neighbour but passing its own filters)
2. Indels are left-aligned, multi-allelic sites split and MNPs atomized into SNVs (`bcftools norm -m-both`, then `-a --atom-overlaps .`), so a Mutect2 MNP such as `CC>TT` matches the two SNVs Strelka2 reports
3. Only variants called by **both** callers are kept (`bcftools isec -n+2`; `filters.consensus.min_callers = 1` keeps every Mutect2 PASS instead, Strelka2-private calls never enter)
4. The consensus VCF preserves Mutect2 FORMAT fields (DP, AD, AF) and snpEff annotations (`-w1`)

**What uses the consensus VCF:**

| Analysis | Input | Why consensus matters |
|----------|-------|----------------------|
| **TMB** | Consensus VCF | Removes false-positive calls that inflate TMB |
| **Mutational Signatures** | Consensus VCF | Cleaner mutation profiles = more accurate signature decomposition |
| **HRD Profiler** | Consensus VCF + ASCAT segments | Mutation features (N[C>T]G, N[C>G]T, microhomology deletions) depend on accurate somatic calls |

**What still uses Mutect2-only:**

| Analysis | Reason |
|----------|--------|
| **Somatic SNV report table** | Needs full snpEff + CIViC annotations (only in Mutect2 path) |
| **CIViC annotation** | Requires snpEff gene/protein annotations from Mutect2 |

The somatic SNV report table includes a **Callers** column indicating whether each variant was confirmed by both callers (`M+S`) or only by Mutect2 (`M`). Variants with VAF < 5% are flagged with a warning badge.

**Fallback:** per pair. A pair without Strelka2 output (e.g. a sarek run without `--tools strelka,manta`) uses Mutect2 PASS for TMB, signatures and HRD; `--skip_strelka_consensus true` does the same for every pair. When Strelka2 output exists, the consensus needs `--gcnv_fasta` (+ `.fai`) for `bcftools norm`; without it the run stops at start with that message.

**Output:** `somatic/consensus/*.caller_stats.tsv` shows variant counts per caller and intersection size.

---

## Gene panels

Every gene list comes from one table, `assets/gene_panels/genes.tsv` (543 genes, one
row per gene; sources in `genes.manifest.txt`):

| Column | Genes | Generated file / use |
|---|---|---|
| `germline` | 105 | `germline_cancer_predisposition.txt`: germline SNV / CNV, LOH second hits |
| `second_hit` | 86 | germline genes where a somatic loss-of-function event completes biallelic inactivation; outside the core: `somatic_track_second_hit.txt` (49) |
| `somatic_core` | 60 | `somatic_snv_panel.txt`: somatic SNV report |
| `tso500` | 523 | `somatic_tso500_panel.txt` (Illumina TruSight Oncology 500 DNA content); outside the core: `somatic_track_tso500.txt` (466) |
| `somatic_cnv` | 38 | `somatic_38_oncogenes.txt`: somatic CNV report |
| `mmr`, `hrr` | 8, 26 | `bin/gene_sets.py` (`MMR_GENES`, `HRR_GENES`): correlation and findings |

Other generated files: `germline_second_hit.txt` (86, the `second_hit` column) and
`somatic_cnv_gatk_genes.txt` (113, germline + somatic CNV genes annotated by the GATK
somatic CNV step, `--cnv_gatk_genes`).

Edit `genes.tsv`, then run `python3 scripts/build_gene_lists.py` (`--check` reports
stale files). Custom lists can be passed to the filtering steps: `--germline_panel`,
`--somatic_snv_panel`, `--somatic_cnv_panel`, `--second_hit_panel`, `--tso500_panel`
(one HGNC symbol per line; `null` disables a track). The same lists reach the report
model and the germline–somatic correlation.

**Somatic SNV tracks.** The same chain (panel + PASS + impact + DP, MANE, gnomAD / Tier IV,
CIViC, AMP) runs per track; files carry the track suffix
(`<pair>.second_hit.somatic_snvs.tsv`, `<pair>.tso500.somatic_snvs.tsv`):

- **core** (60 genes): somatic SNV tables of HTML and PDF, as before.
- **second hit** (germline genes outside the core): only the germline–somatic
  correlation and the findings cards read it, for genes with a germline P/LP variant;
  it is never shown in the somatic tables. Activating, risk-allele and non-tumor germline
  genes are excluded (ALK, CCND1-3, CCNE1, CDK4, CFTR, CTRC, EGFR, HOXB13, KIT, MDM2,
  MDM4, MET, MITF, MYC, PDGFRA, RET, SPINK1).
- **pan-tumor TSO500** (genes outside the core): Tier I/II join the actionable variants
  and the findings (marked TSO500), Tier III go to a separate "Pan-tumor (TSO500)" table
  (collapsed in HTML), Tier IV to the population polymorphisms. RYBP and TERC have no
  MANE Select transcript in MANE v1.2 (the snpEff database), so they are not annotated.
  TERT promoter hotspots (C228T chr5:1295113, C250T chr5:1295135) are
  `upstream_gene_variant` (MODIFIER) and are not called by the impact filter; when they
  lie outside `--panel_bed` (Agilent S07604514 ends ~120 bp before them, depth 0–2) the
  HTML and PDF say "TERT promoter: not assessable", so absence is not read as negative.

Gene context for the TSO500 genes (`assets/gene_annotations/somatic_gene_context.tsv`,
`scripts/build_gene_context.py`): ROLE from CancerMine (CC0) only when the literature is
one-sided (>= 5x citations, >= 20), GENE_SUMMARY from the CIViC gene description (CC0) or
the NCBI Gene RefSeq summary; the 60 core genes keep their curated text. COSMIC CGC and
OncoKB are not used (commercial use needs a licence).

### Germline SNV + CNV panel (105 genes)

**File:** `germline_cancer_predisposition.txt` | **Used by:** Germline SNV filtering (HIGH/MODERATE impact; splice-region / synonymous / intronic with SpliceAI >= 0.2) + Germline CNV reporting (panel minus `filters.germline_cnv.exclude_genes`; `<DEL>`/`<DUP>`, QA>=4)

Sources: ACMG SF v3.2, ESMO guidelines, published WES cancer panels.

<details><summary>Full gene list (click to expand)</summary>

AIP, ALK, APC, AR, ATM, ATR, AXIN2, BAP1, BARD1, BLM, BMPR1A, BRAF, BRCA1, BRCA2, BRIP1, BUB1B, CASR, CCND1, CCND2, CCND3, CCNE1, CDC73, CDH1, CDK4, CDK6, CDKN1B, CDKN2A, CFTR, CHEK2, CTNNA1, CTRC, DICER1, EGFR, ENG, EPCAM, ERBB2, FANCC, FANCD2, FANCM, FGFR1, FGFR2, FGFR3, FH, FLCN, GREM1, HOXB13, KIF1B, KIT, KRAS, LZTR1, MAX, MBD4, MDM2, MDM4, MEN1, MET, MITF, MLH1, MLH3, MRE11, MSH2, MSH3, MSH6, MUTYH, MYC, NBN, NF1, NF2, NTHL1, PALB2, PDGFRA, PIK3CA, PMS2, POLD1, POLE, POT1, PRKAR1A, PTCH1, PTEN, RAD51C, RAD51D, RB1, RECQL, RECQL4, RET, RNF43, RPS20, SDHA, SDHAF2, SDHB, SDHC, SDHD, SMAD4, SMARCA4, SMARCB1, SMARCE1, SPINK1, STK11, SUFU, TMEM127, TP53, TSC1, TSC2, VHL, WRN

</details>

### Somatic SNV panel (60 genes)

**File:** `somatic_snv_panel.txt` | **Used by:** Mutect2 somatic SNV/indel filtering (HIGH/MODERATE impact)

Covers actionable oncogenes, tumor suppressors, and DNA repair genes.

<details><summary>Full gene list (click to expand)</summary>

APC, AR, ARAF, ATM, ATR, B2M, BLM, BRAF, BRCA1, BRCA2, BRIP1, CDH1, CDKN1B, CDKN2B, CDKN2C, CHEK2, CTNNB1, EGFR, ERBB2, FBXW7, FGFR2, GNAS, IDH1, IDH2, KIT, KRAS, MAP2K1, MEN1, MLH1, MLH3, MSH2, MSH6, MTOR, NF1, NOTCH3, NRAS, NTHL1, NTRK1, NTRK3, PALB2, PDGFRA, PDGFRB, PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RET, RNF43, ROS1, SDHA, SDHD, SETD2, SMAD4, SMARCA4, STK11, TP53, WRN

</details>

### Somatic CNV panel (38 genes)

**File:** `somatic_38_oncogenes.txt` | **Used by:** ASCAT/CNVkit somatic CNV reporting (HOMDEL CN=0, AMP CN>=5, 2nd-hit LOH)

<details><summary>Full gene list (click to expand)</summary>

APC, B2M, BLM, BRAF, BRCA1, BRCA2, CTNNB1, EGFR, EPCAM, ERBB2, FBXW7, KRAS, MAP2K1, MEN1, MET, MLH1, MLH3, MSH2, MSH3, MSH6, MTOR, MYC, NOTCH3, NRAS, PALB2, PIK3CA, PIK3R1, PMS2, POLD1, POLE, PTEN, RECQL4, RNF43, SDHA, SDHD, SMAD4, TP53, WRN

</details>

---

## Output

```
clinical_results/
├── reporting/
│   ├── <pair>.report_model.json            ← the report model (every value both reports show)
│   ├── <pair>_clinical_report.pdf          ← clinical report
│   ├── <pair>_molecular_report.html        ← laboratory report
│   └── <pair>_molecular_report.txt         ← plain-text summary
├── germline/
│   ├── snv/filtered/<normal>.filtered.vcf.gz             ← panel + impact + depth
│   ├── snv/extract/
│   │   ├── <normal>.germline_snvs.tsv                    ← candidates kept by the filters
│   │   ├── <normal>.germline_all_filtered.tsv            ← every rejected variant with FILTER_REASON
│   │   └── <normal>.germline_snvs_acmg_annotated.tsv     ← candidates with ACMG class, criteria, spec
│   ├── pharmaco/<normal>.pharmaco.tsv                    ← DPYD / UGT1A1 genotypes
│   └── cnv_gatk/
│       ├── counts/  ploidy/  calls/                      ← gCNV case-mode intermediates
│       ├── postprocessed/<normal>_{denoised.tsv,intervals.vcf.gz,segments.vcf.gz}
│       └── annotated/<normal>_segments_annotated.tsv     ← MANE exons, LOWMAP, LOWMAP_SEG
├── somatic/
│   ├── consensus/<pair>.consensus.vcf.gz(.tbi), <pair>.caller_stats.tsv
│   ├── snv/filtered/<pair>{,.second_hit,.tso500}.filtered.vcf.gz
│   ├── snv/civic/<pair>{,.tso500}.civic.vcf.gz
│   ├── snv/extract/<pair>{,.second_hit,.tso500}.somatic_snvs.tsv   ← Callers, gnomAD, POP_COMMON
│   ├── snv/civic_extract/<pair>{,.tso500}.civic.exploded.tsv       ← CIViC evidence records
│   ├── cnv/ascat_cnvkit/<pair>_somatic_cnv_report.{tsv,txt}        ← HOMDEL / HETDEL / AMP / LOH + ASCAT metrics
│   ├── cnv_gatk/
│   │   ├── counts/  model_segments/<patient>/  called_segments/
│   │   ├── plots/<patient>/{CR_plots,Modeled_segments_plots}/
│   │   └── annotated/<patient>/<tumor>_annotated.tsv               ← MANE genes, CN with ASCAT purity
│   ├── tmb/<pair>.tmb.tsv, tmb/filtered/<pair>.tmb_filter.vcf
│   ├── signatures/filtered_vcfs/<pair>.sigfilter.vcf
│   ├── signatures/<pair>/{matrix_profiles,assignment}/             ← SigProfiler (COSMIC v3.4)
│   └── hrd/<pair>_hrd_predictions.tsv, hrd/hrd_results/
├── qc/gene_coverage/<sample>_gene_coverage.tsv
└── pipeline_info/execution_{timeline,report,trace}_<date>.*, pipeline_dag_<date>.html
```

The MANE / gnomAD / ClinVar / REVEL / SpliceAI annotated VCFs are intermediate (work
directory only).

---

## Reproducibility

- Every resource is a local, dated file (gnomAD v4.1, ClinVar release, REVEL v1.3,
  CIViC snapshot, MANE v1.2, cohort CNV models); no annotation depends on a web service
  at run time. Only the clinical-trial list is a live query (dated in the report).
- The reports record the pipeline commit, the sarek release, the tool versions and the
  resource files; `report_model.json` keeps every value with the thresholds in force.
- `run_full.sh` pins the sarek release; the custom images are versioned.
- Report inputs are staged per pair and a hash of `bin/` and `assets/` is part of them,
  so `-resume` reuses a report only when nothing it depends on changed.

---

## Tests

```bash
python -m pytest -q tests/                     # unit tests of bin/ (synthetic data)
bash tests/stub/make_fixture.sh                # fake sarek output, samplesheet, empty resources
nextflow run . -profile test -stub --outdir results_test   # NXF_SYNTAX_PARSER=v1 with Nextflow >= 26
```

The unit tests (pytest; `reportlab` and `pymupdf` installed) cover ACMG classification,
MSI status, signature calls, AMP tiers, CNV harmonisation, the report model (metrics,
germline CNV segments, signatures, CPIC, identification) and the report helpers, and
check that every report script loads. The stub run executes the whole workflow with the
`stub:` block of every module (no tools, no containers; `bgzip` / `tabix` needed),
which checks channels, joins and report inputs. CI (`.github/workflows/ci.yml`) runs
both on every push, the stub run on Nextflow 24.04.0 and the latest release.

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `Config parsing failed` / `Statements cannot be mixed with script declarations` (Nextflow ≥ 26) | `export NXF_SYNTAX_PARSER=v1` (`run_full.sh` does it) |
| `Session [xxx] lock file already acquired` | Another run from the same launch directory; `rm -f <launch-dir>/.nextflow/cache/*/db/LOCK` only if none is running |
| `Can't lock file: .nextflow/history` | Phase 2 of `run_full.sh` launches from the repository directory: run one at a time, or use `nextflow run` from separate directories |
| A report did not change after a code fix | It does with `-resume` (the code hash is a report input); for other tasks, `bin/` is not part of the task hash: remove their work directory or run without `-resume` |
| Out of memory (GATK gCNV) | `--max_memory` (caps every task) or lower `process_high` in `conf/base.config` (56 GB) |
| gCNV: model version mismatch | Use `broadinstitute/gatk:4.6.2.0` (gcnvkernel 0.9) with a model built by the same version |
| gCNV: interval mismatch | `--gcnv_interval_list` must be the interval list the cohort model was trained on |
| `--gcnv_fasta ... is required for the Mutect2 + Strelka2 consensus` | Set `--gcnv_fasta` (+ `.fai`) or `--skip_strelka_consensus true` |
| `... is a CRAM: set --cram_fasta` | Give the FASTA the CRAMs were written with, or skip the GATK CNV steps |
| CIViC annotation fails | The task fails (it never passes an unannotated VCF). Check `--civic_cache` and that it was built with the image's civicpy version (`scripts/build_civic_cache.sh`) |
| `annotate_cnvs_complete_fixed.R wrote no ...` | The GATK somatic CNV annotation failed; see the task's `.command.err` |
| Clinical trials "API unavailable" | No HTTPS access to clinicaltrials.gov from the host |
| `Join mismatch` in sarek (Mutect2) | `--tools` must be on the sarek command line (`run_full.sh` does it) |

---

## Citation

If you use `nf-sarek-clinical` in your research, please cite:

- **nf-core/sarek**: Garcia M, et al. (2020). *Genome Medicine*.
- **CIViC/civicpy** (v5.4.0): Griffith M, et al. (2017). *Nature Genetics*.
- **SigProfilerAssignment**: Bergstrom EN, et al. (2022). *Cell Genomics*.
- **HRProfiler**: Chopra N, et al. (2024). [github.com/FlorPio/HRProfiler-workflow](https://github.com/FlorPio/HRProfiler-workflow).
- **GATK CNV workflow**: [github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic).
- **ASCAT**: Van Loo P, et al. (2010). *PNAS*.
- **GATK** (CNV, gCNV): Van der Auwera GA, O'Connor BD (2020). *Genomics in the Cloud*. O'Reilly.
- **gnomAD v4**: Chen S, et al. (2024). *Nature*; **ClinVar**: Landrum MJ, et al. (2018). *Nucleic Acids Res*.
- **REVEL**: Ioannidis NM, et al. (2016). *Am J Hum Genet* (data ODbL); **SpliceAI**: Jaganathan K, et al. (2019). *Cell*.
- **MANE**: Morales J, et al. (2022). *Nature*; **ACMG/AMP**: Richards S, et al. (2015). *Genet Med*; Tavtigian SV, et al. (2020). *Hum Mutat*; Pejaver V, et al. (2022). *Am J Hum Genet*.
- **AMP/ASCO/CAP somatic tiers**: Li MM, et al. (2017). *J Mol Diagn*; **CPIC** guidelines (cpicpgx.org).
- **cancerhotspots.org**: Chang MT, et al. (2018). *Cancer Discov*; **CancerMine**: Lever J, et al. (2019). *Nat Methods* (CC0); **CIViC** data CC0.
- **GIAB genome stratifications** v3.3; **COSMIC** mutational signatures v3.4 (cancer.sanger.ac.uk/signatures).

## License

MIT — see [LICENSE](LICENSE). Resource licences: see Citation (REVEL ODbL; CIViC, CancerMine CC0; cancerhotspots.org terms of use).
