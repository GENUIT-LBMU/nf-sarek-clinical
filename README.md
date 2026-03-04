# nf-sarek-clinical

**Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exome outputs.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

![Pipeline overview](metro_map_sarek_clinical.png)

---

## Pipeline modes

`nf-sarek-clinical` can be launched in three ways depending on what has already been processed:

| Mode | Command | When to use |
|------|---------|-------------|
| **Full** | `./run_full.sh` | Start from FASTQ — runs nf-core/sarek then the clinical report |
| **Post-sarek** | `./run_full.sh --skip-sarek` | nf-core/sarek already completed; run clinical report only |
| **Clinical only** | `nextflow run main.nf` | Explicit Nextflow call with custom parameter overrides |

Any module can be skipped independently with `--skip_<module> true` (see [Skip flags](#skip-flags)).

---

## Report sections

Each section below corresponds to a block in the HTML molecular report and a chapter in the clinical PDF.

### 0 · Major Clinical Findings

The first page of the PDF. Automatically populated when any of the following thresholds is crossed:

- **HRD** — HRProfiler probability ≥ threshold (default 0.45). Includes BRCA1/2 germline/somatic hits and detected HRD-associated signatures (SBS3, SBS8).
- **MMR deficiency** — MMR gene P/LP germline variant, biallelic somatic MMR event, MSI score > 10 %, or MMR mutational signatures (SBS6/15/20/21/26/44, ID1/ID2).
- **High TMB** — ≥ 10 mut/Mb.
- **Relevant CNVs** — amplification or homozygous deletion of an oncogene or tumour suppressor in the somatic 38-gene panel.

---

### 1 · Germline SNVs

Hereditary variant analysis restricted to a configurable gene panel (default: 90 ACMG SF v3.2 hereditary cancer genes).

**Filtering pipeline:**
1. dbNSFP annotation via SnpSift — adds 9 in silico predictors, REVEL, ClinVar, gnomAD allele frequencies.
2. Panel filter — only variants in the configured gene list are retained.
3. Impact filter — coding variants (missense, stop, frameshift, splice) only.
4. Depth filter — DP ≥ 5 (configurable with `--germline_min_dp`).
5. gnomAD rescue — variants with no AF in dbNSFP are queried against the Ensembl VEP REST API to recover gnomAD AF.
6. Population filter — gnomAD AF < 0.01 (BA1 benign excluded).
7. ACMG classification — each variant is scored P / LP / VUS / LB / B using a rule set that includes: ClinVar status, 9-tool in silico majority vote, REVEL ≥ 0.7 (PP3), PVS1 for LoF in LoF-intolerant genes, PM2 for ultra-rare variants.

**Output:** `*.germline_snvs.tsv` (P/LP/VUS only), `*.germline_snvs_acmg_annotated.tsv` (all candidates with scores), `*.germline_all_filtered.tsv` (all rejected variants with reason).

---

### 2 · Pharmacogenomics

Identifies actionable germline variants in five drug-metabolism genes relevant to oncology treatment selection:

| Gene | Drug relevance |
|------|---------------|
| **DPYD** | 5-fluorouracil toxicity (DPYD*2A, *13, HapB3, rs67376798) |
| **UGT1A1** | Irinotecan toxicity (UGT1A1*28, *6) |
| **CYP2D6** | Tamoxifen metabolism |
| **CYP2C19** | Multiple drugs including antidepressants used in supportive care |
| **NUDT15** | Thiopurine toxicity |

Variants are extracted from the germline VCF by rsID lookup. **Output:** `*.pharmaco.tsv`.

---

### 3 · Germline CNVs (GATK gCNV)

Germline copy-number calling on the normal sample using a pre-built Panel-of-Normals (PoN). Runs in CASE mode against the cohort PoN — no tumour BAM required.

**Steps:** CollectReadCounts → DetermineGermlineContigPloidy → GermlineCNVCaller → PostprocessGermlineCNVCalls → MANE annotation (exon-level, bundled R script from [FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)).

**Filters applied during annotation:** number of points ≥ 5, quality score ≥ 10.

**Output:** `*_segments_annotated.tsv` — one row per segment with gene/exon overlap, copy-number call, and quality metrics.

---

### 4 · Somatic SNVs

Somatic variant interpretation from Mutect2 VCF output.

**Filtering pipeline:**
1. dbNSFP annotation — same 9-predictor panel as germline.
2. Panel filter — somatic SNV gene list (oncogenes, TSGs, DNA repair genes).
3. PASS + coding impact filter.
4. CIViC annotation — variants are matched against the CIViC clinical variant database ([civicpy](https://civicpy.org)) to retrieve clinical evidence items, drug associations, and evidence level.
5. AMP/ASCO/CAP tiering — each variant is assigned Tier IA / IB / IIC / IID / III based on clinical evidence and functional impact.

**Output:** `*.somatic_snvs.tsv`, `*.civic.exploded.tsv` (one row per CIViC evidence item).

---

### 5 · Somatic CNVs (ASCAT / CNVkit)

Interprets copy-number segments from nf-core/sarek's ASCAT and CNVkit outputs, restricted to a 38-gene actionable oncogene/TSG panel.

**Calls reported:**

| Call | Criterion |
|------|-----------|
| **HOMDEL** | CN = 0 (homozygous deletion) |
| **HETDEL** | CN = 1 (heterozygous loss / haploinsufficiency) |
| **AMP** | CN ≥ threshold (default 5, configurable with `--amp_threshold`) |
| **LOH** | Minor allele CN = 0 with total CN ≥ 2 (ASCAT only) |

Also detects **somatic biallelic events**: gene present in both a germline P/LP variant (Section 1) and a somatic second hit (LOH, HOMDEL, or second somatic SNV).

**Output:** `*_somatic_cnv_report.tsv`.

---

### 6 · Somatic CNVs (GATK)

Allele-specific, tumour-normal paired copy-number analysis using GATK's read-depth and allele-fraction modelling.

**Steps:** CollectReadCounts → DenoiseReadCounts (against somatic PoN) → CollectAllelicCounts → ModelSegments → CallCopyRatioSegments → PlotDenoisedCopyRatios → PlotModeledSegments → MANE annotation (bundled R script from [FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)).

Both tumour and normal BAMs are processed in parallel through steps 1–3; pairing occurs at ModelSegments. Diagnostic CR and modelled-segment plots are published alongside the annotated TSV.

**Output:** `*_annotated.tsv` — per-segment log2 copy ratio, allele fraction, call (AMP/DEL/NEUTRAL/LOH), and overlapping MANE gene/exon.

---

### 7 · TMB (Tumour Mutational Burden)

Counts somatic PASS variants meeting minimum quality thresholds (DP ≥ 10, AF ≥ 0.05, ≥ 10 alt reads, coding impact) and normalises by the effective capture size from the WES BED file.

TMB ≥ 10 mut/Mb triggers the "High TMB" Major Clinical Finding card.

**Output:** `*.tmb.tsv` — raw count, effective Mb, and TMB value.

---

### 8 · Mutational Signatures

Decomposes the somatic SNV catalogue into COSMIC v3.4 reference signatures using SigProfilerMatrixGenerator + SigProfilerAssignment.

**Matrices generated:** SBS96 (single-base substitutions), DBS (doublet-base), ID83 (small insertions/deletions).

**Output:** observed mutation matrices per type + assignment activity files (one row per COSMIC signature per sample) + 4 publication-quality PNG plots (SBS96 spectrum, SBS etiology, ID spectrum, ID etiology).

Signatures with non-zero activity are summarised in the PDF report; MMR-associated signatures (SBS6/15/20/21/26/44, ID1/ID2) feed into the Major Clinical Findings MMR card.

---

### 9 · HRD (Homologous Recombination Deficiency)

Runs [HRProfiler](https://github.com/FlorPio/HRProfiler-workflow) (Chopra et al., 2024) which integrates SBS signatures, indel signatures, and allele-specific CNV features (CNV48 basis) to compute an HRD probability score via a logistic regression model.

Two organ-specific models are available: `BREAST` (default) and `OVARIAN`.

**Output:** `*_hrd_predictions.tsv` — HRD probability and binary classification (HRD+ / HRD−) at the configured threshold.

---

### 10 · Gene Coverage

Per-gene depth statistics (mean, min, max) derived from mosdepth summary files produced by nf-core/sarek. Reported for both normal and tumour samples.

**Output:** `*_gene_coverage.tsv`.

---

### 11 · Germline–Somatic Correlation

Cross-section integration that identifies clinically significant two-hit events and biomarker associations:

- **Two-hit analysis** — P/LP germline variant (Section 1) + somatic second hit (somatic SNV, LOH, or HOMDEL in the same gene).
- **Biallelic somatic events** — two independent somatic hits in the same gene (detected in Section 5, no germline variant required).
- **Signature correlation** — detected SBS/ID mutational signatures correlated with the implicated pathways (MMR, HRD, tobacco, UV).
- **Biomarker matching** — germline hits in genes with active therapeutic trials cross-referenced with the CIViC database.

VUS variants are excluded from two-hit correlation (only P/LP used as germline first hits).

**Output:** `correlation_section.json` — consumed by both the HTML and PDF report generators.

---

## Reference files

All bioinformatics **tools** are bundled inside the Docker images. The **R annotation scripts** for GATK CNV (`exon_annotation.R`, `annotate_cnvs_complete_fixed.R`) are sourced from [FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic) and bundled in `bin/` — no manual installation required.

What you need to provide are the large **reference databases and PoN model files**.
Define all paths once in `full_pipeline.params.yml` (template in the repo root).

---

### 1 · dbNSFP

```bash
# Download from the Liu Lab (free registration):
# https://sites.google.com/site/jpopgen/dbNSFP  →  hg38 version
mkdir -p /path/to/refs/dbNSFP && cd /path/to/refs/dbNSFP
tabix -s 1 -b 2 -e 2 dbNSFP4.x.txt.gz   # index after download
```

### 2 · SnpSift

```bash
# Download snpEff package (includes SnpSift):
# https://pcingola.github.io/SnpEff/  →  snpEff_latest_core.zip
mkdir -p /path/to/refs/snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip -d /path/to/refs/snpEff
# SnpSift.jar → /path/to/refs/snpEff/snpEff/SnpSift.jar
```

### 3 · hg38 reference genome

```bash
# From GATK GCP bucket:
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta     /path/to/refs/hg38/hg38.fa
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai /path/to/refs/hg38/hg38.fa.fai
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict      /path/to/refs/hg38/hg38.dict

# Create main-chromosomes subset (required by somatic GATK CNV):
samtools faidx hg38.fa chr1 chr2 ... chrX chrY > main_chroms/hg38.main_chroms.fa
samtools faidx main_chroms/hg38.main_chroms.fa
gatk CreateSequenceDictionary -R main_chroms/hg38.main_chroms.fa
```

### 4 · WES capture BED

Kit-specific BED file for TMB calculation. Download from your kit vendor (e.g. [Agilent eArray](https://earray.chem.agilent.com/suredesign/)).

### 5 · MANE exon coordinates

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.4.summary.txt.gz
# Convert to 6-column format (chr, start, end, gene_symbol, exon, transcript_id)
# See FlorPio/gatk-cnv-somatic for the conversion script
```

### 6 · GATK CNV support files and PoN (build once)

All interval lists, common SNPs VCF, and PoN construction scripts are provided by:
**[github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)**

```bash
git clone https://github.com/FlorPio/gatk-cnv-somatic /path/to/gatk-cnv-somatic
# Follow that repo's README to generate:
#   exomas.interval_list            (germline FilterIntervals output)
#   preprocessed.interval_list      (somatic PreprocessIntervals output)
#   common_snps_preprocessed.vcf.gz (common SNPs from GATK bundle)
#   ploidy-model/                   (DetermineGermlineContigPloidy COHORT)
#   PoN_CNV/                        (GermlineCNVCaller COHORT)
#   wes-do-gc.pon.hdf5              (CreateReadCountPanelOfNormals)
```

> PoN files are institution-specific (built from your cohort of ≥ 30 normals). Rebuild required if GATK version changes (`broadinstitute/gatk:4.6.2.0`).

### Suggested directory layout

```
/path/to/refs/
├── dbNSFP/
│   ├── dbNSFP4.1a.txt.gz  (+ .tbi)
├── snpEff/
│   └── SnpSift.jar
├── hg38/
│   ├── hg38.fa  (+ .fai, .dict)
│   └── <kit>/<kit>_Regions.bed
└── cnv/
    ├── exons_mane.txt
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

## Running the pipeline

### Option A — Full pipeline from FASTQ (recommended)

#### Step 1 — Create the samplesheet

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
1229,XY,0,1229_N07,2,/path/to/reads/1229_N07_L2_1.fq.gz,/path/to/reads/1229_N07_L2_2.fq.gz
1229,XY,1,1229_T08,1,/path/to/reads/1229_T08_L1_1.fq.gz,/path/to/reads/1229_T08_L1_2.fq.gz
1229,XY,1,1229_T08,2,/path/to/reads/1229_T08_L2_1.fq.gz,/path/to/reads/1229_T08_L2_2.fq.gz
```

> Sample naming convention: `<patient>_N##` for normal, `<patient>_T##` for tumor.

#### Step 2 — Launch

```bash
cd /path/to/nf-sarek-clinical

./run_full.sh \
  --run        <RUN_NAME> \
  --samplesheet /path/to/runs/<RUN_NAME>/samplesheet.csv
```

#### Common flags

```bash
./run_full.sh --run <RUN_NAME> --samplesheet /path/to/samplesheet.csv --resume   # resume after interruption
./run_full.sh --run <RUN_NAME> --skip-sarek                                        # clinical report only
./run_full.sh --run <RUN_NAME> --skip-sarek --resume                               # resume clinical only
```

#### Output layout

```
/path/to/runs/<RUN_NAME>/
├── samplesheet.csv
├── results/            ← nf-core/sarek outputs
├── clinical_work/      ← nf-sarek-clinical work dir
└── clinical_results/   ← HTML + PDF reports  ← final output
```

---

### Option B — Clinical report only (explicit Nextflow call)

```bash
cd /path/to/nf-sarek-clinical

nextflow run main.nf \
  -profile docker \
  -work-dir /path/to/runs/MY_RUN/clinical_work \
  --sarek_outdir  /path/to/runs/MY_RUN/results \
  --run_name      MY_RUN \
  --outdir        /path/to/runs/MY_RUN/clinical_results \
  --dbnsfp        /path/to/refs/dbNSFP/dbNSFP4.1a.txt.gz \
  --snpsift_jar   /path/to/refs/snpEff/SnpSift.jar \
  --panel_bed     /path/to/refs/hg38/<kit>/<kit>_Regions.bed \
  --gcnv_ploidy_model   /path/to/refs/cnv/ploidy-model \
  --gcnv_pon_dir        /path/to/refs/cnv/PoN_CNV \
  --gcnv_interval_list  /path/to/refs/cnv/exomas.interval_list \
  --gcnv_fasta          /path/to/refs/hg38/hg38.fa \
  --gcnv_dict           /path/to/refs/hg38/hg38.dict \
  --gcnv_mane_txt       /path/to/refs/cnv/exons_mane.txt \
  --somatic_pon_hdf5       /path/to/refs/cnv/wes-do-gc.pon.hdf5 \
  --somatic_interval_list  /path/to/refs/cnv/preprocessed.interval_list \
  --common_snps_vcf        /path/to/refs/cnv/common_snps_preprocessed.vcf.gz \
  --somatic_fasta          /path/to/refs/cnv/main_chroms/hg38.main_chroms.fa \
  --somatic_dict           /path/to/refs/cnv/main_chroms/hg38.main_chroms.dict
```

#### Skip flags {#skip-flags}

| Flag | When to use |
|------|------------|
| `--skip_somatic_cnv true` | No ASCAT/CNVkit output |
| `--skip_hrd true` | No ASCAT segments |
| `--skip_germline_cnv true` | No normal BAMs / no gCNV PoN |
| `--skip_somatic_cnv_gatk true` | No tumor BAMs / no somatic PoN |
| `--skip_signatures true` | Too few somatic variants |
| `--skip_coverage true` | No mosdepth results |
| `--skip_pdf_report true` | HTML only |

---

## Requirements

| Tool | Version |
|------|---------|
| **Nextflow** | ≥ 23.04.0 |
| **Docker** | ≥ 20.10 |
| RAM | ≥ 32 GB (GATK gCNV peaks at ~56 GB) |
| Disk | ≥ 200 GB (reference databases + work dir) |

**Sarek prerequisite:** nf-core/sarek ≥ 3.3.0 with tools: `haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff`

---

## Parameters

### Required

| Parameter | Description |
|-----------|-------------|
| `--sarek_outdir` | Completed nf-core/sarek results directory |
| `--run_name` | Run identifier (used in filenames and report header) |
| `--outdir` | Where clinical results will be written |
| `--dbnsfp` | Path to `dbNSFP4.x.txt.gz` (bgzipped + tabix indexed) |
| `--snpsift_jar` | Path to `SnpSift.jar` |

### SNV & TMB

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--germline_panel` | `assets/gene_panels/germline_90_cancer_predisposition.txt` | 90-gene hereditary cancer panel |
| `--somatic_snv_panel` | `assets/gene_panels/somatic_snv_panel.txt` | Somatic SNV gene list |
| `--germline_min_dp` | `5` | Min depth for germline variants |
| `--panel_bed` | — | WES capture BED for TMB Mb calculation |
| `--tmb_min_dp` | `10` | Min tumor depth for TMB |
| `--tmb_min_af` | `0.05` | Min allele fraction for TMB |

### HRD

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hrd_organ` | `BREAST` | Organ model: `BREAST` or `OVARIAN` |
| `--hrd_threshold` | `0.45` | HRD probability threshold |

### Module toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_germline_cnv` | `false` | Skip GATK germline gCNV |
| `--skip_somatic_cnv` | `false` | Skip ASCAT + CNVkit interpretation |
| `--skip_somatic_cnv_gatk` | `false` | Skip GATK somatic CNV |
| `--skip_signatures` | `false` | Skip mutational signatures |
| `--skip_hrd` | `false` | Skip HRD |
| `--skip_pharmaco` | `false` | Skip pharmacogenomics |
| `--skip_coverage` | `false` | Skip gene coverage |
| `--skip_pdf_report` | `false` | Skip clinical PDF |

---

## Output

```
clinical_results/
├── germline/
│   ├── snv/extract/
│   │   ├── *.germline_snvs.tsv
│   │   ├── *.germline_all_filtered.tsv
│   │   └── *.germline_snvs_acmg_annotated.tsv
│   ├── pharmaco/            *.pharmaco.tsv
│   └── cnv_gatk/annotated/  *_segments_annotated.tsv
├── somatic/
│   ├── snv/extract/         *.somatic_snvs.tsv
│   ├── snv/civic_extract/   *.civic.exploded.tsv
│   ├── cnv/                 *_somatic_cnv_report.tsv
│   ├── cnv_gatk/annotated/  *_annotated.tsv
│   ├── tmb/                 *.tmb.tsv
│   ├── signatures/          SBS/ID/DBS matrices + 4 PNG plots
│   └── hrd/                 *_hrd_predictions.tsv
├── qc/gene_coverage/        *_gene_coverage.tsv
└── reporting/
    ├── *_molecular_report.html
    ├── *_molecular_report.txt
    └── *_clinical_report.pdf
```

---

## Gene Panels

Built-in panels in `assets/gene_panels/`:

- **`germline_90_cancer_predisposition.txt`** — 90 hereditary cancer genes (ACMG SF v3.2)
- **`somatic_38_oncogenes.txt`** — 38 actionable somatic genes for CNV reporting
- **`somatic_snv_panel.txt`** — Somatic SNV panel (oncogenes, TSGs, DNA repair)

Custom panel: `--germline_panel /path/to/my_panel.txt` (one gene symbol per line).

---

## Troubleshooting

| Problem | Solution |
|---------|---------|
| Session lock error | `rm -f <pipeline-dir>/.nextflow/cache/*/db/LOCK` |
| Task fails after code fix | Delete stale work dir: `rm -rf clinical_work/XX/XXXX…` then `-resume` |
| Out of memory | Lower `process_high` memory in `conf/base.config` |
| gCNV model version mismatch | Rebuild PoN with same GATK version (`broadinstitute/gatk:4.6.2.0`) |
| CIViC annotation fails | Requires internet on first run; passes through unannotated on failure |
| ACMG scores differ HTML vs PDF | Check `bin/final_report_acmg.py`, `bin/clinical_pdf_report.py`, `modules/local/snpsift_germline/main.nf` Step 4 are in sync |

---

## Citation

If you use `nf-sarek-clinical` in your research, please cite:

- **nf-core/sarek**: Garcia M, et al. (2020). *Genome Medicine*.
- **dbNSFP**: Liu X, et al. (2020). *Genome Medicine*.
- **CIViC/civicpy**: Griffith M, et al. (2017). *Nature Genetics*.
- **SigProfilerAssignment**: Bergstrom EN, et al. (2022). *Cell Genomics*.
- **HRProfiler**: Chopra N, et al. (2024). [github.com/FlorPio/HRProfiler-workflow](https://github.com/FlorPio/HRProfiler-workflow).
- **GATK CNV workflow**: [github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic).
- **ASCAT**: Van Loo P, et al. (2010). *PNAS*.

## License

MIT — see [LICENSE](LICENSE).
