# nf-sarek-clinical

**Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exome outputs.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

![Pipeline overview](metro_map_sarek_clinical.png)

---

## Reference files

All bioinformatics **tools** (GATK, Python, R packages, SnpSift, etc.) are bundled inside the Docker images — no manual tool installation required.

What you **do** need to provide are the large **reference databases and model files** listed below. Nextflow automatically bind-mounts them into the containers at runtime.

Define all paths once in `full_pipeline.params.yml` (template provided in the repo root).

---

### 1 · dbNSFP — variant functional annotation database

Used for germline and somatic SNV annotation (in silico predictors, REVEL, gnomAD AF, ClinVar).

```bash
# Download dbNSFP 4.x from the Liu Lab (requires registration):
# https://sites.google.com/site/jpopgen/dbNSFP
# Choose the hg38 version: dbNSFP4.x.txt.gz

mkdir -p /path/to/refs/dbNSFP
cd /path/to/refs/dbNSFP
# After download, index with tabix:
bgzip -@ 8 dbNSFP4.x.txt      # if downloaded as plain text
tabix -s 1 -b 2 -e 2 dbNSFP4.x.txt.gz
```

---

### 2 · SnpSift — VCF filtering and extraction

SnpSift is part of the **snpEff** package.

```bash
# Download snpEff (includes SnpSift):
# https://pcingola.github.io/SnpEff/  →  "Download" → snpEff_latest_core.zip

mkdir -p /path/to/refs/snpEff
cd /path/to/refs/snpEff
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
# SnpSift.jar will be at /path/to/refs/snpEff/snpEff/SnpSift.jar
```

---

### 3 · hg38 reference genome

Required by GATK germline and somatic CNV modules.

```bash
mkdir -p /path/to/refs/hg38
cd /path/to/refs/hg38

# Option A — from GATK Google Cloud bucket (fast):
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta        hg38.fa
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai    hg38.fa.fai
gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict         hg38.dict

# Option B — from Broad FTP:
# https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg38/v0

# Create main-chromosomes subset (required for somatic GATK CNV):
samtools view -H hg38.fa | grep "^@SQ" | awk '{print $2}' | grep -E "SN:chr[0-9XYM]+" \
  | sed 's/SN://' > main_chroms.txt
samtools faidx hg38.fa $(cat main_chroms.txt | tr '\n' ' ') > main_chroms/hg38.main_chroms.fa
samtools faidx main_chroms/hg38.main_chroms.fa
gatk CreateSequenceDictionary -R main_chroms/hg38.main_chroms.fa
```

---

### 4 · WES capture BED

Your sequencing kit-specific BED file. Used to calculate the capture size for TMB (mut/Mb).

- **Agilent SureSelect**: download from [Agilent eArray](https://earray.chem.agilent.com/suredesign/) (free registration)
- **Twist / KAPA / other**: contact the kit manufacturer for the target BED

Place it at `/path/to/refs/hg38/<kit>/<kit>_Regions.bed`.

---

### 5 · MANE exon coordinates

Used to annotate germline CNV segments with gene/exon names.

```bash
# Download MANE summary from NCBI:
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.4.summary.txt.gz
gunzip MANE.GRCh38.v1.4.summary.txt.gz

# Convert to the 5-column format expected by the pipeline
# (chr, start, end, gene_symbol, exon, transcript_id):
# A conversion script is provided in assets/scripts/mane_to_exons.py (TODO: add to repo)
# The expected format is tab-separated, no header:
# chr1  69091  70008  OR4F5  1  NM_001005484
```

---

### 6 · GATK somatic CNV support files

The complete workflow for building these files (interval lists, common SNPs VCF, PoN) is available at:
**[github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)**

That repository also contains the two R annotation scripts required by this pipeline:
- `exon_annotation.R` — germline CNV exon annotation
- `annotate_cnvs_complete_fixed.R` — somatic CNV gene annotation

```bash
git clone https://github.com/FlorPio/gatk-cnv-somatic /path/to/gatk-cnv-somatic
# Follow the instructions in that repo to generate:
#   preprocessed.interval_list
#   common_snps_preprocessed.vcf.gz
# and copy the R scripts to your refs/cnv/ directory.
cp /path/to/gatk-cnv-somatic/exon_annotation.R          /path/to/refs/cnv/
cp /path/to/gatk-cnv-somatic/annotate_cnvs_complete_fixed.R /path/to/refs/cnv/
```

---

### 7 · GATK CNV Panel of Normals (build once from your cohort)

PoN files are **institution-specific** — they must be built from your own cohort of normal samples (typically ≥ 30) using the same sequencing kit and GATK version (`broadinstitute/gatk:4.6.2.0`).

Follow the step-by-step instructions in **[github.com/FlorPio/gatk-cnv-somatic](https://github.com/FlorPio/gatk-cnv-somatic)** to generate:

| File | GATK tool |
|------|-----------|
| `ploidy-model/` | `DetermineGermlineContigPloidy` (COHORT mode) |
| `PoN_CNV/` | `GermlineCNVCaller` (COHORT mode) |
| `wes-do-gc.pon.hdf5` | `CreateReadCountPanelOfNormals` |
| `exomas.interval_list` | `FilterIntervals` on germline read counts |

> **Note:** If you rebuild the PoN with a different GATK version you must rerun all CNV calling for existing samples.

---

### Suggested directory layout

```
/path/to/refs/
├── dbNSFP/
│   ├── dbNSFP4.1a.txt.gz
│   └── dbNSFP4.1a.txt.gz.tbi
├── snpEff/
│   └── SnpSift.jar
├── hg38/
│   ├── hg38.fa  (+ .fai, .dict, BWA indexes)
│   └── <kit>/
│       └── <kit>_Regions.bed
└── cnv/
    ├── exons_mane.txt
    ├── exomas.interval_list            ← germline FilterIntervals output
    ├── preprocessed.interval_list      ← somatic PreprocessIntervals output
    ├── common_snps_preprocessed.vcf.gz (+ .tbi)
    ├── exon_annotation.R               ← from FlorPio/gatk-cnv-somatic
    ├── annotate_cnvs_complete_fixed.R  ← from FlorPio/gatk-cnv-somatic
    ├── ploidy-model/                   ← built from cohort
    ├── PoN_CNV/                        ← built from cohort
    ├── wes-do-gc.pon.hdf5              ← built from cohort
    └── main_chroms/
        ├── hg38.main_chroms.fa  (+ .fai)
        └── hg38.main_chroms.dict
```

### Configure your paths

Edit `full_pipeline.params.yml` (copy from the template in the repo root) and set each path to your local files. This file is read automatically by `run_full.sh`.

---

## Running the pipeline

### Option A — Full pipeline from FASTQ (recommended)

`run_full.sh` orchestrates nf-core/sarek (alignment + variant calling) followed by nf-sarek-clinical (clinical report).

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
# Resume after interruption
./run_full.sh --run <RUN_NAME> --samplesheet /path/to/samplesheet.csv --resume

# Clinical report only (sarek already completed)
./run_full.sh --run <RUN_NAME> --skip-sarek

# Resume clinical report only
./run_full.sh --run <RUN_NAME> --skip-sarek --resume
```

#### Output layout

```
/path/to/runs/<RUN_NAME>/
├── samplesheet.csv          ← input
├── results/                 ← nf-core/sarek outputs (BAMs, VCFs, CNVs, QC)
├── clinical_work/           ← nf-sarek-clinical work dir
└── clinical_results/        ← HTML + PDF clinical reports  ← final output
```

---

### Option B — Clinical report only

Use when nf-core/sarek was run separately. Requires Nextflow ≥ 23.04.0 and Docker ≥ 20.10.

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
  --gcnv_exon_r         /path/to/refs/cnv/exon_annotation.R \
  --somatic_pon_hdf5       /path/to/refs/cnv/wes-do-gc.pon.hdf5 \
  --somatic_interval_list  /path/to/refs/cnv/preprocessed.interval_list \
  --common_snps_vcf        /path/to/refs/cnv/common_snps_preprocessed.vcf.gz \
  --somatic_fasta          /path/to/refs/cnv/main_chroms/hg38.main_chroms.fa \
  --somatic_dict           /path/to/refs/cnv/main_chroms/hg38.main_chroms.dict \
  --somatic_annotate_r     /path/to/refs/cnv/annotate_cnvs_complete_fixed.R
```

#### Skip flags (when data is missing)

| Flag | When to use |
|------|------------|
| `--skip_somatic_cnv true` | No ASCAT/CNVkit output |
| `--skip_hrd true` | No ASCAT segments |
| `--skip_germline_cnv true` | No normal BAMs |
| `--skip_somatic_cnv_gatk true` | No tumor BAMs |
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
│   │   ├── *.germline_snvs.tsv                   — variants passing all filters (P/LP/VUS in report)
│   │   ├── *.germline_all_filtered.tsv            — all rejected variants with FILTER_REASON column
│   │   └── *.germline_snvs_acmg_annotated.tsv     — all candidates with ACMG_class, criteria, insil scores
│   ├── pharmaco/            *.pharmaco.tsv
│   └── cnv_gatk/annotated/  *_segments_annotated.tsv
├── somatic/
│   ├── snv/extract/         *.somatic_snvs.tsv
│   ├── snv/civic_extract/   *.civic.exploded.tsv
│   ├── cnv/                 *_somatic_cnv_report.tsv
│   ├── cnv_gatk/annotated/  *_annotated.tsv
│   ├── tmb/                 *.tmb.tsv
│   ├── signatures/
│   │   ├── matrix_profiles/output/   SBS/, ID/, DBS/ observed matrices
│   │   └── assignment/               COSMIC v3.4 activities + 4 PNG plots
│   └── hrd/                 *_hrd_predictions.tsv
├── qc/gene_coverage/        *_gene_coverage.tsv
└── reporting/
    ├── *_molecular_report.html
    ├── *_molecular_report.txt
    └── *_clinical_report.pdf
```

### Report section order (PDF)

| # | Section | Key filter |
|---|---------|-----------|
| 0 | **Major Clinical Findings** | HRD/MMR/CNV triggers |
| 1 | **Germline–Somatic Correlation** | P/LP only; second hit or biomarker match |
| 2 | **Germline SNVs** | ACMG P/LP/VUS |
| 3 | **Somatic SNVs** | All from filtered TSV |
| 4 | **Somatic CNVs (ASCAT/CNVkit)** | 38-gene panel |
| 5 | **TMB · Pharmacogenomics · Gene Coverage · HRD** | — |
| 6 | **Mutational Signatures** | SBS96 + ID83 spectra + etiology plots |
| 7 | **Germline/Somatic CNVs (GATK)** | All annotated segments |

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
| Task fails after code fix | Delete the stale work dir: `rm -rf clinical_work/XX/XXXX…` then `-resume` |
| Out of memory | Lower `process_high` memory in `conf/base.config` |
| gCNV model version mismatch | Rebuild PoN with same GATK version (`broadinstitute/gatk:4.6.2.0`) |
| CIViC annotation fails | Requires internet on first run; step passes through unannotated on failure |
| ACMG scores differ between HTML and PDF | Check all three ACMG implementations are in sync: `bin/final_report_acmg.py`, `bin/clinical_pdf_report.py`, `modules/local/snpsift_germline/main.nf` Step 4 |

---

## Citation

If you use `nf-sarek-clinical` in your research, please cite:

- **nf-core/sarek**: Garcia M, et al. (2020). *Genome Medicine*.
- **dbNSFP**: Liu X, et al. (2020). *Genome Medicine*.
- **CIViC/civicpy**: Griffith M, et al. (2017). *Nature Genetics*.
- **SigProfilerAssignment**: Bergstrom EN, et al. (2022). *Cell Genomics*.
- **HRProfiler**: Chopra N, et al. (2024). [github.com/FlorPio/HRProfiler-workflow](https://github.com/FlorPio/HRProfiler-workflow).
- **ASCAT**: Van Loo P, et al. (2010). *PNAS*.

## License

MIT — see [LICENSE](LICENSE).
