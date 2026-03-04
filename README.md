# nf-sarek-clinical

**Clinical interpretation postprocessing pipeline for nf-core/sarek WES cancer exome outputs.**

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)

![Pipeline overview](metro_map_sarek_clinical.png)

---

## What this pipeline does

`nf-sarek-clinical` takes the output of [nf-core/sarek](https://nf-co.re/sarek) (tumour–normal WES) and produces a clinical-grade molecular report in HTML and PDF format. It integrates germline and somatic variant analysis, copy-number calling, mutational signatures, pharmacogenomics, and HRD scoring into a single automated workflow. Also can you run in a full mode from fastq to report.

**Final deliverables per patient:**

- `*_molecular_report.html` — interactive report with all findings
- `*_clinical_report.pdf` — printable PDF with a Major Clinical Findings summary page
- Structured TSV files for each analysis module (variant tables, CNV tables, TMB, HRD)

---

## Report sections

| # | Section | What it analyses |
|---|---------|-----------------|
| 0 | **Major Clinical Findings** | Auto-populated summary page: HRD, MMR deficiency, high TMB, relevant CNVs, and pathogenic germline SNV and CNV. |
| 1 | **Germline–Somatic Correlation** | Two-hit events (P/LP germline + somatic second hit), biallelic somatic, signature correlation |
| 2 | **Germline SNVs** | 90-gene hereditary cancer panel; ACMG P/LP/VUS classification with 9-tool in silico scoring, REVEL, ClinVar, gnomAD |
| 3 | **Pharmacogenomics** | DPYD, UGT1A1 — drug-metabolism variants from germline VCF |
| 4 | **Germline CNVs (GATK gCNV)** | Copy-number calls on normal BAM using a pre-built cohort PoN; MANE gene/exon annotation |
| 5 | **Somatic SNVs** | Mutect2 VCF filtered by somatic panel; AMP/ASCO/CAP tiering with CIViC evidence |
| 6 | **Somatic CNVs (ASCAT/CNVkit)** | HOMDEL / HETDEL / AMP / LOH on 38-gene actionable panel from sarek CNV outputs |
| 7 | **Somatic CNVs (GATK)** | Allele-specific paired CNV (tumour + normal BAMs); MANE annotation; diagnostic plots |
| 8 | **TMB** | Tumour Mutational Burden in mut/Mb from PASS Mutect2 variants over WES capture size |
| 9 | **Mutational Signatures** | SigProfilerAssignment COSMIC v3.4 SBS96/DBS/ID83 decomposition + 4 plots |
| 10 | **HRD** | HRProfiler logistic model (SBS + indel + CNV48 features) → probability score + classification |
| 11 | **Gene Coverage** | Per-gene mean/min/max depth from mosdepth (normal + tumour) |

---

## Dependencies

### Software (install once on the host)

| Tool | Minimum version | Install |
|------|----------------|---------|
| **Nextflow** | 23.04.0 | `curl -s https://get.nextflow.io \| bash` |
| **Docker** | 20.10 | [docs.docker.com/get-docker](https://docs.docker.com/get-docker/) |
| **Java** | 11 (for Nextflow) | `sudo apt install default-jre` |

> All bioinformatics tools (GATK, Python, R, SigProfiler, SnpSift, etc.) run inside Docker containers — **no manual tool installation is needed**.

### nf-core/sarek (pulled automatically)

nf-core/sarek ≥ 3.3.0 is pulled automatically by Nextflow in Mode A. It must have been run with at minimum:

```
--tools haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff
```

### Hardware

| Resource | Minimum | Notes |
|----------|---------|-------|
| RAM | 32 GB | GATK gCNV peaks at ~56 GB |
| Disk | 200 GB | Reference databases + work directory |
| CPU | 8 cores | More = faster parallel tasks |

---

## First-time setup

### Step 1 — Clone the repository

```bash
git clone https://github.com/GENUIT-LBMU/nf-sarek-clinical.git /path/to/nf-sarek-clinical
cd /path/to/nf-sarek-clinical
```

> For the private GitHub repo, configure Git credentials first:
> ```bash
> git config --global credential.helper store
> # or use SSH key added to the GENUIT-LBMU org
> ```

### Step 2 — Configure reference file paths

Copy and edit the params template:

```bash
cp full_pipeline.params.yml my_lab_params.yml
```

Edit `my_lab_params.yml` with your local reference file paths (see [Reference files](#reference-files) below for where to obtain each file):

```yaml
# SNV annotation
dbnsfp:          /data/refs/dbNSFP/dbNSFP4.1a.txt.gz
snpsift_jar:     /data/refs/snpEff/SnpSift.jar

# WES capture BED (for TMB calculation)
panel_bed:       /data/refs/hg38/agilent/S07604514_Regions.bed
gene_coords_bed: /data/refs/hg38/agilent/S07604514_Padded.withGenes.bed

# Germline GATK CNV
gcnv_ploidy_model:  /data/refs/cnv/ploidy-model
gcnv_pon_dir:       /data/refs/cnv/PoN_CNV
gcnv_interval_list: /data/refs/cnv/exomas.interval_list
gcnv_fasta:         /data/refs/hg38/hg38.fa
gcnv_dict:          /data/refs/hg38/hg38.dict
gcnv_mane_txt:      /data/refs/cnv/exons_mane.txt

# Somatic GATK CNV
somatic_pon_hdf5:        /data/refs/cnv/wes-do-gc.pon.hdf5
somatic_interval_list:   /data/refs/cnv/preprocessed.interval_list
common_snps_vcf:         /data/refs/cnv/common_snps_preprocessed.vcf.gz
somatic_fasta:           /data/refs/cnv/main_chroms/hg38.main_chroms.fa
somatic_dict:            /data/refs/cnv/main_chroms/hg38.main_chroms.dict
```

> The R annotation scripts (`exon_annotation.R`, `annotate_cnvs_complete_fixed.R`) are **bundled in `bin/`** — no path needed.

### Step 3 — Configure sarek parameters (Mode A only)

Edit `conf/sarek.config` with your lab's nf-core/sarek reference paths. These are constant across all runs:

```nextflow
params {
    genome       = 'GATK.GRCh38'
    wes          = true
    save_output_as_bam = true
    intervals    = '/data/refs/hg38/agilent/S07604514_Regions_Annotated.bed'
    ascat_genome = 'hg38'
    ascat_alleles = '/data/refs/hg38/ascat/ascat_alleles.zip'
    ascat_loci    = '/data/refs/hg38/ascat/ascat_loci.zip'
    ascat_loci_gc = '/data/refs/hg38/ascat/ascat_loci_gc.zip'
    ascat_loci_rt = '/data/refs/hg38/ascat/ascat_loci_rt.zip'
}
```

Also edit line 23 of `run_full.sh` to set your runs base directory:

```bash
RUNS_BASE="/your/runs/directory"   # default: /home/user/runs
```

---

## Mode A — Full pipeline (FASTQ → clinical report)

Use this mode when starting from raw sequencing reads. `run_full.sh` orchestrates both phases sequentially: **Phase 1** runs nf-core/sarek (alignment, variant calling, QC) and **Phase 2** runs nf-sarek-clinical (clinical interpretation and report).

### Step 1 — Create the run directory and samplesheet

```bash
mkdir -p /path/to/runs/20260109_1229
```

Create `samplesheet.csv` following nf-core/sarek format. **Sample naming is critical**: use `<patient>_N##` for the normal and `<patient>_T##` for the tumour.

```csv
patient,sex,status,sample,lane,fastq_1,fastq_2
1229,XY,0,1229_N07,2,/path/to/reads/1229_N07_L2_1.fq.gz,/path/to/reads/1229_N07_L2_2.fq.gz
1229,XY,1,1229_T08,1,/path/to/reads/1229_T08_L1_1.fq.gz,/path/to/reads/1229_T08_L1_2.fq.gz
1229,XY,1,1229_T08,2,/path/to/reads/1229_T08_L2_1.fq.gz,/path/to/reads/1229_T08_L2_2.fq.gz
```

Fields:
- `patient` — unique patient identifier (used in output filenames)
- `sex` — `XX` (female) or `XY` (male); used by ASCAT for sex-chromosome copy-number
- `status` — `0` = normal, `1` = tumour
- `sample` — sample name; must follow `<patient>_N##` / `<patient>_T##`
- `lane` — sequencing lane number (integer); multiple lanes per sample are merged automatically
- `fastq_1` / `fastq_2` — absolute paths to gzipped FASTQ files

### Step 2 — Run

```bash
cd /path/to/nf-sarek-clinical

./run_full.sh \
    --run        20260109_1229 \
    --samplesheet /path/to/runs/20260109_1229/samplesheet.csv
```

Nextflow will create the output directories automatically:

```
/path/to/runs/20260109_1229/
├── samplesheet.csv
├── run.config             ← auto-generated per-run Nextflow config
├── work/                  ← sarek intermediate files
├── results/               ← nf-core/sarek outputs (BAMs, VCFs, QC)
├── clinical_work/         ← nf-sarek-clinical intermediate files
└── clinical_results/      ← final HTML + PDF reports  ← HERE
```

### Step 3 — Resume after interruption

If the run stops (server reboot, OOM, etc.), resume from where it left off:

```bash
./run_full.sh \
    --run        20260109_1229 \
    --samplesheet /path/to/runs/20260109_1229/samplesheet.csv \
    --resume
```

---

## Mode B — Post-sarek (sarek already completed → clinical report)

Use this mode when nf-core/sarek has already been run and you only need the clinical interpretation.

### Required sarek output structure

The pipeline expects the following directories inside `<sarek_outdir>/`:

```
results/
├── variant_calling/
│   ├── mutect2/<pair>/         ← somatic SNV VCFs
│   ├── haplotypecaller/<normal>/ ← germline VCF
│   ├── ascat/<pair>/           ← ASCAT segments + metrics
│   ├── cnvkit/<pair>/          ← CNVkit .cns / .cnr files
│   └── msisensorpro/<pair>/    ← MSI score file
├── preprocessing/recalibrated/ ← BAM + BAI files (for GATK CNV)
└── reports/mosdepth/           ← coverage summaries
```

> BAMs are required only for GATK CNV modules. If BAMs are absent, use `--skip_germline_cnv true --skip_somatic_cnv_gatk true`.

### Option B.1 — Using run_full.sh (local clone)

The simplest approach when working from the cloned repo directory:

```bash
cd /path/to/nf-sarek-clinical

./run_full.sh \
    --run        20260109_1229 \
    --skip-sarek
```

This reads reference paths from `full_pipeline.params.yml` (in the pipeline dir) and uses sarek results from `/path/to/runs/20260109_1229/results/`.

Resume if interrupted:

```bash
./run_full.sh --run 20260109_1229 --skip-sarek --resume
```

### Option B.2 — Pulling directly from GitHub (recommended for external users)

No local clone needed. Configure Nextflow SCM credentials once:

```bash
mkdir -p ~/.nextflow
cat > ~/.nextflow/scm <<'EOF'
providers {
    github {
        user     = 'your-github-username'
        password = 'your-personal-access-token'
    }
}
EOF
```

Then run:

```bash
nextflow run GENUIT-LBMU/nf-sarek-clinical \
    -profile docker \
    -work-dir /path/to/runs/20260109_1229/clinical_work \
    -params-file /path/to/my_lab_params.yml \
    --sarek_outdir /path/to/runs/20260109_1229/results \
    --run_name     20260109_1229 \
    --outdir       /path/to/runs/20260109_1229/clinical_results
```

Resume if interrupted:

```bash
nextflow run GENUIT-LBMU/nf-sarek-clinical \
    -profile docker \
    -work-dir /path/to/runs/20260109_1229/clinical_work \
    -params-file /path/to/my_lab_params.yml \
    -resume \
    --sarek_outdir /path/to/runs/20260109_1229/results \
    --run_name     20260109_1229 \
    --outdir       /path/to/runs/20260109_1229/clinical_results
```

---

## Skip flags

Use these when a data type is missing or a module is not needed:

| Flag | Default | When to use |
|------|---------|-------------|
| `--skip_germline_cnv true` | `false` | No normal BAMs, or no gCNV PoN built yet |
| `--skip_somatic_cnv true` | `false` | No ASCAT / CNVkit output from sarek |
| `--skip_somatic_cnv_gatk true` | `false` | No tumour BAMs, or no somatic PoN |
| `--skip_hrd true` | `false` | No ASCAT segments (required by HRProfiler) |
| `--skip_signatures true` | `false` | Too few somatic variants (< ~50) |
| `--skip_coverage true` | `false` | No mosdepth output from sarek |
| `--skip_pharmaco true` | `false` | Pharmacogenomics not needed |
| `--skip_pdf_report true` | `false` | Generate HTML only (faster) |

Example — skip both GATK CNV modules and HRD:

```bash
./run_full.sh --run 20260109_1229 --skip-sarek \
    --skip_germline_cnv true \
    --skip_somatic_cnv_gatk true \
    --skip_hrd true
```

Or with `nextflow run`:

```bash
nextflow run GENUIT-LBMU/nf-sarek-clinical \
    ... \
    --skip_germline_cnv true \
    --skip_somatic_cnv_gatk true \
    --skip_hrd true
```

---

## Reference files

All Docker images contain the bioinformatics tools. The files below must be downloaded once and stored locally. Configure their paths in `my_lab_params.yml`.

### dbNSFP — variant functional annotation

Download from the **Liu Lab** (free registration required):
[sites.google.com/site/jpopgen/dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) → choose hg38 version.

```bash
mkdir -p /data/refs/dbNSFP && cd /data/refs/dbNSFP
# After download:
tabix -s 1 -b 2 -e 2 dbNSFP4.x.txt.gz
```

### SnpSift — VCF filtering and field extraction

Part of the **snpEff** package:

```bash
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
unzip snpEff_latest_core.zip -d /data/refs/snpEff
# Jar → /data/refs/snpEff/snpEff/SnpSift.jar
```

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

Download from your kit vendor. Used to calculate the effective Mb for TMB.

- **Agilent SureSelect**: [earray.chem.agilent.com/suredesign](https://earray.chem.agilent.com/suredesign/) (free registration)
- **Twist / KAPA / Roche**: request from the manufacturer

### MANE exon coordinates

```bash
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/MANE.GRCh38.v1.4.summary.txt.gz
gunzip MANE.GRCh38.v1.4.summary.txt.gz
# Convert to 6-column TSV (chr, start, end, gene_symbol, exon, transcript_id)
# Conversion script available in FlorPio/gatk-cnv-somatic
```

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

> PoN files are **institution-specific**: built from your own cohort of ≥ 30 normal samples with the same kit and the same GATK version (`broadinstitute/gatk:4.6.2.0`). Rebuild if GATK version changes.

> The R annotation scripts (`exon_annotation.R`, `annotate_cnvs_complete_fixed.R`) from that repo are already **bundled in `bin/`** — no manual copy needed.

### ASCAT resources (Mode A only)

Required for ASCAT allele-specific copy-number in nf-core/sarek. Download from the **ASCAT GitHub**:
[github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES](https://github.com/VanLoo-lab/ascat/tree/master/ReferenceFiles/WES)

Download the four zip files for hg38 (`alleles`, `loci`, `loci_gc`, `loci_rt`) and set their paths in `conf/sarek.config`.

### Suggested directory layout

```
/data/refs/
├── dbNSFP/
│   ├── dbNSFP4.1a.txt.gz
│   └── dbNSFP4.1a.txt.gz.tbi
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

## Parameters reference

### Required

| Parameter | Description |
|-----------|-------------|
| `--sarek_outdir` | Path to completed nf-core/sarek results directory |
| `--run_name` | Run identifier used in output filenames and report header |
| `--outdir` | Where clinical results will be written |
| `--dbnsfp` | `dbNSFP4.x.txt.gz` (bgzipped + tabix indexed) |
| `--snpsift_jar` | `SnpSift.jar` from the snpEff package |

### SNV & TMB

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--germline_panel` | `assets/gene_panels/germline_90_cancer_predisposition.txt` | Gene list for germline SNV filter (one symbol per line) |
| `--somatic_snv_panel` | `assets/gene_panels/somatic_snv_panel.txt` | Gene list for somatic SNV filter |
| `--germline_min_dp` | `5` | Minimum read depth for germline variants |
| `--panel_bed` | — | WES capture BED for TMB calculation (mut/Mb) |
| `--tmb_min_dp` | `10` | Minimum tumour depth for TMB variants |
| `--tmb_min_af` | `0.05` | Minimum allele fraction for TMB variants |

### HRD

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--hrd_organ` | `BREAST` | Organ model: `BREAST` or `OVARIAN` |
| `--hrd_threshold` | `0.45` | Probability threshold for HRD+ classification |

### Module toggles

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--skip_germline_cnv` | `false` | Skip GATK germline gCNV |
| `--skip_somatic_cnv` | `false` | Skip ASCAT + CNVkit interpretation |
| `--skip_somatic_cnv_gatk` | `false` | Skip GATK somatic CNV |
| `--skip_signatures` | `false` | Skip mutational signatures |
| `--skip_hrd` | `false` | Skip HRD analysis |
| `--skip_pharmaco` | `false` | Skip pharmacogenomics |
| `--skip_coverage` | `false` | Skip gene coverage |
| `--skip_pdf_report` | `false` | Skip clinical PDF (HTML only) |

---

## Output

```
clinical_results/
├── germline/
│   ├── snv/extract/
│   │   ├── *.germline_snvs.tsv                  ← P/LP/VUS variants (in report)
│   │   ├── *.germline_all_filtered.tsv           ← all rejected with FILTER_REASON
│   │   └── *.germline_snvs_acmg_annotated.tsv   ← all candidates + ACMG scores
│   ├── pharmaco/
│   │   └── *.pharmaco.tsv
│   └── cnv_gatk/annotated/
│       └── *_segments_annotated.tsv
├── somatic/
│   ├── snv/extract/         *.somatic_snvs.tsv
│   ├── snv/civic_extract/   *.civic.exploded.tsv
│   ├── cnv/                 *_somatic_cnv_report.tsv
│   ├── cnv_gatk/annotated/  *_annotated.tsv
│   ├── tmb/                 *.tmb.tsv
│   ├── signatures/
│   │   ├── matrix_profiles/output/   SBS96/, DBS/, ID/ matrices
│   │   └── assignment/               COSMIC v3.4 activities + 4 PNG plots
│   └── hrd/                 *_hrd_predictions.tsv
├── qc/gene_coverage/        *_gene_coverage.tsv
└── reporting/
    ├── *_molecular_report.html
    ├── *_molecular_report.txt
    └── *_clinical_report.pdf
```

---

## Gene panels

Built-in panels bundled in `assets/gene_panels/`:

- **`germline_90_cancer_predisposition.txt`** — 90 hereditary cancer genes (ACMG SF v3.2)
- **`somatic_38_oncogenes.txt`** — 38 actionable genes for somatic CNV reporting
- **`somatic_snv_panel.txt`** — Somatic SNV panel (oncogenes, TSGs, DNA repair)

To use a custom panel: `--germline_panel /path/to/my_panel.txt` (one HGNC gene symbol per line).

---

## Troubleshooting

| Problem | Solution |
|---------|---------|
| `Session [xxx] lock file already acquired` | `rm -f <pipeline-dir>/.nextflow/cache/*/db/LOCK` |
| Task fails after a code fix | Delete stale work dir: `rm -rf clinical_work/XX/XXXX…` then re-run with `-resume` |
| Out of memory (GATK gCNV) | Lower `process_high` memory in `conf/base.config` (default 56 GB) |
| gCNV: `model version mismatch` | Rebuild PoN with same GATK version (`broadinstitute/gatk:4.6.2.0`) |
| CIViC annotation fails | Requires internet on first run; passes through unannotated on network failure |
| ACMG classification differs HTML vs PDF | Check that `bin/final_report_acmg.py`, `bin/clinical_pdf_report.py`, and `modules/local/snpsift_germline/main.nf` Step 4 are in sync |
| `Join mismatch` error in sarek (Mutect2) | `--tools` must be passed on the CLI in `run_full.sh`, not only in `conf/sarek.config` (already handled) |

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
