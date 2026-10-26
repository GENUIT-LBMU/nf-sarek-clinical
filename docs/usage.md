# Usage

## Prerequisites

### 1. Run nf-core/sarek first

`nf-sarek-clinical` is a **downstream** pipeline. You must first complete an nf-core/sarek run with the following tools:

```bash
nextflow run nf-core/sarek \
  -profile docker \
  --input samplesheet.csv \
  --outdir /data/runs/MY_RUN/results \
  --genome GATK.GRCh38 \
  --tools haplotypecaller,mutect2,cnvkit,ascat,msisensorpro,snpeff \
  --save_output_as_bam \
  --wes \
  --intervals /data/hg38/WES_panel.bed
```

The `--save_output_as_bam` flag is required for the GATK CNV modules (if enabled).

### 2. Prepare databases

| Database | Download | Index |
|----------|----------|-------|
| Local gnomAD v4.1, ClinVar, REVEL v1.3 (germline SNVs) | Built with `scripts/build_gnomad_panel.sh`, `scripts/build_clinvar_panel.sh`, `scripts/build_revel_panel.py` (see README) | written by the scripts |
| SnpSift.jar | Included in [snpEff package](https://pcingola.github.io/SnpEff/) | — |
| WES BED | From capture kit vendor | optional: `bgzip && tabix` |

## Running the pipeline

### Minimal run (Docker)

```bash
nextflow run nf-sarek-clinical/main.nf \
  -profile docker \
  --sarek_outdir /data/runs/MY_RUN/results \
  --run_name     MY_RUN_2024 \
  --snpsift_jar  /opt/snpEff/SnpSift.jar \
  --gnomad_tsv   /data/refs/gnomad/gnomad_v4.1_panel.tsv.gz \
  --clinvar_tsv  /data/refs/clinvar/clinvar_<DATE>_panel.tsv.gz \
  --revel_tsv    /data/refs/revel/revel_v1.3_panel.tsv.gz \
  --outdir       /results/MY_RUN_2024_clinical
```

### With all features

```bash
nextflow run nf-sarek-clinical/main.nf \
  -profile docker \
  --sarek_outdir  /data/runs/MY_RUN/results \
  --run_name      MY_RUN_2024 \
  --snpsift_jar   /opt/snpEff/SnpSift.jar \
  --gnomad_tsv    /data/refs/gnomad/gnomad_v4.1_panel.tsv.gz \
  --clinvar_tsv   /data/refs/clinvar/clinvar_<DATE>_panel.tsv.gz \
  --revel_tsv     /data/refs/revel/revel_v1.3_panel.tsv.gz \
  --panel_bed     /data/hg38/S07604514_Regions.bed \
  --gene_coords_bed /data/hg38/S07604514_Regions_Annotated.bed \
  --outdir        /results/MY_RUN_2024_clinical \
  --genome_build  GRCh38
```

### Skip heavy modules for faster runs

```bash
# Skip mutational signatures (requires SigProfiler conda env)
--skip_signatures true

# Skip CIViC annotation (requires internet access)
# CIViC is part of the somatic SNV module; to skip:
--skip_somatic_snv true
```

### Resuming a run

```bash
nextflow run nf-sarek-clinical/main.nf -resume [other params...]
```

## Custom gene panels

Replace the built-in panels with your own (one gene symbol per line, `#` for comments):

```bash
--germline_panel   /path/to/my_germline_genes.txt   # for SNV filter + LOH
--somatic_cnv_panel /path/to/my_somatic_cnv.txt     # for CNVkit HOMDEL/HETDEL/AMP
--somatic_snv_panel /path/to/my_somatic_snv.txt     # for Mutect2 SNV filter
```

## Samplesheet naming convention

This pipeline auto-discovers samples from the Sarek output directory. It expects the standard Sarek naming conventions:

- Normal samples: `{PATIENT_ID}_N{XX}` (e.g. `682_N01`)
- Tumor samples:  `{PATIENT_ID}_T{XX}` (e.g. `682_T01`)
- Pairs: `{TUMOR}_vs_{NORMAL}` (e.g. `682_T01_vs_682_N01`)

This matches the samplesheet format used by `launch_nt_V6.sh`.

## Output files

See [output.md](output.md) for a full description of all output files.
