# Output Description

## `results/germline/snv/`

| File | Description |
|------|-------------|
| `filtered/*.filtered.vcf.gz` | VCF filtered to: 90-gene cancer predisposition panel + HIGH/MODERATE impact + DP ≥ 5 |
| `extract/*.germline_snvs.tsv` | Tab-delimited table ready for Excel/R, with all annotation fields |

**Key TSV columns:** CHROM, POS, REF, ALT, ANN[*].GENE, ANN[*].HGVS_C / HGVS_P (MANE Select; `MANE_PLUS_CLINICAL`), GEN[0].GT, GEN[0].AD, ZYGOSITY, gnomAD_AF / gnomAD_FAF95 / gnomAD_AF_source (local gnomAD v4.1), ClinVar_sig / ClinVar_stars (local ClinVar), REVEL_score (local REVEL v1.3), SpliceAI_DS, ClinVar_PS1 / ClinVar_PM5, ACMG_class / ACMG_criteria / ACMG_spec (`*_acmg_annotated.tsv`)

## `results/germline/pharmaco/`

| File | Description |
|------|-------------|
| `*.pharmaco.tsv` | Pharmacogenomics variants (DPYD *2A/*13/HapB3, UGT1A1 *28). Relevant for fluoropyrimidine (5-FU, capecitabine) and irinotecan toxicity risk assessment. |

## `results/somatic/snv/`

| File | Description |
|------|-------------|
| `filtered/*.filtered.vcf.gz` | PASS + somatic panel (60+ genes) + HIGH/MODERATE + tumor DP ≥ 10 (MANE annotation and local gnomAD fields when configured) |
| `extract/*.somatic_snvs.tsv` | Somatic SNV table: MANE HGVS (+ `MANE_PLUS_CLINICAL`), VAF, callers, local gnomAD (`GNOMAD_*`, `gnomAD_AF`, `gnomAD_AF_source`, `gnomAD_FAF95`, `gnomAD_NHOMALT`, `gnomAD_POP_AF`, `POP_COMMON` = 1 for AMP Tier IV). |
| `extract/*.second_hit.somatic_snvs.tsv` | Second-hit track: somatic SNVs in germline second-hit genes outside the core (correlation only) |
| `extract/*.tso500.somatic_snvs.tsv` | Pan-tumor track: somatic SNVs in TSO500 genes outside the core (same columns) |
| `filtered/*.{second_hit,tso500}.filtered.vcf.gz`, `civic/*.tso500.civic.vcf.gz`, `civic_extract/*.tso500.civic.exploded.tsv` | Same steps for the tracks |
| `civic/*.civic.vcf.gz` | VCF annotated with CIViC clinical evidence |
| `civic_extract/*.civic.exploded.tsv` | CIViC evidence exploded to one row per evidence item |

**CIViC TSV fields include:** CIViC Variant Name, Molecular Profile, Evidence Level, Significance, Disease, Therapies, AMP Category, FDA Companion Test.

## `results/somatic/cnv/ascat_cnvkit/`

| File | Description |
|------|-------------|
| `*_somatic_cnv_report.tsv` | Master table of all somatic CNV events |
| `*_somatic_cnv_report.txt` | Human-readable clinical summary |

**TSV columns:** sample, event_type, gene, source, chrom, start, end, cn, log2, minor, loh_type, purity, ploidy, WGD, GI, LOH, goodness_of_fit.

**Event types:**
- `HOMDEL` (CN=0, CNVkit) — Homozygous deletion
- `HETDEL` (CN=1, CNVkit) — Heterozygous deletion
- `AMP` (CN≥5, CNVkit) — Amplification
- `2nd_hit_LOH` (ASCAT) — Loss of heterozygosity in germline cancer genes

**Gene panels:**
- CNVkit events: filtered to 38-gene somatic panel (KRAS, BRAF, EGFR, ERBB2, TP53, BRCA1/2, MLH1, MSH2, etc.)
- ASCAT LOH: filtered to 90-gene germline cancer predisposition panel

## `results/somatic/tmb/`

| File | Description |
|------|-------------|
| `*.tmb.tsv` | Tumor Mutational Burden per sample |

**Columns:** sample, n_variants, panel_bases, panel_mb, TMB_mut_per_mb

**Filtering criteria:** FILTER=PASS + HIGH/MODERATE impact + tumor DP≥10 + AF≥0.05 + ALT reads>10.

**Clinical TMB thresholds (FDA):** Low < 10 mut/Mb; High ≥ 10 mut/Mb (pembrolizumab TMB-H indication).

## `results/somatic/signatures/`

| Directory | Description |
|-----------|-------------|
| `filtered_vcfs/` | PASS-filtered VCFs per sample (SigProfiler input) |
| `<PAIR>/matrix_profiles/output/` | SBS/DBS/ID mutational profiles and plots for the pair |
| `<PAIR>/assignment/<PAIR>_SBS96/` | COSMIC v3.4 SBS decomposition + activities + plots |
| `<PAIR>/assignment/<PAIR>_DBS/` | COSMIC DBS decomposition |
| `<PAIR>/assignment/<PAIR>_ID/` | COSMIC INDEL decomposition |

SigProfiler runs once per tumor/normal pair, with its input VCF staged inside the task.

## `results/qc/gene_coverage/`

| File | Description |
|------|-------------|
| `*_gene_coverage.tsv` | Per-gene coverage statistics (mean, min, max, n_regions) for each sample |

Coverage is calculated from mosdepth regional BED files, aggregated per gene symbol. Useful for identifying genes with insufficient coverage (mean < 30X) before reporting variants.

## `results/reporting/`

| File | Description |
|------|-------------|
| `*_molecular_report.html` | Interactive clinical HTML report with KPI dashboard and data tables |
| `*_molecular_report.txt` | Plain-text version of the report |

**Dashboard KPIs:**
- **TMB** (mut/Mb): Tumor Mutational Burden — FDA threshold: ≥10 mut/Mb = TMB-High
- **MSI** (%): Microsatellite Instability fraction — MSIsensor-pro score ≥20% typically MSI-H
- **LOH**: Loss of Heterozygosity fraction (ASCAT)
- **GI**: Genomic Instability index (ASCAT)
- **WGD**: Whole Genome Doubling detected (ASCAT; 0=No, 1=Yes)
- **Ploidy**: Tumour ploidy (ASCAT)
- **Purity**: Tumour cell fraction % (ASCAT)
- **Coverage**: Mean on-target coverage for normal and tumour samples (mosdepth)
