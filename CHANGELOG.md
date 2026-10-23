# Changelog

All notable changes to `nf-sarek-clinical` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Local gnomAD v4.1 population frequencies for germline SNVs (`--gnomad_tsv`).
  `scripts/build_gnomad_panel.sh` builds a ~13 MB subset of gnomAD v4.1 joint
  (exomes+genomes) over the MANE exons (±100 bp) of the panel genes; the new
  `GNOMAD_ANNOTATE` module adds `GNOMAD_*` INFO fields (AF, FAF95 grpmax, ancestry
  group, admixed-American AF, homozygote counts, gnomAD FILTER).
- Canonical `gnomAD_AF` / `gnomAD_AF_source` columns in germline TSVs, used by
  BA1/BS1/PM2, the HTML/PDF reports and the germline–somatic correlation.
- Germline re-annotation with MANE Select transcripts (`--snpeff_mane_data`,
  new `SNPEFF_MANE` module, snpEff 5.2 + `GRCh38.mane.1.2.refseq`). Sarek's
  snpEff `-canon` keeps the longest-CDS transcript, so germline HGVS did not
  follow clinical nomenclature for ~10% of candidate variants (BRCA1 p.Ser1634Gly vs
  p.Ser1613Gly, MUTYH p.Gln335His vs p.Gln338His, MEN1, MET, POLD1, CASR, EPCAM).
- Local ClinVar germline classifications (`--clinvar_tsv`).
  `scripts/build_clinvar_panel.sh` builds a panel subset of the NCBI ClinVar
  GRCh38 VCF (contigs renamed, left-normalised); the new `CLINVAR_ANNOTATE`
  module annotates by allele, so indels are covered. Canonical `ClinVar_sig`,
  `ClinVar_stars`, `ClinVar_revstat`, `ClinVar_id`, `ClinVar_conf`, `ClinVar_dn`,
  `ClinVar_source` columns are read by the ACMG engine, the germline filters, the
  HTML/PDF reports (classification with stars) and the correlation. When set,
  the local release is the only ClinVar source (no dbNSFP 4.1a fallback).
- `ZYGOSITY` column in germline TSVs (HET, HOM, HOM/HEMI on non-PAR chrX/chrY,
  HET(VAF alta)), shown as "Cig." in the HTML and PDF germline tables.
- `QC_FLAG` column in germline TSVs (`DP<20`, `VAF_baja`, and the VCF FILTER
  when not PASS: `CNN_SNP`, `CNN_INDEL`, `LowQual`), shown inside the VAF cell
  of the HTML and PDF germline tables. Non-PASS variants are flagged, not removed.

### Changed
- With `--gnomad_tsv` set, no patient coordinates are sent to Ensembl or gnomAD
  web APIs. Unknown frequency (outside the resource regions, no dbNSFP value) no
  longer triggers PM2; sites failing gnomAD QC in exomes and genomes keep their AF
  for BA1/BS1 but never support PM2.
- Germline multi-allelic sites (more than one ALT, e.g. GT 1/2 in the AR CAG/GGC
  repeats) are no longer reported: they are moved to `germline_all_filtered.tsv`
  with `FILTER_REASON = multiallelic (GT=…; ALT=…)`. Such sites are not reportable
  in a predisposition test and their VAF / ClinVar match were wrong.
- Germline panel: RAD51D and MRE11 (HGNC symbol, not MRE11A) added (103 -> 105
  genes); both were already in the HRR gene set but never reached the germline
  SNV filter. The file is renamed `germline_cancer_predisposition.txt` (no gene
  count in the name). The local gnomAD v4.1 and ClinVar resources are rebuilt
  over the new regions (124 genes) — rebuild them after updating.
- ClinVar review-star threshold (`filters.acmg.clinvar_min_stars`, default 2).
  Below it a ClinVar class is no longer final: P/LP go to the ACMG score with a
  VUS floor (except BA1) and a `ClinVar:LP(1*)` review criterion; B/LB are no
  longer removed by the germline filter and are scored. A ClinVar VUS with >= 2
  stars can now be brought down to LB/B by BS1/BA1 (e.g. MYC p.Asn26Ser, AF 3.2%).
  The `clinvar_benign` filter reason now shows the stars.

### Fixed
- Germline variants with VAF > `vaf_max` (0.75) were discarded as
  `vaf_too_high`, hiding homozygous/hemizygous variants in recessive genes
  (MUTYH, NTHL1, MSH3, MBD4, FANC*, CMMRD) from reports and correlation.
  High VAF is no longer a rejection reason.
- Germline low-VAF floor (0.25 HIGH / 0.35 MODERATE) discarded real
  heterozygotes (e.g. VAF 0.34 at DP 59). Now a single 0.20 floor, applied only
  at DP >= 20 (`vaf_floor_min_dp`), never to ClinVar P/LP (non-conflicting).
- Germline multi-hit filter paired variants anywhere in the same exon and
  always removed in-frame deletions, discarding real variants (AR GGN repeat at
  VAF 0.62). Now only variants <= 30 bp apart pair, and each member is removed
  only if it has VAF < 0.40, FILTER != PASS or GQ < 20 (`multi_hit_min_gq`);
  the reason names the partner variant and distance.
- Multi-value dbNSFP frequencies (`0.53,0.53`) were read as AF 0 by the ACMG
  engine (PM2 applied to common variants).

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
  - `germline_103_cancer_predisposition.txt` (103 genes, ACMG SF v3.2 + ESMO)
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
