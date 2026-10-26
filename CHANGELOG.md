# Changelog

All notable changes to `nf-sarek-clinical` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- One gene table, `assets/gene_panels/genes.tsv` (543 genes; manifest of sources),
  from which `scripts/build_gene_lists.py` writes every panel file and
  `bin/gene_sets.py` (MMR / HRR); `--check` flags stale files. Includes the
  Illumina TruSight Oncology 500 DNA content (523 genes, data sheet M-GL-00173
  v9.0, symbols updated to HGNC; 9 data-sheet typos fixed against an independent
  list).
- Somatic SNV second-hit track (`--second_hit_panel`): germline genes outside the
  core somatic panel (49) where a somatic loss-of-function completes biallelic
  inactivation. Only the germline–somatic correlation and the findings cards read
  it; before, 64 of 105 germline genes (RAD51C/D, MUTYH, MSH3, CDKN2A, VHL, …) could
  never get a somatic SNV second hit.
- Pan-tumor somatic SNV track (`--tso500_panel`): TSO500 genes outside the core
  (466). Tier I/II join the actionable variants and findings (marked TSO500),
  Tier III go to a separate table (collapsed in HTML), Tier IV to the population
  polymorphisms. Tracks reuse the same processes (`conf/modules.config`,
  `task.ext.suffix`). 1128: 88 TSO500 VUS; IMAX1141: 15.
- Gene context for the TSO500 genes (`scripts/build_gene_context.py`): ROLE from
  CancerMine (CC0) only when one-sided (>= 5x, >= 20 citations; 0 opposite roles on
  the 60 curated genes), summary from CIViC (CC0) or NCBI Gene.
- TERT promoter (C228T / C250T): "not assessable" note in HTML and PDF when the
  hotspots lie outside `--panel_bed` (Agilent S07604514 stops ~120 bp short).
- Germline CNV low-mappability flag (`--gcnv_lowmap_bed`, built by
  `scripts/build_lowmap_panel.sh` from GIAB stratifications v3.3
  `GRCh38_lowmappabilityall`, panel MANE exons ±1 kb). The annotated TSV gains
  `LOWMAP` / `LOWMAP_SEG`; segments with >= `lowmap_min_frac` (0.5) of their panel
  exon bases in these regions go to a separate "low-mappability — orthogonal
  confirmation required" table in HTML and PDF and are not used by the
  correlation. CHEK2 exon 14 HOMDEL (pseudogene artefact, 5 of 7 normals, shown as
  "High" and as a first hit) and a FANCD2 HOMDEL in 1229 move there; MLH1
  exon 4–6 deletion in 711 is unchanged.
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
- Non-coding splice candidates (`filters.germline_snv.noncoding_splice`, needs
  SpliceAI): splice-region (+3..+8 / -3..-20, donor 5th base, polypyrimidine
  tract), synonymous and intronic variants of panel genes now reach SpliceAI and
  are kept when `SpliceAI_DS >= 0.2` or ClinVar is P/LP (any stars); otherwise
  they are filtered as `spliceai_below_threshold`. Before, only HIGH / MODERATE
  variants passed the panel filter, so 504 ClinVar >= 2 stars P/LP non-coding
  variants of the panel (e.g. MSH2 c.942+3A>T) could never be reported.

### Changed
- Germline CNV gene list: `exon_annotation.R` reads `--genes` (the germline
  panel, `params.germline_panel`) instead of a fixed list of 91 genes, minus
  `filters.germline_cnv.exclude_genes` (16 somatic oncogenes: AR, BRAF,
  CCND1-3, CCNE1, CDK6, EGFR, ERBB2, FGFR1-3, KRAS, MDM2, MDM4, MYC). MRE11 is
  now annotated (the old list used MRE11A, absent from MANE); CYP2D6 and PRSS1
  are out. The R stops if a panel gene is missing from the MANE exon table.
- Germline CNV events are the gCNV ALT `<DEL>` / `<DUP>` (copy number against
  the contig ploidy) instead of `CN != 2`, in the HTML, PDF and correlation.
  With `CN != 2`, the single male chrX copy read as a deletion and a male chrX
  CN=2 duplication was hidden. Section titles say "deletions / duplications".
- REVEL PP3 / BP4 only when the MANE effect is missense: a synonymous or
  intronic MANE variant could carry the REVEL score of another transcript
  (`REVEL_v1.3_other_tx`; 4 such rows in 1128, all filtered by BA1).
- Germline filter (Step 2): BA1 uses the rule set of the ACMG engine (new
  column `gnomAD_BA1`) and never removes a ClinVar P/LP variant with >=
  `clinvar_min_stars` stars (founder pathogenic variants).
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
- ACMG points follow ClinGen SVI (Tavtigian 2020): PM2 is Supporting (+1);
  REVEL uses the calibrated cut-offs of Pejaver 2022 (PP3 0.644 / 0.773 / 0.932,
  BP4 0.290 / 0.183 / 0.016); the 9-tool panel only counts when REVEL is missing,
  at Supporting; LB is -1..-6 and B <= -7. PVS1 now depends on the gene mechanism
  (`filters.acmg.pvs1_not_applicable`: no PVS1 in oncogenes, POLD1/POLE, HOXB13,
  EPCAM, KIF1B, RECQL, ...) and on the variant type (last exon -> PVS1_Strong /
  PVS1_Moderate, start_lost -> PVS1_Moderate, stop_lost -> PM4). Rare missense
  variants without functional evidence no longer reach LP, and a new null variant
  is LP instead of P. Report methodology text and criteria colours updated.
- BA1 / BS1 use the gnomAD v4.1 max filtering allele frequency across ancestry
  groups (new canonical `gnomAD_FAF95` column) instead of the overall AF, in the
  germline filter and the ACMG engine; PM2_Supporting needs AF and FAF95 below
  `pm2_af`. Variants common in one ancestry group (BRCA2 p.Ile2490Thr, CASR
  p.Thr445Ala, TSC1 p.Lys587Arg, RECQL4 p.Arg522Cys: FAF95 10-20%, mostly admixed
  American) are BA1 even when the overall AF is < 1%. `filters.acmg.frequency_metric`
  ('faf95' / 'af'). Reports show "AF / FAF".
- Local REVEL v1.3 resource (`--revel_tsv`, `scripts/build_revel_panel.py`,
  new `REVEL_ANNOTATE` module): the score of the MANE Select transcript replaces
  the dbNSFP 4.1a value (first transcript). Canonical `REVEL_score` /
  `REVEL_source` columns.
- SpliceAI 1.3.1 on the germline candidates (`--spliceai`, new `SPLICEAI`
  module, `SpliceAI_DS` / `SpliceAI_detail` columns). ClinGen splicing rules
  (Walker 2023): SpliceAI >= 0.2 gives PP3_Supporting (the stronger of REVEL and
  SpliceAI counts), a missense BP4 needs SpliceAI <= 0.1; no PP3 / BP4 for null
  variants. The uncalibrated 9-tool dbNSFP panel is no longer scored. The
  germline ACMG classification no longer depends on dbNSFP when `--revel_tsv`
  is set (dbNSFP still runs for pharmacogenomics and display). HTML report: new
  SpliceAI column.
- ACMG PS1 / PM5 from ClinVar at the same residue (`--clinvar_protein_index`,
  `scripts/build_clinvar_protein_index.py`; new `ClinVar_PS1` / `ClinVar_PM5`
  columns). References: ClinVar P/LP missense with >= `clinvar_min_stars` stars,
  not at exon edges. PS1 (+4) needs SpliceAI < 0.2; PM5 (+2) only without PS1.
  Leave-one-out on ClinVar missense >= 2 stars (germline genes): P/LP classified
  LP rise from 0 to 612 of 2685; no benign variant reaches LP (0 of 1718).
- ClinGen VCEP gene specifications (CSpec registry, Released, 2026-10-23) in
  `bin/acmg_shared.py` (`VCEP_SPECS`, `VCEP_PM1`; `vcep_specs = true`):
  gene-specific BA1 / BS1 / BS1_Supporting (-1) / PM2_Supporting for APC, ATM,
  BRAF, BRCA1, BRCA2, CDH1, DICER1, ENG, KRAS, MLH1, MSH2, MSH6, PALB2, PMS2,
  PTEN, TP53, VHL, and PM1 (+2 / +1; PP3 + PM1 at most 4) for TP53, PTEN, DICER1,
  ENG, VHL, KRAS, BRAF (TP53 / VHL hotspot counts from cancerhotspots.org v2,
  `scripts/build_vcep_hotspots.py`). New `ACMG_spec` column (rule set used).
  `scripts/build_gnomad_panel.sh` adds grpmax AF / AC / AN (`GNOMAD_AF_GRPMAX`,
  `GNOMAD_AC_GRPMAX`, `GNOMAD_AN_GRPMAX`; older resources use the general rules).
  Leave-one-out on ClinVar >= 2 stars in the 17 genes (30,919 variants): B/LB
  classified B/LB rise from 1,291 to 2,379 of 14,850; P/LP classified P/LP from
  13,830 to 13,938 of 16,069 (PM1: KRAS, TP53, VHL, BRAF, PTEN); P/LP classified
  B/LB from 22 to 27 (ClinVar keeps their class in the pipeline).

### Fixed
- GATK somatic CNV: segments without heterozygous SNPs (MINOR_AF NA) are
  classified by copy number instead of "UNKNOWN"; AMP+LOH counts as an
  amplification (IMAX1141: MYC CN >= 5 now reported) but not as LOH (a CN 3 gain has
  a minor allele fraction ~0.33 with both alleles present); the R script reads its
  genes from `--cnv_gatk_genes` (germline + somatic CNV panels from genes.tsv) instead
  of a fixed list that missed 8 of the 38 somatic CNV genes and spelled MRE11 as
  MRE11A; the EXON column is filled (it read a column that does not exist).
- A patient with two tumors and one normal kept only one tumor in the GATK somatic
  CNV (join by patient): every tumor is now paired with the normal.
- A pair without ASCAT or without CNVkit no longer stops the run: without ASCAT the
  somatic CNV report uses CNVkit only, without CNVkit it is skipped (both logged).
- One somatic CNV harmonisation (`bin/cnv_shared.py`) for HTML and PDF instead of
  two inline copies; the unused `bin/harmonize_cnvs.py` is removed.
- MANE exon table: `exons_mane_v1.2.txt` (`scripts/build_mane_exons.py`, MANE v1.2
  Select + Plus Clinical) recommended for `--gcnv_mane_txt`; germline and somatic
  CNV events of 4 normals / 4 tumors unchanged.
- Somatic population frequency from the local gnomAD v4.1 resource and AMP
  Tier IV. The somatic Mutect2 VCF goes through `GNOMAD_ANNOTATE` (`--gnomad_tsv`)
  and the somatic TSV gains `gnomAD_AF`, `gnomAD_AF_source`, `gnomAD_FAF95`,
  `gnomAD_NHOMALT`, `gnomAD_POP_AF` and `POP_COMMON`. A variant with FAF95 (else AF)
  >= `filters.somatic_snv.max_pop_af` (0.01) is Tier IV (AMP/ASCO/CAP 2017,
  MAF >= 1%): listed in a separate table in HTML and PDF and not used as a somatic
  finding, double hit or second hit. 1128 BRCA2 p.Thr1915Met (gnomAD AF 2.65%, 679
  homozygotes, VAF 0.125) moves from Tier III to Tier IV; no other variant of 6 test
  patients reaches 1%. Somatic tables show a gnomAD column. Without `--gnomad_tsv`
  the run warns and no variant is Tier IV.
- dbNSFP annotation removed from the somatic branch (`SNPSIFT_DBNSFP_SOMATIC`):
  no somatic step read its fields; saves 2–4 min per pair. The somatic TSV no
  longer has `dbNSFP_*` columns and `somatic/snv/dbnsfp/` is no longer published.
- MANE Select nomenclature everywhere. The somatic Mutect2 VCF is now
  re-annotated by `SNPEFF_MANE` too (IMAX1141 BRAF p.Lys641Glu → p.Lys601Glu;
  8 of 69 somatic HGVS change in 6 test patients: BRAF, CTNNB1, POLD1, ROS1,
  CHEK2). The MANE database also holds MANE Plus Clinical transcripts and snpEff
  listed them first for BRAF, FGFR2, KRAS, MUTYH and SMARCA4, so germline HGVS and
  the ClinVar protein index (PS1 / PM5) used them (BRAF V640, MUTYH G396D; the
  MUTYH p.Gln338His above is NM_001128425.2, MANE Select is p.Gln310His).
  `SNPEFF_MANE` lists them after MANE Select (`--mane_plus_clinical`,
  `assets/mane/mane_v1.2_plus_clinical.txt`) and keeps their HGVS in
  `INFO/MANE_PLUS_CLINICAL`; `scripts/build_clinvar_protein_index.py` takes
  `--mane-plus-clinical` (rebuilt index: 247 rows of BRAF, FGFR2, GNAS and MUTYH
  renumbered, one CDKN2A p14ARF-only row dropped). HTML and PDF variant tables
  show the transcript under HGVS_c and the MANE Plus Clinical HGVS under HGVS_p
  when it differs.
- Reports mixed patients when a run had more than one pair: the germline and
  somatic tables were attached with `.combine()` (every pair × every table, N²
  report tasks) and the PDF took the CIViC table from a loose channel in
  arrival order; seen with two pairs: two PDFs labelled with one pair, one of
  them with the other's CIViC table, and none for the second pair. Each pair
  now gets one HTML and one PDF task with its own ASCAT, HRD, somatic, CIViC
  and (through its normal) germline files, and the reports read those staged
  files instead of searching the outdir by name prefix. A pair without ASCAT
  output still gets its reports (NA metrics, warning in the log).
- Mutational signatures could report a previous run's result: the input VCFs
  were copied to a shared `${workDir}/mutsig_vcfs`, where SigProfiler also wrote
  its `input/`/`output/`, and a `-resume` run reused them (seen: activities of
  the previous run and a published SBS96 matrix of 6372 mutations for a VCF of
  3192). With several pairs in a run, the PDF took the first sample's row and
  the plots were drawn for the first sample only. SigProfiler now runs once per
  pair with its VCF staged in the task; outputs go to
  `somatic/signatures/<pair>/`, and every reader picks its pair by exact name.
- AMP tier: one implementation (`bin/amp_tier.py`) for the HTML table, the PDF
  table and the PDF findings card, which disagreed (level A resistance was IA in
  HTML and IIA in PDF; HTML also produced "IID"/"IV"). Only accepted CIViC
  evidence with direction SUPPORTS sets the tier — DOES_NOT_SUPPORT records were
  counted as support, and submitted evidence weighed like accepted. The tier now
  depends on the tumor type (new `--tumor_type` / `--tumor_types`): level A/B in
  the patient's tumor type or pan-cancer → IA/IB, C → IIC, D/E → IID; A/B from
  another tumor type → IIC. Without a tumor type the highest tier is IIC.
  Evidence is shown grouped by significance with level, disease, evidence
  against and the number of submitted records; `%2C` in disease names decoded.
- Four latent somatic bugs, no effect on the current runs:
  - consensus used `bcftools isec -n=N` ("exactly N"): with
    `min_callers = 1` it kept only the Mutect2-private calls. Now `-n+N`, and
    any value other than 1 or 2 stops the run;
  - the CNV gene BED was read with its first interval as column names (lost);
    a header is now detected from non-numeric start/end, and a headerless
    file is no longer rejected;
  - the HTML report took the SBS96 activities as ID signatures whenever the
    output path or run name contained "ID" (`-path "*ID*"`); HTML, PDF and the
    plot step now match the `<run>_SBS96|_DBS|_ID/Assignment_Solution` directory;
  - the somatic SNV table took the first snpEff ANN for gene, impact and HGVS,
    not the one that met the panel/impact filter. Each row now uses the first
    ANN with a panel gene and HIGH/MODERATE impact.
- Somatic steps read the tumor at a fixed VCF position (`GEN[1]` in TMB and the
  SNV table, a `grep '_T'` on sample names in signatures). Mutect2 writes its
  samples in alphabetical order, so a tumor named before its normal would have
  had TMB and VAF taken from the normal. New `MUTECT2_SAMPLE_ORDER` step
  rewrites each Mutect2 VCF as normal,tumor from its `##normal_sample` /
  `##tumor_sample` header lines (error if missing), and signatures take the
  tumor from `##tumor_sample`. Signatures no longer need CNVkit output: the
  tumor name used to come from the CNVkit files, and a pair without them was
  dropped silently.
- Somatic consensus and the signature/HRD input kept Mutect2 records whose
  FILTER merely contained `haplotype` (`bcftools view -f PASS,haplotype`), so
  `haplotype;orientation`, `clustered_events;haplotype`, `weak_evidence`…
  got through. Both now keep FILTER exactly `PASS` or exactly `haplotype`, the
  rule TMB and the somatic SNV table already used: 23–91 artefact records
  fewer per patient in the consensus.
- Somatic consensus dropped every Mutect2 MNP (and the matching Strelka2 SNVs):
  Mutect2 writes adjacent substitutions as one record (CC>TT), Strelka2 as
  separate SNVs, and the exact-allele `bcftools isec` never matched them. Both
  callers are now atomized (`bcftools norm -a --atom-overlaps .`) before the
  intersection; +3 to +14 consensus variants per patient, none lost. TMB and
  signatures count each base of the MNP, as Strelka2 does.
- CIViC annotation was empty in every run since 2026-01-09: civicpy 5.2.0 could
  not unpickle the nightly cache (`Can't get attribute 'ClingenCode'`) and the
  task copied the unannotated VCF, so every somatic variant read as Tier III.
  Now civicpy 5.4.0 (image `civicpy:5.4.0`), a frozen dated cache passed with
  `--civic_cache` (`scripts/build_civic_cache.sh`; no download, no staleness
  refresh at run time), and any civicpy failure or missing CIVIC header fails
  the task (`CIVIC_EXTRACT` too). Evidence status is now `accepted` +
  `submitted` (before, `-i submitted` alone left accepted evidence out).
- Somatic GATK CNV with CRAM input: CollectReadCounts and CollectAllelicCounts
  now read a CRAM with `--cram_fasta` (BAM keeps `somatic_fasta`). Since the
  germline CNV CRAM fix the somatic steps also received CRAMs and failed with
  "incompatible contigs"; before it they silently never ran on CRAM-only Sarek
  output. A CRAM sample without `--cram_fasta` stops the run with a clear error.
  1128 from CRAM: read/allelic counts, denoised CR, called segments and the
  annotated TSV identical to the historic BAM run.
- Germline CNV report texts: the HTML note said "NP >= 5 and QA >= 10" and the
  PDF "105-gene panel" and "CN != 2", while the filter is QA >=
  `filters.germline_cnv.min_qa` (4), no NP minimum, ALT `<DEL>`/`<DUP>`, panel
  minus somatic oncogenes. Notes, "no events" boxes and methodology tables now
  print the threshold actually applied.
- Germline CNV readers: a non-integer QA (`.`, empty, text) made
  `germline_somatic_correlate.py` crash, and since the reports call it with
  `|| true` the whole germline-somatic section disappeared silently. QA is now
  read as a number (`4.0` = 4) in correlation, HTML and PDF; a non-numeric QA
  skips the segment (below threshold) and the correlation logs a warning.
- Germline CNV input: Sarek `recalibrated/` is searched for `*.recal.bam` and
  `*.recal.cram` (BAM used when both exist). Before, only BAM was searched, so
  with CRAM output gCNV never ran and the run ended SUCCESS with no germline
  CNV. A CRAM is read with `--cram_fasta` (the FASTA Sarek wrote it with,
  plus `.fai` and `.dict`); the run fails if a normal is CRAM without it, if
  an index is missing, or if no normal BAM/CRAM is found.
- Germline CNV (`MANE_ANNOTATE_GERMLINE_CNV`): an `exon_annotation.R` failure
  now fails the task. Before, the error was only logged and the raw segments
  VCF lines were copied without a header, so every reader saw CN=2 and the
  reports said "no germline CNV". The task now copies only
  `<sample>_segments.txt` and checks the columns CHROM, POS, END, CN, QA, gen.
- `--skip_signatures true` now skips `MUTSIG_SIGPROFILER` (before, it only
  stopped the report from waiting for it, and the process still ran).
- SigProfiler image `sigprofiler:0.1.8`: SigProfilerMatrixGenerator 1.3.6 ->
  1.4.0. The GRCh38 checksums shipped with 1.3.6 no longer match the files
  served by the Alexandrov lab FTP, so a freshly built image failed with
  "The specified genome GRCh38 has not been installed".
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
