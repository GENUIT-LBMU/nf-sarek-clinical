/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SAREK-CLINICAL MAIN WORKFLOW
    ==============================
    Orchestrates all postprocessing modules after nf-core/sarek completes.

    Expected Sarek output structure:
      <sarek_outdir>/
        annotation/
          haplotypecaller/<sample>/*.vcf.gz       (snpEff annotated germline)
          mutect2/<pair>/*mutect2.filtered_snpEff.ann.vcf.gz
        variant_calling/
          cnvkit/<pair>/*.somatic.call.cns  (or *.call.cns)
          ascat/<pair>/*.metrics.txt  +  *segments.txt
          msisensorpro/<pair>/<pair>
        preprocessing/
          recalibrated/<sample>/<sample>.recal.bam  +  .bam.bai
        reports/
          mosdepth/<sample>/<sample>.recal.regions.bed.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GERMLINE_POSTPROCESS  } from '../subworkflows/local/germline/main'
include { SOMATIC_POSTPROCESS   } from '../subworkflows/local/somatic/main'
include { CNV_GERMLINE_GATK     } from '../subworkflows/local/cnv_germline_gatk/main'
include { CNV_SOMATIC_GATK      } from '../subworkflows/local/cnv_somatic_gatk/main'
include { GENE_COVERAGE         } from '../modules/local/gene_coverage/main'
include { REPORT_MODEL          } from '../modules/local/report_model/main'
include { FINAL_REPORT          } from '../modules/local/final_report/main'
include { CLINICAL_PDF_REPORT   } from '../modules/local/clinical_pdf_report/main'

workflow SAREK_CLINICAL {

    // ── Reference file channels ───────────────────────────────────────────────
    def ch_snpsift_jar   = params.snpsift_jar ? file(params.snpsift_jar)          : error("--snpsift_jar required")
    def ch_panel_bed     = params.panel_bed   ? file(params.panel_bed)   : []
    def ch_gene_coords   = params.gene_coords_bed ? file(params.gene_coords_bed) : []

    // Local gnomAD panel resource (required for germline SNVs, optional for somatic Tier IV).
    def gnomad_regions = params.gnomad_regions ?: (params.gnomad_tsv ? params.gnomad_tsv.toString().replaceAll(/\.tsv\.gz$/, '_regions.bed.gz') : null)
    def ch_gnomad        = params.gnomad_tsv
        ? [ file(params.gnomad_tsv,            checkIfExists: true),
            file("${params.gnomad_tsv}.tbi",   checkIfExists: true),
            file(gnomad_regions,               checkIfExists: true),
            file("${gnomad_regions}.tbi",      checkIfExists: true) ]
        : []
    // Germline SNV classification uses only local, versioned resources (no dbNSFP,
    // no web APIs): gnomAD (frequency), ClinVar and REVEL are required
    if (!params.skip_germline_snv) {
        if (!params.gnomad_tsv)  error("--gnomad_tsv is required for germline SNVs (scripts/build_gnomad_panel.sh)")
        if (!params.clinvar_tsv) error("--clinvar_tsv is required for germline SNVs (scripts/build_clinvar_panel.sh)")
        if (!params.revel_tsv)   error("--revel_tsv is required for germline SNVs (scripts/build_revel_panel.py)")
    }
    if (!params.gnomad_tsv && !params.skip_somatic_snv) {
        log.warn "--gnomad_tsv not set: somatic variants have no population frequency, so common polymorphisms are not set to AMP Tier IV."
    }

    // snpEff MANE database (optional). Empty list => keep sarek's snpEff -canon annotation.
    def ch_snpeff_mane   = params.snpeff_mane_data
        ? [ file(params.snpeff_mane_data, checkIfExists: true), params.snpeff_mane_db,
            file(params.mane_plus_clinical, checkIfExists: true) ]
        : []
    if (params.snpeff_mane_data && !file("${params.snpeff_mane_data}/${params.snpeff_mane_db}/snpEffectPredictor.bin").exists()) {
        error("--snpeff_mane_data: ${params.snpeff_mane_data}/${params.snpeff_mane_db}/snpEffectPredictor.bin not found")
    }
    if (!params.snpeff_mane_data && !(params.skip_germline_snv && params.skip_somatic_snv)) {
        log.warn "--snpeff_mane_data not set: germline and somatic HGVS use sarek's snpEff -canon transcript (longest CDS), which may differ from MANE Select."
    }
    if (params.snpeff_mane_data && !params.snpeff_mane_db.contains('mane.1.2.') && file(params.mane_plus_clinical).name == 'mane_v1.2_plus_clinical.txt') {
        log.warn "--snpeff_mane_db ${params.snpeff_mane_db} is not MANE v1.2: set --mane_plus_clinical to the MANE Plus Clinical list of that release."
    }

    // Local ClinVar panel resource (required for germline SNVs).
    def ch_clinvar       = params.clinvar_tsv
        ? [ file(params.clinvar_tsv,           checkIfExists: true),
            file("${params.clinvar_tsv}.tbi",  checkIfExists: true) ]
        : []

    // Local REVEL v1.3 panel resource (required for germline SNVs).
    def ch_revel         = params.revel_tsv
        ? [ file(params.revel_tsv,             checkIfExists: true),
            file("${params.revel_tsv}.tbi",    checkIfExists: true) ]
        : []

    // SpliceAI on germline candidates (optional). Empty list => no splicing evidence.
    def spliceai_fasta = params.spliceai_fasta ?: params.gcnv_fasta
    def ch_spliceai      = (params.spliceai && spliceai_fasta)
        ? [ file(spliceai_fasta,               checkIfExists: true),
            file("${spliceai_fasta}.fai",      checkIfExists: true) ]
        : []
    if (params.spliceai && !spliceai_fasta && !params.skip_germline_snv) {
        log.warn "--spliceai true but no --spliceai_fasta / --gcnv_fasta: SpliceAI skipped (no PP3/BP4 splicing evidence)."
    }

    // ClinVar protein index for PS1 / PM5 (optional). [] => no PS1 / PM5.
    def ch_protein_index = params.clinvar_protein_index
        ? file(params.clinvar_protein_index, checkIfExists: true)
        : []
    if (!params.clinvar_protein_index && !params.skip_germline_snv) {
        log.warn "--clinvar_protein_index not set: no ACMG PS1 / PM5. Build it with scripts/build_clinvar_protein_index.py."
    }

    def ch_germline_panel    = file(params.germline_panel)
    def ch_somatic_cnv_panel = file(params.somatic_cnv_panel)
    def ch_somatic_snv_panel = file(params.somatic_snv_panel)
    def ch_second_hit_panel  = params.second_hit_panel ? file(params.second_hit_panel, checkIfExists: true) : []
    def ch_tso500_panel      = params.tso500_panel     ? file(params.tso500_panel,     checkIfExists: true) : []
    def ch_pharmaco_variants = file(params.pharmaco_variants)

    def sarek_out = file(params.sarek_outdir)

    // ── Auto-discover Sarek outputs ───────────────────────────────────────────

    // Germline HaplotypeCaller VCFs (snpEff annotated)
    def ch_hc_vcfs = Channel.fromPath(
        "${sarek_out}/annotation/haplotypecaller/**/*.vcf.gz"
    ).map { vcf ->
        def sample = vcf.parent.name
        [sample, vcf]
    }

    // Somatic Mutect2 VCFs (snpEff annotated)
    def ch_mutect2_vcfs = Channel.fromPath(
        "${sarek_out}/annotation/mutect2/**/*mutect2.filtered_snpEff.ann.vcf.gz"
    ).map { vcf ->
        def pair = vcf.parent.name
        [pair, vcf]
    }

    // Somatic Strelka2 VCFs (SNVs + indels, from variant_calling — not snpEff annotated).
    // When sarek was run without --tools strelka, these dirs don't exist; Channel.fromPath
    // with checkIfExists:false naturally emits an empty channel — no .ifEmpty needed.
    // (.ifEmpty { Channel.empty() } emits a Channel-object as a value, which breaks the
    // downstream .filter in SOMATIC_POSTPROCESS.)
    def ch_strelka2_snvs = Channel.fromPath(
        "${sarek_out}/variant_calling/strelka/**/*.somatic_snvs.vcf.gz",
        checkIfExists: false
    ).map { vcf -> [vcf.parent.name, vcf] }

    def ch_strelka2_indels = Channel.fromPath(
        "${sarek_out}/variant_calling/strelka/**/*.somatic_indels.vcf.gz",
        checkIfExists: false
    ).map { vcf -> [vcf.parent.name, vcf] }

    // CNVkit .somatic.call.cns (fallback to .call.cns)
    def ch_cns = Channel.fromPath(
        "${sarek_out}/variant_calling/cnvkit/**/*.{somatic.call.cns,call.cns}"
    ).map { cns ->
        def pair   = cns.parent.name
        def parts  = pair.split('_vs_')
        def tumor  = parts.size() >= 2 ? parts[0] : pair
        def normal = parts.size() >= 2 ? parts[1] : 'UNKNOWN'
        [pair, tumor, normal, cns]
    }.groupTuple(by: 0)
    .map { pair, tumors, normals, cnss ->
        def best = cnss.find { it.name.contains('somatic.call.cns') } ?: cnss[0]
        [pair, tumors[0], normals[0], best]
    }

    // ASCAT metrics + segments
    def ch_ascat_metrics  = Channel.fromPath(
        "${sarek_out}/variant_calling/ascat/**/*.metrics.txt"
    ).map { f -> [f.parent.name, f] }

    def ch_ascat_segments = Channel.fromPath(
        "${sarek_out}/variant_calling/ascat/**/*segment*.txt"
    ).map { f -> [f.parent.name, f] }

    def ch_ascat = ch_ascat_metrics.join(ch_ascat_segments, remainder: true)
        .map { pair, met, seg -> [pair, seg ?: file('NO_ASCAT_SEG'), met ?: file('NO_ASCAT_MET')] }

    // Tumor sample names (for mutsig tumor column detection)

    // MSIsensor-pro result (file named as the pair) and mosdepth summary per sample (reports)
    def ch_msi = Channel.fromPath(
        "${sarek_out}/variant_calling/msisensorpro/*/*", checkIfExists: false
    ).filter { f -> f.name == f.parent.name }.map { f -> [f.name, f] }
    def ch_mosdepth_summary = Channel.fromPath(
        "${sarek_out}/reports/mosdepth/*/*.mosdepth.summary.txt", checkIfExists: false
    ).map { f -> [f.parent.name, f] }

    // mosdepth regions BED per sample
    def ch_mosdepth = Channel.fromPath(
        "${sarek_out}/reports/mosdepth/**/*.recal.regions.bed.gz"
    ).map { f ->
        def sample = f.parent.name
        [sample, f]
    }

    // Sarek recalibrated alignments (needed for GATK CNV modules). BAM is used
    // when present (pipeline started from BAMs, or save_output_as_bam); CRAM
    // only for samples without a BAM. The index must sit next to the file.
    // Resolved here (not in channel operators) so a validation error is reported alone.
    def recal_alns = files("${sarek_out}/preprocessing/recalibrated/**/*.recal.{bam,cram}")
        .groupBy { it.parent.name }
        .collect { sample, alns ->
            def aln = alns.find { it.name.endsWith('.bam') } ?: alns[0]
            def idx = file(aln.toString() + (aln.name.endsWith('.bam') ? '.bai' : '.crai'))
            if (!idx.exists()) error("Index not found for ${aln}: expected ${idx}")
            [sample, aln, idx]
        }
        .sort { it[0] }
    def ch_bams = Channel.fromList(recal_alns)

    // ── GERMLINE POSTPROCESSING (SNV + Pharmaco) ──────────────────────────────
    if (!params.skip_germline_snv || !params.skip_pharmaco) {
        GERMLINE_POSTPROCESS(
            ch_hc_vcfs,
            ch_snpsift_jar,
            ch_germline_panel,
            ch_pharmaco_variants,
            ch_gnomad,
            ch_clinvar,
            ch_snpeff_mane,
            ch_revel,
            ch_spliceai,
            ch_protein_index
        )
    }

    // ── CRAM reference for the GATK CNV steps ───────────────────────────────
    // A CRAM is decoded with the FASTA Sarek wrote it with (--cram_fasta, plus
    // .fai and .dict). Without --cram_fasta (BAM input) each step stages its own
    // FASTA as a placeholder.
    def cram_ref = { fasta, fai, dict ->
        if (!params.cram_fasta) return [fasta, fai, dict]
        def c_fai  = file("${params.cram_fasta}.fai")
        def c_dict = file(params.cram_fasta.toString().replaceAll(/\.(fa|fasta)$/, '') + '.dict')
        if (!c_fai.exists() || !c_dict.exists())
            error("--cram_fasta needs ${c_fai} and ${c_dict} next to it")
        [file(params.cram_fasta), c_fai, c_dict]
    }
    def check_cram = { alns, skip_flag ->
        alns.each { sample, aln, idx ->
            if (aln.name.endsWith('.cram') && !params.cram_fasta)
                error("${aln} is a CRAM: set --cram_fasta to the FASTA Sarek wrote it with (or ${skip_flag} true)")
        }
    }

    // ── GERMLINE CNV – GATK gCNV ─────────────────────────────────────────────
    if (!params.skip_germline_cnv) {
        // Validate required params
        if (!params.gcnv_ploidy_model)  error("--gcnv_ploidy_model is required when --skip_germline_cnv false")
        if (!params.gcnv_pon_dir)       error("--gcnv_pon_dir is required when --skip_germline_cnv false")
        if (!params.gcnv_interval_list) error("--gcnv_interval_list is required when --skip_germline_cnv false")
        if (!params.gcnv_fasta)         error("--gcnv_fasta is required when --skip_germline_cnv false")
        if (!params.gcnv_dict)          error("--gcnv_dict is required when --skip_germline_cnv false")
        if (!params.gcnv_mane_txt)      error("--gcnv_mane_txt is required when --skip_germline_cnv false")
        if (!params.gcnv_lowmap_bed)    error("--gcnv_lowmap_bed is required when --skip_germline_cnv false (build it with scripts/build_lowmap_panel.sh)")
        if (!file(params.gcnv_lowmap_bed).exists()) error("--gcnv_lowmap_bed not found: ${params.gcnv_lowmap_bed}")

        def gcnv_fasta = file(params.gcnv_fasta)
        def gcnv_fai   = file("${params.gcnv_fasta}.fai")

        // Only process normal samples for germline CNV
        def normal_alns = recal_alns.findAll { sample, aln, idx -> sample ==~ /.*_N\d*$/ }
        if (!normal_alns)
            error("No normal BAM/CRAM found in ${sarek_out}/preprocessing/recalibrated; use --skip_germline_cnv true")
        check_cram(normal_alns, '--skip_germline_cnv')
        ch_normal_bams = Channel.fromList(normal_alns)
        def (cram_fasta, cram_fai, cram_dict) = cram_ref(gcnv_fasta, gcnv_fai, file(params.gcnv_dict))

        CNV_GERMLINE_GATK(
            ch_normal_bams,
            file(params.gcnv_interval_list),
            file(params.gcnv_ploidy_model),
            file(params.gcnv_pon_dir),
            gcnv_fasta,
            gcnv_fai,
            file(params.gcnv_dict),
            file(params.gcnv_mane_txt),
            cram_fasta,
            cram_fai,
            cram_dict,
            ch_germline_panel,
            params.filters.germline_cnv.exclude_genes,
            file(params.gcnv_lowmap_bed)
        )
    }

    // ── SOMATIC POSTPROCESSING (SNV + Consensus + ASCAT/CNVkit + TMB + Signatures + HRD) ─
    // Reference FASTA for bcftools norm in consensus module
    // --skip_strelka_consensus: every pair uses Mutect2 PASS only (no Strelka2 intersection)
    if (params.skip_strelka_consensus) {
        ch_strelka2_snvs   = Channel.empty()
        ch_strelka2_indels = Channel.empty()
    }
    def has_strelka = !params.skip_strelka_consensus &&
        files("${sarek_out}/variant_calling/strelka/**/*.somatic_snvs.vcf.gz").size() > 0
    if (has_strelka && !params.gcnv_fasta)
        error("--gcnv_fasta (reference FASTA + .fai) is required for the Mutect2 + Strelka2 consensus (bcftools norm), or set --skip_strelka_consensus true")
    def consensus_fasta     = params.gcnv_fasta ? file(params.gcnv_fasta) : file('NO_FASTA')
    def consensus_fasta_fai = params.gcnv_fasta ? file("${params.gcnv_fasta}.fai") : file('NO_FASTA_FAI')

    def ch_hrd_tsv = Channel.empty()
    if (!params.skip_somatic_snv || !params.skip_somatic_cnv || !params.skip_tmb || !params.skip_signatures || !params.skip_hrd) {
        // CIViC annotation runs inside SOMATIC_POSTPROCESS with a frozen cache (no network)
        if (!params.civic_cache) error("--civic_cache is required for somatic postprocessing (build it with scripts/build_civic_cache.sh)")
        if (!file(params.civic_cache).exists()) error("--civic_cache not found: ${params.civic_cache}")
        if (params.tumor_types && !file(params.tumor_types).exists())
            error("--tumor_types not found: ${params.tumor_types}")
        if (!params.tumor_type && !params.tumor_types)
            log.warn "No --tumor_type / --tumor_types: AMP tiers are capped at IIC (CIViC evidence treated as from another tumor type)"
        def min_callers = params.filters.consensus.min_callers
        if (!(min_callers in [1, 2]))
            error("filters.consensus.min_callers must be 1 (all Mutect2 PASS) or 2 (Mutect2 and Strelka2), got ${min_callers}")
        SOMATIC_POSTPROCESS(
            ch_mutect2_vcfs,
            ch_strelka2_snvs,
            ch_strelka2_indels,
            ch_cns,
            ch_ascat,
            ch_gnomad,
            ch_snpsift_jar,
            ch_somatic_snv_panel,
            ch_second_hit_panel,
            ch_tso500_panel,
            ch_somatic_cnv_panel,
            ch_germline_panel,
            ch_gene_coords,
            ch_panel_bed,
            consensus_fasta,
            consensus_fasta_fai,
            file(params.civic_cache),
            ch_snpeff_mane
        )
        if (!params.skip_hrd) {
            ch_hrd_tsv = SOMATIC_POSTPROCESS.out.hrd_tsv
        }
    }

    // ── SOMATIC CNV – GATK PAIRED WORKFLOW ───────────────────────────────────
    if (!params.skip_somatic_cnv_gatk) {
        // Validate required params
        if (!params.somatic_pon_hdf5)      error("--somatic_pon_hdf5 is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_interval_list) error("--somatic_interval_list is required when --skip_somatic_cnv_gatk false")
        if (!params.common_snps_vcf)       error("--common_snps_vcf is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_fasta)         error("--somatic_fasta is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_dict)          error("--somatic_dict is required when --skip_somatic_cnv_gatk false")
        if (!params.gcnv_mane_txt)         error("--gcnv_mane_txt is required when --skip_somatic_cnv_gatk false")

        def somatic_fasta = file(params.somatic_fasta)
        def somatic_fai   = file("${params.somatic_fasta}.fai")
        def snps_tbi      = file("${params.common_snps_vcf}.tbi")
        check_cram(recal_alns, '--skip_somatic_cnv_gatk')
        def (cram_fasta, cram_fai, cram_dict) = cram_ref(somatic_fasta, somatic_fai, file(params.somatic_dict))

        CNV_SOMATIC_GATK(
            ch_bams,
            file(params.somatic_interval_list),
            file(params.somatic_pon_hdf5),
            file(params.common_snps_vcf),
            snps_tbi,
            somatic_fasta,
            somatic_fai,
            file(params.somatic_dict),
            cram_fasta,
            cram_fai,
            cram_dict,
            file(params.gcnv_mane_txt),
            ch_ascat_metrics.map { pair, met -> [pair.split('_vs_')[0], met] }   // [tumor, metrics]
        )
    }

    // ── GENE COVERAGE ─────────────────────────────────────────────────────────
    if (!params.skip_coverage) {
        GENE_COVERAGE(
            ch_mosdepth,
            ch_germline_panel,
            params.gcnv_mane_txt ? file(params.gcnv_mane_txt) : file('NO_MANE')
        )
    }

    // ── REPORT INPUTS: one tuple per pair, every file looked up by its key ───
    // Pairs come from Mutect2 and ASCAT. Each lookup channel is collected
    // into a map (key -> file), so the reports start once every upstream
    // result exists, a missing file becomes a placeholder (NO_*) and a normal
    // shared by two tumors is found for both. Keys: pair, tumor or normal sample.
    // Every report input is staged, so -resume reruns a report when any of
    // its inputs changes (no lookup in --outdir).
    def somatic_ran  = !params.skip_somatic_snv || !params.skip_somatic_cnv || !params.skip_tmb || !params.skip_signatures || !params.skip_hrd
    def germline_ran = !params.skip_germline_snv || !params.skip_pharmaco
    def by_key = { ch -> ch.toList().map { l -> [l.collectEntries { [(it[0]): it[1]] }] } }
    def som    = { name -> somatic_ran ? SOMATIC_POSTPROCESS.out."${name}" : Channel.empty() }

    def ch_pairs = ch_mutect2_vcfs.map { pair, _vcf -> pair }
        .mix(ch_ascat_metrics.map { pair, _met -> pair })
        .unique()

    def lookups = [
        ascat:     ch_ascat_metrics,
        hrd:       ch_hrd_tsv,
        germ:      germline_ran ? GERMLINE_POSTPROCESS.out.snv_tsv : Channel.empty(),
        som:       som('somatic_snv_tsv'),
        civic:     som('civic_tsv'),
        som2:      som('second_hit_tsv'),
        tso:       som('tso500_tsv'),
        tso_civic: som('tso500_civic_tsv'),
        tmb:       som('tmb_tsv'),
        msi:       ch_msi,
        mosdepth:  ch_mosdepth_summary,
        cnv:       som('cnv_tsv'),
        cnv_gatk:  params.skip_somatic_cnv_gatk ? Channel.empty()
                       : CNV_SOMATIC_GATK.out.annotated_tsv.map { _pid, tumor, tsv -> [tumor, tsv] },
        gcnv:      params.skip_germline_cnv ? Channel.empty() : CNV_GERMLINE_GATK.out.annotated_tsv,
        pgx:       germline_ran ? GERMLINE_POSTPROCESS.out.pharmaco : Channel.empty(),
        gene_cov:  params.skip_coverage ? Channel.empty() : GENE_COVERAGE.out.tsv,
        sig:       som('signatures_dir'),
        sig_prof:  som('signatures_profiles'),
    ]
    def names = lookups.keySet() as List
    def ch_lookup = names.drop(1).inject(by_key(lookups[names[0]])) { acc, n -> acc.combine(by_key(lookups[n])) }
        .map { maps -> [names, maps].transpose().collectEntries() }

    def ch_report_input = ch_pairs
        .combine(ch_lookup)
        .map { pair, m ->
            def tumor  = pair.contains('_vs_') ? pair.split('_vs_')[0] : pair
            def normal = pair.contains('_vs_') ? pair.split('_vs_')[1] : ''
            if (!m.ascat[pair]) log.warn "${pair}: no ASCAT metrics, the report shows NA for purity, ploidy, LOH, GI and WGD"
            [pair,
             m.ascat[pair]        ?: file('NO_ASCAT'),
             m.hrd[pair]          ?: file('NO_HRD'),
             m.germ[normal]       ?: file('NO_GERMLINE'),
             m.som[pair]          ?: file('NO_SOMATIC'),
             m.civic[pair]        ?: file('NO_CIVIC'),
             m.som2[pair]         ?: file('NO_SECOND_HIT'),
             m.tso[pair]          ?: file('NO_TSO500'),
             m.tso_civic[pair]    ?: file('NO_TSO500_CIVIC'),
             m.tmb[pair]          ?: file('NO_TMB'),
             m.msi[pair]          ?: file('NO_MSI'),
             m.mosdepth[normal]   ?: file('NO_COV_N'),
             m.mosdepth[tumor]    ?: file('NO_COV_T'),
             m.cnv[pair]          ?: file('NO_CNV_SOMATIC'),
             m.cnv_gatk[tumor]    ?: file('NO_CNV_GATK_SOMATIC'),
             m.gcnv[normal]       ?: file('NO_CNV_GATK_GERMLINE'),
             m.pgx[normal]        ?: file('NO_PHARMACO'),
             m.gene_cov[normal]   ?: file('NO_GENE_COV_N'),
             m.gene_cov[tumor]    ?: file('NO_GENE_COV_T'),
             m.sig[pair]          ?: file('NO_SIG'),
             m.sig_prof[pair]     ?: file('NO_SIG_PROFILES')]
        }

    // Versions and resources shown in both reports (bin/report_shared.py).
    // code_md5 covers the scripts and assets the reports read from projectDir,
    // which are not part of the task hash: a change reruns the reports on -resume.
    def git = { List a ->
        try { (['git', '-C', projectDir.toString()] + a).execute().text.trim() } catch (e) { '' }
    }
    def commit = workflow.commitId ? workflow.commitId.take(7) : git(['rev-parse', '--short', 'HEAD'])
    def dirty  = !workflow.commitId && git(['status', '--porcelain', '--untracked-files=no'])
    def code_files = (files("${projectDir}/bin/*.{py,R}") + files("${projectDir}/assets/**/*.{txt,tsv}"))
        .sort { it.toString() }
    def base = { p -> p ? file(p.toString()).name : null }
    def sarek_release = file("${sarek_out}/pipeline_info/sarek_release.txt")   // written by run_full.sh
    def report_info = groovy.json.JsonOutput.toJson([
        pipeline:  "${workflow.manifest.name} ${workflow.manifest.version}".toString(),
        commit:    commit ? (dirty ? "${commit} + local changes".toString() : commit) : 'unknown',
        sarek:     params.sarek_version ?: (sarek_release.exists() ? sarek_release.text.trim() : ''),
        nextflow:  nextflow.version.toString(),
        resources: [
            'gnomAD':                 base(params.gnomad_tsv),
            'ClinVar':                base(params.clinvar_tsv),
            'ClinVar protein index':  base(params.clinvar_protein_index),
            'REVEL':                  base(params.revel_tsv),
            'CIViC (civicpy cache)':  base(params.civic_cache),
            'snpEff database (MANE)': params.snpeff_mane_db ?: null,
            'MANE exons':             base(params.gcnv_mane_txt),
            'MANE CDS (TMB)':         base(params.mane_cds_bed),
            'Capture BED':            base(params.panel_bed),
            'gCNV ploidy model':      params.skip_germline_cnv ? null : base(params.gcnv_ploidy_model),
            'gCNV model (cohort)':    params.skip_germline_cnv ? null : base(params.gcnv_pon_dir),
            'gCNV intervals':         params.skip_germline_cnv ? null : base(params.gcnv_interval_list),
            'Somatic CNV PoN (GATK)': params.skip_somatic_cnv_gatk ? null : base(params.somatic_pon_hdf5),
            'COSMIC signatures':      'COSMIC v3.4 (SigProfilerAssignment)',
        ],
        code_md5:  code_files.collect { it.text }.join('\n').md5(),
    ])
    def sarek_versions = file("${sarek_out}/pipeline_info/nf_core_sarek_software_mqc_versions.yml")
    if (!sarek_versions.exists()) sarek_versions = file('NO_SAREK_VERSIONS')

    // Analyses switched off in this run: the reports say "not analysed"
    // instead of "no result" (bin/report_shared.py)
    def not_run = [
        germline_snv:     params.skip_germline_snv,
        pharmaco:         !germline_ran || params.skip_pharmaco,
        germline_cnv:     params.skip_germline_cnv,
        somatic_snv:      !somatic_ran,
        somatic_cnv:      !somatic_ran,
        somatic_cnv_gatk: params.skip_somatic_cnv_gatk,
        tmb:              !somatic_ran || params.skip_tmb,
        signatures:       !somatic_ran || params.skip_signatures,
        hrd:              !somatic_ran || params.skip_hrd,
        coverage:         params.skip_coverage,
    ].findAll { it.value }.keySet().join(',')

    // ── REPORT MODEL → laboratory HTML + clinical PDF ────────────────────────
    // One model per pair holds every value and interpretation; the two reports
    // only render it (bin/report_model.py), so they cannot diverge.
    if (!params.skip_report || !params.skip_pdf_report) {
        REPORT_MODEL(
            ch_report_input,
            params.run_name,
            not_run,
            report_info,
            sarek_versions,
            params.samplesheet ? file(params.samplesheet, checkIfExists: true) : file('NO_SAMPLESHEET')
        )
        // the renderers also need the signature directories (plots)
        def ch_render_input = REPORT_MODEL.out.model
            .join(ch_report_input.map { row -> [row[0], row[-2], row[-1]] })
        if (!params.skip_report)     FINAL_REPORT(ch_render_input)
        if (!params.skip_pdf_report) CLINICAL_PDF_REPORT(ch_render_input)
    }
}
