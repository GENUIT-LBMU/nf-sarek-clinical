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
include { FINAL_REPORT          } from '../modules/local/final_report/main'
include { CLINICAL_PDF_REPORT   } from '../modules/local/clinical_pdf_report/main'

workflow SAREK_CLINICAL {

    // ── Reference file channels ───────────────────────────────────────────────
    ch_dbnsfp        = params.dbnsfp      ? file(params.dbnsfp)               : error("--dbnsfp required")
    ch_dbnsfp_tbi    = params.dbnsfp      ? file("${params.dbnsfp}.tbi")      : error("--dbnsfp required")
    ch_snpsift_jar   = params.snpsift_jar ? file(params.snpsift_jar)          : error("--snpsift_jar required")
    ch_panel_bed     = params.panel_bed   ? file(params.panel_bed)   : []
    ch_gene_coords   = params.gene_coords_bed ? file(params.gene_coords_bed) : []

    ch_germline_panel    = file(params.germline_panel)
    ch_somatic_cnv_panel = file(params.somatic_cnv_panel)
    ch_somatic_snv_panel = file(params.somatic_snv_panel)
    ch_pharmaco_variants = file(params.pharmaco_variants)

    sarek_out = file(params.sarek_outdir)

    // ── Auto-discover Sarek outputs ───────────────────────────────────────────

    // Germline HaplotypeCaller VCFs (snpEff annotated)
    ch_hc_vcfs = Channel.fromPath(
        "${sarek_out}/annotation/haplotypecaller/**/*.vcf.gz"
    ).map { vcf ->
        def sample = vcf.parent.name
        [sample, vcf]
    }

    // Somatic Mutect2 VCFs (snpEff annotated)
    ch_mutect2_vcfs = Channel.fromPath(
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
    ch_strelka2_snvs = Channel.fromPath(
        "${sarek_out}/variant_calling/strelka/**/*.somatic_snvs.vcf.gz",
        checkIfExists: false
    ).map { vcf -> [vcf.parent.name, vcf] }

    ch_strelka2_indels = Channel.fromPath(
        "${sarek_out}/variant_calling/strelka/**/*.somatic_indels.vcf.gz",
        checkIfExists: false
    ).map { vcf -> [vcf.parent.name, vcf] }

    // CNVkit .somatic.call.cns (fallback to .call.cns)
    ch_cns = Channel.fromPath(
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
    ch_ascat_metrics  = Channel.fromPath(
        "${sarek_out}/variant_calling/ascat/**/*.metrics.txt"
    ).map { f -> [f.parent.name, f] }

    ch_ascat_segments = Channel.fromPath(
        "${sarek_out}/variant_calling/ascat/**/*segment*.txt"
    ).map { f -> [f.parent.name, f] }

    ch_ascat = ch_ascat_metrics.join(ch_ascat_segments, remainder: true)
        .map { pair, met, seg -> [pair, seg ?: file('NO_FILE'), met ?: file('NO_FILE')] }

    // Tumor sample names (for mutsig tumor column detection)
    ch_tumor_samples = ch_cns.map { pair, tumor, normal, cns -> [pair, tumor] }

    // mosdepth regions BED per sample
    ch_mosdepth = Channel.fromPath(
        "${sarek_out}/reports/mosdepth/**/*.recal.regions.bed.gz"
    ).map { f ->
        def sample = f.parent.name
        [sample, f]
    }

    // Sarek recalibrated BAMs (needed for GATK CNV modules)
    ch_bams = Channel.fromPath(
        "${sarek_out}/preprocessing/recalibrated/**/*.recal.bam"
    ).map { bam ->
        def sample = bam.parent.name
        def bai    = bam.parent.resolve("${sample}.recal.bam.bai")
        [sample, bam, bai]
    }

    // ── GERMLINE POSTPROCESSING (SNV + Pharmaco) ──────────────────────────────
    if (!params.skip_germline_snv || !params.skip_pharmaco) {
        GERMLINE_POSTPROCESS(
            ch_hc_vcfs,
            ch_dbnsfp,
            ch_dbnsfp_tbi,
            ch_snpsift_jar,
            ch_germline_panel,
            ch_pharmaco_variants
        )
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

        def gcnv_fasta = file(params.gcnv_fasta)
        def gcnv_fai   = file("${params.gcnv_fasta}.fai")

        // Only process normal samples for germline CNV
        ch_normal_bams = ch_bams.filter { sample, bam, bai -> sample ==~ /.*_N\d+$/ }

        CNV_GERMLINE_GATK(
            ch_normal_bams,
            file(params.gcnv_interval_list),
            file(params.gcnv_ploidy_model),
            file(params.gcnv_pon_dir),
            gcnv_fasta,
            gcnv_fai,
            file(params.gcnv_dict),
            file(params.gcnv_mane_txt)
        )
    }

    // ── SOMATIC POSTPROCESSING (SNV + Consensus + ASCAT/CNVkit + TMB + Signatures + HRD) ─
    // Reference FASTA for bcftools norm in consensus module
    def consensus_fasta     = params.gcnv_fasta ? file(params.gcnv_fasta) : file('NO_FILE')
    def consensus_fasta_fai = params.gcnv_fasta ? file("${params.gcnv_fasta}.fai") : file('NO_FILE')

    ch_hrd_tsv = Channel.empty()
    if (!params.skip_somatic_snv || !params.skip_somatic_cnv || !params.skip_tmb || !params.skip_signatures || !params.skip_hrd) {
        SOMATIC_POSTPROCESS(
            ch_mutect2_vcfs,
            ch_strelka2_snvs,
            ch_strelka2_indels,
            ch_cns,
            ch_ascat,
            ch_dbnsfp,
            ch_dbnsfp_tbi,
            ch_snpsift_jar,
            ch_somatic_snv_panel,
            ch_somatic_cnv_panel,
            ch_germline_panel,
            ch_gene_coords,
            ch_panel_bed,
            ch_tumor_samples,
            consensus_fasta,
            consensus_fasta_fai
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

        CNV_SOMATIC_GATK(
            ch_bams,
            file(params.somatic_interval_list),
            file(params.somatic_pon_hdf5),
            file(params.common_snps_vcf),
            snps_tbi,
            somatic_fasta,
            somatic_fai,
            file(params.somatic_dict),
            file(params.gcnv_mane_txt)
        )
    }

    // ── GENE COVERAGE ─────────────────────────────────────────────────────────
    if (!params.skip_coverage) {
        GENE_COVERAGE(
            ch_mosdepth,
            ch_germline_panel
        )
    }

    // ── ALL-UPSTREAM-DONE GATE ────────────────────────────────────────────────
    // Ensures FINAL_REPORT and CLINICAL_PDF_REPORT run only after every
    // upstream process (germline SNVs, somatic SNVs, signatures, GATK CNVs,
    // gene coverage) has finished writing its output files.
    // Each active subworkflow contributes a trigger element; collect() waits
    // for all of them, then maps to a single 'ready' value for combine().
    def ch_gate = Channel.of('base')

    if (!params.skip_germline_snv || !params.skip_pharmaco) {
        ch_gate = ch_gate.mix(GERMLINE_POSTPROCESS.out.snv_tsv.map { 'g' })
    }
    if (!params.skip_somatic_snv || !params.skip_somatic_cnv || !params.skip_tmb || !params.skip_signatures || !params.skip_hrd) {
        ch_gate = ch_gate.mix(SOMATIC_POSTPROCESS.out.somatic_snv_tsv.map { 's' })
        if (!params.skip_signatures) {
            ch_gate = ch_gate.mix(SOMATIC_POSTPROCESS.out.signatures_dir.map { 'sig' })
        }
    }
    if (!params.skip_germline_cnv) {
        // Use annotated_tsv (not segments_vcf) so the gate fires AFTER
        // MANE_ANNOTATE_GERMLINE_CNV publishes its output to disk.
        // Using segments_vcf caused a race condition where FINAL_REPORT
        // started before the annotated TSV was available.
        ch_gate = ch_gate.mix(CNV_GERMLINE_GATK.out.annotated_tsv.map { 'gcnv' })
    }
    if (!params.skip_somatic_cnv_gatk) {
        ch_gate = ch_gate.mix(CNV_SOMATIC_GATK.out.annotated_tsv.map { 'scnv' })
    }
    if (!params.skip_coverage) {
        ch_gate = ch_gate.mix(GENE_COVERAGE.out.tsv.map { 'cov' })
    }

    ch_all_done = ch_gate.collect().map { 'ready' }

    // ── FINAL REPORT ─────────────────────────────────────────────────────────
    if (!params.skip_report) {
        // Join ascat_metrics with hrd_tsv (already synced with HRD process).
        // Then combine with ch_all_done to enforce "run last" ordering.
        // Collect germline/somatic TSVs as proper Nextflow inputs so cache
        // is invalidated whenever their content changes.
        ch_germ_tsv   = params.skip_germline_snv ? Channel.value(file('NO_FILE'))
                                                  : GERMLINE_POSTPROCESS.out.snv_tsv.map { _s, f -> f }
        ch_som_tsv    = params.skip_somatic_snv  ? Channel.value(file('NO_FILE'))
                                                  : SOMATIC_POSTPROCESS.out.somatic_snv_tsv.map { _s, f -> f }
        ch_civic_tsv  = params.skip_somatic_snv  ? Channel.value(file('NO_FILE'))
                                                  : SOMATIC_POSTPROCESS.out.civic_tsv.map { _s, f -> f }

        ch_report_input = ch_ascat_metrics
            .join(ch_hrd_tsv, remainder: true)
            .map { pair, met, hrd -> [pair, met, hrd ?: file('NO_FILE')] }
            .combine(ch_all_done)
            .map { pair, met, hrd, _ -> [pair, met, hrd] }
            .combine(ch_germ_tsv)
            .combine(ch_som_tsv)
            .map { pair, met, hrd, gtsv, stsv -> [pair, met, hrd, gtsv, stsv] }

        FINAL_REPORT(
            ch_report_input,
            params.sarek_outdir,
            params.run_name
        )
    }

    // ── CLINICAL PDF REPORT ───────────────────────────────────────────────────
    if (!params.skip_pdf_report) {
        ch_pdf_input = ch_ascat_metrics
            .join(ch_hrd_tsv, remainder: true)
            .map { pair, met, hrd -> [pair, met, hrd ?: file('NO_FILE')] }
            .combine(ch_all_done)
            .map { pair, met, hrd, _ -> [pair, met, hrd] }
            .combine(ch_germ_tsv)
            .combine(ch_som_tsv)
            .map { pair, met, hrd, gtsv, stsv -> [pair, met, hrd, gtsv, stsv] }

        CLINICAL_PDF_REPORT(
            ch_pdf_input,
            ch_civic_tsv,
            params.sarek_outdir,
            params.run_name
        )
    }
}
