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

workflow SAREK_CLINICAL {

    // ── Reference file channels ───────────────────────────────────────────────
    ch_dbnsfp        = params.dbnsfp      ? file(params.dbnsfp)      : error("--dbnsfp required")
    ch_snpsift_jar   = params.snpsift_jar ? file(params.snpsift_jar) : error("--snpsift_jar required")
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

    // CNVkit .somatic.call.cns (fallback to .call.cns)
    ch_cns = Channel.fromPath(
        "${sarek_out}/variant_calling/cnvkit/**/*.{somatic.call.cns,call.cns}"
    ).map { cns ->
        def pair   = cns.parent.name
        def parts  = pair.tokenize('_vs_')
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
        if (!params.gcnv_exon_r)        error("--gcnv_exon_r is required when --skip_germline_cnv false")

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
            file(params.gcnv_exon_r),
            file(params.gcnv_mane_txt)
        )
    }

    // ── SOMATIC POSTPROCESSING (SNV + ASCAT/CNVkit + TMB + Signatures) ───────
    if (!params.skip_somatic_snv || !params.skip_somatic_cnv || !params.skip_tmb || !params.skip_signatures) {
        SOMATIC_POSTPROCESS(
            ch_mutect2_vcfs,
            ch_cns,
            ch_ascat,
            ch_dbnsfp,
            ch_snpsift_jar,
            ch_somatic_snv_panel,
            ch_somatic_cnv_panel,
            ch_germline_panel,
            ch_gene_coords,
            ch_panel_bed,
            ch_tumor_samples
        )
    }

    // ── SOMATIC CNV – GATK PAIRED WORKFLOW ───────────────────────────────────
    if (!params.skip_somatic_cnv_gatk) {
        // Validate required params
        if (!params.somatic_pon_hdf5)      error("--somatic_pon_hdf5 is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_interval_list) error("--somatic_interval_list is required when --skip_somatic_cnv_gatk false")
        if (!params.common_snps_vcf)       error("--common_snps_vcf is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_fasta)         error("--somatic_fasta is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_dict)          error("--somatic_dict is required when --skip_somatic_cnv_gatk false")
        if (!params.somatic_annotate_r)    error("--somatic_annotate_r is required when --skip_somatic_cnv_gatk false")
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
            file(params.somatic_annotate_r),
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

    // ── FINAL REPORT ─────────────────────────────────────────────────────────
    if (!params.skip_report) {
        FINAL_REPORT(
            ch_ascat_metrics,
            params.sarek_outdir,
            params.run_name
        )
    }
}
