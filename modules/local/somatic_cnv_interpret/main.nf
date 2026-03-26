/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SOMATIC CNV INTERPRETATION MODULE (ASCAT + CNVkit)
    ====================================================
    Uses parse_somatic_cnv.py (in bin/) to cross-reference:
      - CNVkit .somatic.call.cns → HOMDEL, HETDEL, AMP events (38-gene somatic panel)
      - ASCAT segments + metrics  → LOH 2nd hit events (90-gene germline panel)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process SOMATIC_CNV_INTERPRET {
    tag "$pair"
    label 'process_low'
    label 'python'

    publishDir "${params.outdir}/somatic/cnv/ascat_cnvkit", mode: params.publish_dir_mode

    input:
    tuple val(pair),
          val(tumor), val(normal),
          path(cns_file),          // CNVkit .somatic.call.cns (or .call.cns)
          path(ascat_segments),    // ASCAT *segments.txt  (may be empty string)
          path(ascat_metrics)      // ASCAT *metrics.txt   (may be empty string)
    path germline_panel
    path somatic_panel
    path gene_coords_bed           // WES BED for gene coordinates (optional)

    output:
    tuple val(pair), path("${pair}_somatic_cnv_report.tsv"), emit: tsv
    tuple val(pair), path("${pair}_somatic_cnv_report.txt"), emit: txt

    stub:
    """
    printf "pair\\tevent_type\\tgene\\tcn\\n" > ${pair}_somatic_cnv_report.tsv
    printf "Somatic CNV Report: ${pair}\\nNo events (stub run)\\n" > ${pair}_somatic_cnv_report.txt
    """

    script:
    def ascat_seg_arg  = ascat_segments  ? "--ascat_segments ${ascat_segments}" : ""
    def ascat_met_arg  = ascat_metrics   ? "--ascat_metrics ${ascat_metrics}"   : ""
    def coords_arg     = gene_coords_bed ? "--genecoords ${gene_coords_bed}"    : ""
    """
    parse_somatic_cnv.py \\
        --cns          ${cns_file} \\
        ${ascat_seg_arg} \\
        ${ascat_met_arg} \\
        ${coords_arg} \\
        --germline     ${germline_panel} \\
        --somatic_genes ${somatic_panel} \\
        --sample       ${pair} \\
        --out_tsv      ${pair}_somatic_cnv_report.tsv \\
        --out_txt      ${pair}_somatic_cnv_report.txt \\
        --amp_threshold ${params.filters.somatic_cnv.amp_threshold}
    """
}
