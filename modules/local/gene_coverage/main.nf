/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENE COVERAGE MODULE
    =====================
    Input:  mosdepth *.recal.regions.bed.gz files (from Sarek preprocessing reports)
    Output: per-gene mean/min/max coverage TSV per sample
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process GENE_COVERAGE {
    tag "$sample"
    label 'process_low'
    label 'r_base'

    publishDir "${params.outdir}/qc/gene_coverage", mode: params.publish_dir_mode

    input:
    tuple val(sample), path(regions_bed_gz)
    path genes_list

    output:
    tuple val(sample), path("${sample}_gene_coverage.tsv"), emit: tsv

    stub:
    """
    printf "gene\\tmean_coverage\\tmin_coverage\\tmax_coverage\\tn_regions\\n" > ${sample}_gene_coverage.tsv
    printf "BRCA1\\t120.5\\t45\\t350\\t23\\n" >> ${sample}_gene_coverage.tsv
    """

    script:
    """
    # Decompress mosdepth regions BED
    gunzip -c ${regions_bed_gz} > ${sample}.recal.regions.bed

    Rscript - <<'REOF'
    library(data.table)

    bed_file   <- "${sample}.recal.regions.bed"
    genes_file <- "${genes_list}"
    out_file   <- "${sample}_gene_coverage.tsv"

    # Load mosdepth BED: chrom, start, end, name (gene), coverage
    bed <- fread(bed_file, header=FALSE, sep="\\t",
                 col.names=c("chrom","start","end","gene","coverage"))

    # Load gene list
    genes_of_interest <- readLines(genes_file)
    genes_of_interest <- trimws(genes_of_interest[!grepl("^#", genes_of_interest)])
    genes_of_interest <- toupper(genes_of_interest[nchar(genes_of_interest) > 0])

    bed[["gene"]] <- toupper(bed[["gene"]])

    # Filter and summarise
    bed_filt <- bed[gene %in% genes_of_interest]

    if (nrow(bed_filt) == 0) {
      message("WARNING: No regions matched genes of interest in ", bed_file)
      result <- data.table(gene=character(), mean_coverage=numeric(),
                           min_coverage=numeric(), max_coverage=numeric(), n_regions=integer())
    } else {
      result <- bed_filt[, .(
        mean_coverage = round(mean(coverage, na.rm=TRUE), 2),
        min_coverage  = min(coverage, na.rm=TRUE),
        max_coverage  = max(coverage, na.rm=TRUE),
        n_regions     = .N
      ), by=gene][order(gene)]
    }

    fwrite(result, out_file, sep="\\t")
    message("Gene coverage written to ", out_file, " (", nrow(result), " genes)")
    REOF
    """
}
