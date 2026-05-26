/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HRD PROFILER MODULE
    ====================
    Runs HRProfiler for Homologous Recombination Deficiency (HRD) analysis.
    Uses filtered somatic SNVs (VCF) + ASCAT copy-number segments.

    Reference: Chopra et al., 2024 (github.com/FlorPio/HRProfiler-workflow)
    Docker:    florpio/hrprofiler:1.0-full

    Input CNV format conversion:
      Sarek ASCAT output: sample  chr  startpos  endpos  nMajor  nMinor
      HRProfiler ASCAT:   Sample  chr  startpos  endpos  total.copy.number.inTumour  nMajor  nMinor
      (chr column is also prefixed with "chr" if not already present)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process HRD_PROFILER {
    tag "$pair"
    label 'process_medium'
    label 'hrprofiler'

    publishDir "${params.outdir}/somatic/hrd", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(filtered_vcf), path(ascat_segs)
    val   organ
    val   genome

    output:
    tuple val(pair), path("${pair}_hrd_predictions.tsv"), emit: predictions
    path "hrd_results/",                                  emit: results_dir, optional: true

    stub:
    """
    echo -e "sample\thrd.prob\tprediction" > ${pair}_hrd_predictions.tsv
    echo -e "${pair}\t0.42\tHRD-"          >> ${pair}_hrd_predictions.tsv
    mkdir -p hrd_results
    """

    script:
    """
    # HRProfiler concatenates RESULT_DIR + "output" (no os.path.join) so the dir
    # must end with "/" and the output/ subdir must exist before the run.
    mkdir -p snv_input cnv_input hrd_results/output hrd_results/logs

    # Stage VCF – use pair name as sample ID (SigProfilerMatrixGenerator derives
    # the sample name from the filename without extension)
    cp ${filtered_vcf} snv_input/${pair}.vcf

    # Convert ASCAT segments to HRProfiler ASCAT format
    # Sarek output header: sample  chr  startpos  endpos  nMajor  nMinor
    # HRProfiler expects:  Sample  chr  startpos  endpos  total.copy.number.inTumour  nMajor  nMinor
    # Also adds "chr" prefix when missing (sarek uses numeric chromosome IDs)
    awk '
    NR==1 {
        print "Sample\\tchr\\tstartpos\\tendpos\\ttotal.copy.number.inTumour\\tnMajor\\tnMinor"
        next
    }
    {
        chr   = (\$2 ~ /^chr/) ? \$2 : "chr"\$2
        total = \$5 + \$6
        print "${pair}\\t" chr "\\t" \$3 "\\t" \$4 "\\t" total "\\t" \$5 "\\t" \$6
    }' ${ascat_segs} > cnv_input/${pair}.segments.txt

    lines=\$(wc -l < cnv_input/${pair}.segments.txt)
    echo "[HRD] ${pair}: \${lines} segment lines (including header)"

    # Guard: skip HRProfiler if VCF has too few somatic variants. HRProfiler
    # requires several SNVs to build the 96-context SBS matrix used by its
    # classifier; with <5 the internal model crashes silently (exit 1) without
    # a traceback. Emit NA stub and skip the call instead of aborting.
    n_vars=\$(grep -cv '^#' snv_input/${pair}.vcf || true)
    echo "[HRD] ${pair}: \${n_vars} somatic variants in VCF"
    if [[ "\${n_vars}" -lt 5 ]]; then
        echo "[HRD] WARNING: \${n_vars} somatic variants (<5 required); emitting NA stub" >&2
        echo -e "sample\\thrd.prob\\tprediction" > ${pair}_hrd_predictions.tsv
        echo -e "${pair}\\tNA\\tNA"              >> ${pair}_hrd_predictions.tsv
        exit 0
    fi

    # Run HRProfiler (HRD.py is in bin/ and staged into PATH by Nextflow)
    # Note: trailing slash on --output-dir is required — HRProfiler does
    # RESULT_DIR + "output" string concatenation (not os.path.join)
    # `|| true` so a silent crash falls through to the NA-stub fallback below
    # instead of aborting the whole workflow.
    HRD.py \\
        --snv-dir    snv_input \\
        --cnv-dir    cnv_input \\
        --output-dir hrd_results/ \\
        --organ      ${organ} \\
        --genome     ${genome} \\
        --cnv-file-type ASCAT \\
        --hrd-threshold ${params.hrd_threshold} || true

    # Locate predictions output (filename varies by HRProfiler version)
    pred="\$(find hrd_results -maxdepth 4 -type f \\
        \\( -name "*hrd_predictions*" -o -name "*HRD_predictions*" -o -name "*prediction*" \\) \\
        2>/dev/null | grep -iv '\\.png' | head -1)"

    if [[ -n "\$pred" ]]; then
        cp "\$pred" ${pair}_hrd_predictions.tsv
        echo "[HRD] ${pair}: predictions written to ${pair}_hrd_predictions.tsv"
    else
        echo "[HRD] WARNING: no predictions file found in hrd_results/; emitting stub" >&2
        echo -e "sample\\thrd.prob\\tprediction" > ${pair}_hrd_predictions.tsv
        echo -e "${pair}\\tNA\\tNA"              >> ${pair}_hrd_predictions.tsv
    fi
    """
}
