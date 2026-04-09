/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MUTATIONAL SIGNATURES MODULE (SigProfiler)
    ============================================
    Input:  Mutect2 filtered VCFs (PASS, tumor sample, AD_alt >= min)
    Steps:
      1. bcftools filter: PASS + tumor column + AD_alt >= alt_min
      2. SigProfilerMatrixGenerator  → SBS/DBS/ID profiles + per-sample plots
      3. SigProfilerAssignment       → COSMIC v3.4 decomposition (SBS-96, DBS, ID)
    Output: matrix profiles, assignment results, plots (PDF/PNG)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process MUTSIG_FILTER {
    tag "$pair"
    label 'process_low'
    label 'bcftools'

    publishDir "${params.outdir}/somatic/signatures/filtered_vcfs", mode: params.publish_dir_mode

    input:
    tuple val(pair), val(tumor_sample), path(vcf)
    val alt_min

    output:
    tuple val(pair), path("${pair}.sigfilter.vcf"), emit: vcf

    stub:
    """
    printf "##fileformat=VCFv4.2\\n#CHROM\\tPOS\\tID\\tREF\\tALT\\tQUAL\\tFILTER\\tINFO\\n" > ${pair}.sigfilter.vcf
    """

    script:
    """
    # Dynamically resolve tumor sample name from VCF header
    # (Mutect2 T/N VCFs list: normal first, tumor last)
    TUMOR_IN_VCF=\$(bcftools query -l ${vcf} | grep -i '_T' | head -1 || bcftools query -l ${vcf} | tail -1)

    bcftools view -f PASS,haplotype -s \${TUMOR_IN_VCF} ${vcf} \\
        | bcftools filter -i "FORMAT/AD[0:1] >= ${alt_min}" \\
        -Oz -o ${pair}.sigfilter.vcf.gz

    tabix -f -p vcf ${pair}.sigfilter.vcf.gz
    gunzip -c ${pair}.sigfilter.vcf.gz > ${pair}.sigfilter.vcf

    n=\$(grep -vc '^#' ${pair}.sigfilter.vcf || echo 0)
    echo "[MUTSIG] ${pair}: \${n} variants after filtering (AD_alt >= ${alt_min})"
    """
}


process MUTSIG_SIGPROFILER {
    tag "$run_name"
    label 'process_high'
    label 'sigprofiler'

    publishDir "${params.outdir}/somatic/signatures", mode: params.publish_dir_mode

    input:
    val run_name
    path vcf_dir    // directory containing all *.sigfilter.vcf files
    val genome_build

    output:
    path "matrix_profiles/", emit: profiles
    path "assignment/",       emit: assignment

    stub:
    """
    mkdir -p matrix_profiles assignment
    touch matrix_profiles/${run_name}_SBS96.txt
    touch assignment/${run_name}_SBS96_Activities.txt
    """

    script:
    """
    mkdir -p matrix_profiles assignment

    python3 - "${vcf_dir}" assignment matrix_profiles "${genome_build}" "${run_name}" <<'PYEOF'
import os, shutil, sys
from SigProfilerAssignment import Analyzer as Analyze
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from SigProfilerMatrixGenerator import install as spm_install

# Ensure the genome is installed (no-op if already present)
try:
    spm_install.install(sys.argv[4] if len(sys.argv) > 4 else "GRCh38", rsync=False, bash=True)
except Exception as _e:
    print(f"Note: genome install returned: {_e}", file=sys.stderr)

vcf_dir      = sys.argv[1]
assign_dir   = sys.argv[2]
profiles_dir = sys.argv[3]
genome       = sys.argv[4]   # GRCh37 | GRCh38
proj         = sys.argv[5]

# A) MatrixGenerator: SBS/DBS/ID profiles per sample
# NOTE: SigProfilerMatrixGenerator writes to {vcf_dir}/{mg_project}/output/
#       (NOT to profiles_dir). We run it, then copy the output to profiles_dir.
mg_project  = f"{proj}_MG"
os.makedirs(profiles_dir, exist_ok=True)

matrices = matGen.SigProfilerMatrixGeneratorFunc(
    mg_project, genome, vcf_dir,
    exome=True, plot=True, chrom_based=False
)

# Copy MatrixGenerator output from where it actually lands to profiles_dir.
# SigProfilerMatrixGenerator writes to {vcf_dir}/output/ (no project subdir).
mg_src = os.path.join(vcf_dir, "output")
if os.path.exists(mg_src):
    mg_dst = os.path.join(profiles_dir, "output")
    if os.path.exists(mg_dst):
        shutil.rmtree(mg_dst)
    shutil.copytree(mg_src, mg_dst)
    print(f"MatrixGenerator output copied → {mg_dst}")
else:
    print(f"WARNING: MatrixGenerator output not found at {mg_src}", file=sys.stderr)

# B) SigProfilerAssignment for SBS-96, DBS, ID
os.makedirs(assign_dir, exist_ok=True)

def run_fit(context_type, tag, collapse_sbs96=True):
    out = os.path.join(assign_dir, f"{proj}_{tag}")
    if os.path.exists(out):
        shutil.rmtree(out)
    try:
        Analyze.cosmic_fit(
            samples=vcf_dir,
            output=out,
            input_type="vcf",
            context_type=context_type,
            genome_build=genome,
            cosmic_version=3.4,
            exome=True,
            make_plots=True,
            sample_reconstruction_plots="pdf",
            collapse_to_SBS96=collapse_sbs96,
            verbose=False
        )
    except Exception as e:
        print(f"WARNING: {tag} assignment failed: {e}", file=sys.stderr)
    return out

run_fit("96",    "SBS96", collapse_sbs96=True)
run_fit("DINUC", "DBS",   collapse_sbs96=False)
run_fit("ID",    "ID",    collapse_sbs96=False)
print("SigProfilerAssignment done →", assign_dir)
PYEOF

    # Resolve sample name from the Activities.txt (first non-header data row, col 0)
    SBS96_ACT="\$(find assignment -name "Assignment_Solution_Activities.txt" -path "*SBS96*" 2>/dev/null | head -1)"
    if [ -n "\${SBS96_ACT}" ]; then
        SAMPLE_NAME="\$(awk 'NR==2{print \$1; exit}' "\${SBS96_ACT}")"
    else
        SAMPLE_NAME="${run_name}"
    fi

    SBS96_DIR="\$(find assignment -mindepth 1 -maxdepth 1 -type d -name "*SBS96*" 2>/dev/null | head -1 || true)"
    ID_DIR="\$(find assignment -mindepth 1 -maxdepth 1 -type d -name "*_ID" 2>/dev/null | head -1 || true)"

    MATRIX_OUT="\$([ -d matrix_profiles/output ] && echo "matrix_profiles/output" || echo "")"

    generate_signature_plots.py \\
        \${SBS96_DIR:+--sbs96_dir "\${SBS96_DIR}"} \\
        \${ID_DIR:+--id_dir "\${ID_DIR}"} \\
        --sample "\${SAMPLE_NAME}" \\
        \${MATRIX_OUT:+--matrix_dir "\${MATRIX_OUT}"} \\
        || echo "WARNING: generate_signature_plots.py returned non-zero" >&2
    """
}
