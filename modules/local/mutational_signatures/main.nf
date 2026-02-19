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
    bcftools view -f PASS -s ${tumor_sample} ${vcf} \\
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

vcf_dir      = sys.argv[1]
assign_dir   = sys.argv[2]
profiles_dir = sys.argv[3]
genome       = sys.argv[4]   # GRCh37 | GRCh38
proj         = sys.argv[5]

# A) MatrixGenerator: SBS/DBS/ID profiles per sample
mg_project  = f"{proj}_MG"
mg_out_root = os.path.join(profiles_dir, "output", mg_project)
os.makedirs(profiles_dir, exist_ok=True)
if os.path.exists(mg_out_root):
    shutil.rmtree(mg_out_root)

matrices = matGen.SigProfilerMatrixGeneratorFunc(
    mg_project, genome, vcf_dir,
    exome=True, plot=True, chrom_based=False
)
print("MatrixGenerator done →", mg_out_root)

# B) SigProfilerAssignment for SBS-96, DBS, ID
os.makedirs(assign_dir, exist_ok=True)

def run_fit(context_type, tag):
    out = os.path.join(assign_dir, f"{proj}_{tag}")
    if os.path.exists(out):
        shutil.rmtree(out)
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
        verbose=False
    )
    return out

run_fit("96",    "SBS96")
run_fit("DINUC", "DBS")
run_fit("ID",    "ID")
print("SigProfilerAssignment done →", assign_dir)
PYEOF
    """
}
