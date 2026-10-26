#!/usr/bin/env bash
# Builds a fake sarek output tree and empty resource files for the -stub test
# (conf/test.config). Stub blocks only create their outputs, so the files only
# need to exist with the names and layout the workflow looks for.
#   bash tests/stub/make_fixture.sh [DIR]   (default: tests/stub/fixture)
set -euo pipefail
D="${1:-$(dirname "$0")/fixture}"
rm -rf "$D"
S="$D/sarek_results"
T=T1; N=N1; P="${T}_vs_${N}"

mk() { mkdir -p "$(dirname "$1")"; : > "$1"; }

mk "$S/annotation/haplotypecaller/$N/$N.haplotypecaller.filtered_snpEff.ann.vcf.gz"
mk "$S/annotation/mutect2/$P/$P.mutect2.filtered_snpEff.ann.vcf.gz"
mk "$S/variant_calling/strelka/$P/$P.strelka.somatic_snvs.vcf.gz"
mk "$S/variant_calling/strelka/$P/$P.strelka.somatic_indels.vcf.gz"
mk "$S/variant_calling/cnvkit/$P/$T.somatic.call.cns"
mk "$S/variant_calling/msisensorpro/$P/$P"
for s in $T $N; do
    mk "$S/reports/mosdepth/$s/$s.recal.regions.bed.gz"
    mk "$S/reports/mosdepth/$s/$s.recal.mosdepth.summary.txt"
done
mkdir -p "$S/variant_calling/ascat/$P"
printf 'purity\tploidy\n0.6\t2.0\n' > "$S/variant_calling/ascat/$P/$P.metrics.txt"
mk "$S/variant_calling/ascat/$P/$P.segments.txt"
mkdir -p "$S/pipeline_info"; echo "3.5.1" > "$S/pipeline_info/sarek_release.txt"

# sarek samplesheet with the optional identification columns (reports)
cat > "$D/samplesheet.csv" <<CSV
patient,sex,status,sample,lane,fastq_1,fastq_2,patient_name,birth_date,medical_record,ordering_physician,institution,diagnosis,accession_id,specimen_type,specimen_site,collection_date,received_date,tumor_content
P1,XX,0,$N,L1,n_1.fq.gz,n_2.fq.gz,Test Patient,1970-01-01,MR-1,Dr. Test,Test Hospital,Colon adenocarcinoma,A-2,Blood,,2026-01-01,2026-01-02,
P1,XX,1,$T,L1,t_1.fq.gz,t_2.fq.gz,Test Patient,1970-01-01,MR-1,Dr. Test,Test Hospital,Colon adenocarcinoma,A-1,FFPE,Sigmoid colon,2026-01-01,2026-01-02,60
CSV

R="$D/resources"
for f in snpSift.jar civicpy_cache.pkl genome.fa genome.fa.fai \
         gnomad_panel.tsv.gz gnomad_panel.tsv.gz.tbi gnomad_panel_regions.bed.gz gnomad_panel_regions.bed.gz.tbi \
         clinvar_panel.tsv.gz clinvar_panel.tsv.gz.tbi revel_panel.tsv.gz revel_panel.tsv.gz.tbi; do
    mk "$R/$f"
done
echo "fixture: $D"
