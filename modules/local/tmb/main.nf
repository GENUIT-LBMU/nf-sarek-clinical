/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    TUMOR MUTATIONAL BURDEN (TMB) MODULE
    =====================================
    Input:  snpEff-annotated Mutect2 VCFs from Sarek
    TMB = non-synonymous somatic variants (HIGH / MODERATE: missense, nonsense,
    frameshift, in-frame, splice, start/stop) per Mb of coding territory.
    Steps:
      1. SnpSift filter: PASS + HIGH/MODERATE + DP >= min_dp + AF >= min_af
         + AD_alt >= filters.tmb.min_alt_reads (synonymous variants are not counted)
      2. Count the variants inside capture BED ∩ MANE Select CDS / its size in Mb
    Output: tmb_per_sample.tsv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

process TMB_FILTER {
    tag "$pair"
    label 'process_low'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/tmb/filtered", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path snpsift_jar
    val min_dp
    val min_af
    val min_alt_reads

    output:
    tuple val(pair), path("${pair}.tmb_filter.vcf"), emit: vcf

    stub:
    """
    printf "##fileformat=VCFv4.2\\n" > ${pair}.tmb_filter.vcf
    """

    script:
    """
    FILTER_EXPR="( FILTER = 'PASS' | FILTER = 'haplotype' ) \\
        & (ANN[*].IMPACT = 'HIGH' | ANN[*].IMPACT = 'MODERATE') \\
        & ( GEN[1].DP >= ${min_dp} ) \\
        & ( GEN[1].AF >= ${min_af} ) \\
        & ( GEN[1].AD[1] >= ${min_alt_reads} )"

    java -jar ${snpsift_jar} filter "\${FILTER_EXPR}" ${vcf} \\
        > ${pair}.tmb_filter.vcf
    """
}


process TMB_CALCULATE {
    tag "$pair"
    label 'process_single'
    label 'snpsift'

    publishDir "${params.outdir}/somatic/tmb", mode: params.publish_dir_mode

    input:
    tuple val(pair), path(vcf)
    path panel_bed  // capture BED; pass [] if not available
    path cds_bed    // MANE Select CDS BED (--mane_cds_bed); pass [] if not available

    output:
    tuple val(pair), path("${pair}.tmb.tsv"), emit: tsv

    stub:
    """
    printf "sample\\tn_variants\\tpanel_bases\\tpanel_mb\\tTMB_mut_per_mb\\tterritory\\n" > ${pair}.tmb.tsv
    printf "${pair}\\t0\\tNA\\tNA\\tNA\\tNA\\n" >> ${pair}.tmb.tsv
    """

    script:
    // TMB = non-synonymous somatic variants (HIGH / MODERATE, from TMB_FILTER)
    // in the coding territory / its size in Mb. Territory = capture BED ∩ MANE
    // Select CDS, computed here for each run; without the CDS BED the whole
    // capture is used (territory "capture", logged).
    def bed = panel_bed && panel_bed.name != 'NO_FILE' ? panel_bed : ''
    def cds = cds_bed && cds_bed.name != 'NO_FILE' ? cds_bed : ''
    """
    python3 - "${vcf}" "${bed}" "${cds}" "${pair}" <<'PYEOF'
import bisect, collections, sys

vcf, bed, cds, pair = sys.argv[1:5]

def load(path):
    iv = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            f = line.split('\\t')
            iv[f[0]].append((int(f[1]), int(f[2])))
    merged = {}
    for c, v in iv.items():
        m = []
        for s, e in sorted(v):
            if m and s <= m[-1][1]:
                m[-1][1] = max(m[-1][1], e)
            else:
                m.append([s, e])
        merged[c] = m
    return merged

def intersect(a, b):
    out = {}
    for c in a.keys() & b.keys():
        x, y, i, j, m = a[c], b[c], 0, 0, []
        while i < len(x) and j < len(y):
            s, e = max(x[i][0], y[j][0]), min(x[i][1], y[j][1])
            if s < e:
                m.append([s, e])
            if x[i][1] < y[j][1]:
                i += 1
            else:
                j += 1
        out[c] = m
    return out

if bed and cds:
    region, territory = intersect(load(bed), load(cds)), 'capture_x_MANE_CDS'
elif bed:
    region, territory = load(bed), 'capture'
    sys.stderr.write('WARNING: --mane_cds_bed not set: TMB denominator is the whole capture BED (UTR and flanks included)\\n')
else:
    region, territory = None, 'NA'
starts = {c: [s for s, _e in v] for c, v in (region or {}).items()}

def inside(c, pos):
    v = region.get(c)
    if not v:
        return False
    k = bisect.bisect_right(starts[c], pos - 1) - 1
    return k >= 0 and v[k][0] < pos <= v[k][1]

n = 0
with open(vcf) as fh:
    for line in fh:
        if line.startswith('#'):
            continue
        f = line.split('\\t', 6)
        if region is None or inside(f[0], int(f[1])):
            n += len(f[4].split(','))          # one count per ALT allele
bases = sum(e - s for v in (region or {}).values() for s, e in v) if region is not None else None
mb = bases / 1e6 if bases else None
tmb = '%.4f' % (n / mb) if mb else 'NA'
with open(pair + '.tmb.tsv', 'w') as out:
    out.write('sample\\tn_variants\\tpanel_bases\\tpanel_mb\\tTMB_mut_per_mb\\tterritory\\n')
    out.write('\\t'.join([pair, str(n), str(bases) if bases else 'NA',
                          '%.6f' % mb if mb else 'NA', tmb, territory]) + '\\n')
sys.stderr.write('[TMB] %s: n=%d  territory=%s  Mb=%s  TMB=%s\\n' % (pair, n, territory, '%.2f' % mb if mb else 'NA', tmb))
PYEOF
    """
}
