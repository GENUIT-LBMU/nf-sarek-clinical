"""
cnv_shared.py — harmonised somatic CNV events (CNVkit+ASCAT and GATK) for the
HTML report and the PDF report. Single source: edit only here.

  CNVkit+ASCAT (parse_somatic_cnv.py TSV): HOMDEL, AMP (CN >= amp_threshold) in
    the somatic CNV panel; 2nd_hit_LOH in the germline panel.
  GATK (annotate_cnvs_complete_fixed.R TSV), EVENT_TYPE:
    DELETION with CN 0       -> HOMDEL
    AMPLIFICATION, AMP+LOH   -> AMP when CN >= amp_threshold (AMP+LOH is the
                                common focal amplification of one allele)
    CN-LOH                   -> LOH cnLOH
    DEL+LOH                  -> LOH HETDEL+LOH
  AMP+LOH is not used as LOH: a one-copy gain (CN 3, 2:1) has an expected
  minor allele fraction ~0.33, below the R threshold 0.35, with both alleles
  present (IMAX1278: 14 such genes, none LOH in ASCAT). Gain-LOH comes from
  ASCAT (LOH_major_amp), which estimates the minor copy number.
  Confidence HIGH when both sources report the gene, else MODERATE.
"""
import csv
import os


def load_panel(path):
    genes = set()
    if path and os.path.isfile(path):
        with open(path) as fh:
            for line in fh:
                g = line.strip()
                if g and not g.startswith('#'):
                    genes.add(g.upper())
    return genes


def _cn(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return -1


def _read_ascat_cnvkit(path, somatic_panel, germline_panel, amp_threshold):
    homdel, amp, loh = {}, {}, {}
    if not path or not os.path.isfile(path):
        return homdel, amp, loh
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene  = (row.get('gene') or '').strip().upper()
            etype = (row.get('event_type') or '').strip()
            cn    = _cn(row.get('cn', -1))
            chrom = row.get('chrom', '')
            if etype == 'HOMDEL' and gene in somatic_panel:
                homdel[gene] = {'cn': 0, 'chrom': chrom, 'log2': row.get('log2', '')}
            elif etype == 'AMP' and gene in somatic_panel and cn >= amp_threshold:
                if gene not in amp or cn > amp[gene]['cn']:
                    amp[gene] = {'cn': cn, 'chrom': chrom, 'log2': row.get('log2', '')}
            elif etype == '2nd_hit_LOH' and gene in germline_panel:
                loh[gene] = {'cn': cn, 'loh_type': (row.get('loh_type') or '').strip(), 'chrom': chrom}
    return homdel, amp, loh


def _read_gatk(path, somatic_panel, germline_panel, amp_threshold):
    homdel, amp, loh = {}, {}, {}
    if not path or not os.path.isfile(path):
        return homdel, amp, loh
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            gene  = (row.get('GENE') or '').strip().upper()
            etype = (row.get('EVENT_TYPE') or '').strip()
            cn    = _cn(row.get('CN_ABSOLUTE', -1))
            chrom = row.get('CONTIG', '')
            if etype == 'DELETION' and cn == 0 and gene in somatic_panel:
                homdel[gene] = {'cn': 0, 'chrom': chrom}
            if etype in ('AMPLIFICATION', 'AMP+LOH') and cn >= amp_threshold and gene in somatic_panel:
                if gene not in amp or cn > amp[gene]['cn']:
                    amp[gene] = {'cn': cn, 'chrom': chrom}
            if etype in ('CN-LOH', 'DEL+LOH') and gene in germline_panel:
                loh_type = 'cnLOH' if etype == 'CN-LOH' else 'HETDEL+LOH'
                loh[gene] = {'cn': cn, 'loh_type': loh_type, 'chrom': chrom}
    return homdel, amp, loh


def harmonize_cnv(ascat_cnvkit_path, gatk_somatic_path, somatic_panel_path,
                  germline_panel_path, amp_threshold=5):
    """{'homdel': [...], 'amp': [...], 'loh': [...]}; each row has gene, cn_p
    (CNVkit+ASCAT or GATK), cn_s (GATK), chrom, tools, conf; loh rows loh_type."""
    somatic_panel  = load_panel(somatic_panel_path)
    germline_panel = load_panel(germline_panel_path)
    a_homdel, a_amp, a_loh = _read_ascat_cnvkit(ascat_cnvkit_path, somatic_panel, germline_panel, amp_threshold)
    g_homdel, g_amp, g_loh = _read_gatk(gatk_somatic_path, somatic_panel, germline_panel, amp_threshold)

    def tools(a, g, first):
        return '+'.join(t for t in [first if a else '', 'GATK' if g else ''] if t)

    homdel_list, amp_list, loh_list = [], [], []
    for gene in sorted(set(a_homdel) | set(g_homdel)):
        a, g = a_homdel.get(gene), g_homdel.get(gene)
        homdel_list.append({'gene': gene, 'cn_p': 0, 'cn_s': g['cn'] if g else '',
                            'chrom': (a or g)['chrom'], 'tools': tools(a, g, 'CNVkit+ASCAT'),
                            'conf': 'HIGH' if (a and g) else 'MODERATE'})
    for gene in sorted(set(a_amp) | set(g_amp)):
        a, g = a_amp.get(gene), g_amp.get(gene)
        amp_list.append({'gene': gene, 'cn_p': a['cn'] if a else g['cn'], 'cn_s': g['cn'] if g else '',
                         'chrom': (a or g)['chrom'], 'tools': tools(a, g, 'CNVkit+ASCAT'),
                         'conf': 'HIGH' if (a and g) else 'MODERATE'})
    for gene in sorted(set(a_loh) | set(g_loh)):
        a, g = a_loh.get(gene), g_loh.get(gene)
        loh_list.append({'gene': gene, 'loh_type': a['loh_type'] if a else (g['loh_type'] if g else 'LOH'),
                         'cn_p': a['cn'] if a else '', 'cn_s': g['cn'] if g else '',
                         'chrom': (a or g)['chrom'], 'tools': tools(a, g, 'ASCAT'),
                         'conf': 'HIGH' if (a and g) else 'MODERATE'})
    # HETDEL+LOH first (most severe), then cnLOH, then the rest
    _order = {'HETDEL+LOH': 0, 'cnLOH': 1}
    loh_list.sort(key=lambda r: (_order.get(r['loh_type'], 2), r['gene']))
    return {'homdel': homdel_list, 'amp': amp_list, 'loh': loh_list}


def loh_relevant_genes(germline_tsv='', somatic_tsvs=()):
    """Genes whose LOH is clinically relevant: a germline ACMG P/LP variant, or a
    somatic SNV (core panel or second-hit track TSVs; AMP Tier IV polymorphisms
    excluded). The PDF lists only the LOH of these genes and the HTML marks them."""
    from acmg_shared import acmg_class, getter_for_dict
    genes = set()
    if germline_tsv and os.path.isfile(germline_tsv):
        with open(germline_tsv) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                if acmg_class(getter_for_dict(row)) in ('P', 'LP'):
                    genes.add((row.get('ANN[*].GENE') or '').strip().upper())
    for path in somatic_tsvs:
        if path and os.path.isfile(path):
            with open(path) as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    if (row.get('POP_COMMON') or '').strip() != '1':
                        genes.add((row.get('ANN[*].GENE') or '').strip().upper())
    genes.discard('')
    return genes
