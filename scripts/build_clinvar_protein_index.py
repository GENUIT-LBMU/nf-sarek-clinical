#!/usr/bin/env python3
"""
build_clinvar_protein_index.py — build the ClinVar protein index used for
ACMG PS1 / PM5 (SNPSIFT_EXTRACT_GERMLINE, --clinvar_protein_index).

Output (in --outdir):
  clinvar_<DATE>_protein_index.tsv.gz   one row per ClinVar P/LP missense variant
      GENE AA_POS AA_REF AA_ALT HGVS_P CLASS STARS VCV CHROM POS REF ALT SPLICE_REGION
      (MANE Select protein coordinates; SPLICE_REGION=1 at exon edges, where the
      pathogenicity may come from splicing: not used as PS1 / PM5 reference)
  clinvar_<DATE>_protein_index.manifest.txt

Input: the ClinVar panel subset annotated with the snpEff MANE database, e.g.

  # 1. panel subset, chr names, left-normalised, split (as build_clinvar_panel.sh)
  docker run --rm -u $(id -u):$(id -g) -v /data:/data -w $PWD \\
      quay.io/biocontainers/bcftools:1.20--h8b25389_0 bash -c "
      bcftools view -R regions.nochr.bed /data/refs/clinvar/clinvar_20260822.vcf.gz -Ou \\
      | bcftools annotate --rename-chrs chr_map.txt -x ^INFO/CLNSIG,INFO/CLNREVSTAT -Ou \\
      | bcftools norm -f /data/refs/hg38/hg38.fa -c x -m -any -Ou \\
      | bcftools sort -Ou | bcftools norm -d exact -Ov -o clinvar_panel.vcf"
  # 2. MANE annotation (same database as --snpeff_mane_data)
  docker run --rm -u $(id -u):$(id -g) -v /data:/data -w $PWD \\
      quay.io/biocontainers/snpeff:5.2--hdfd78af_1 snpEff -Xmx8g ann -nodownload \\
      -noStats -dataDir /data/refs/snpeff/data GRCh38.mane.1.2.refseq clinvar_panel.vcf \\
      > clinvar_panel.mane.vcf
  # 3. index
  python3 scripts/build_clinvar_protein_index.py --vcf clinvar_panel.mane.vcf \\
      --mane-plus-clinical assets/mane/mane_v1.2_plus_clinical.txt --outdir /data/refs/clinvar

The database also holds MANE Plus Clinical transcripts, which snpEff may list
first (BRAF NM_001374258.1: p.Val640Glu): the ANN of a transcript in
--mane-plus-clinical is never used, so every row is in MANE Select coordinates,
as the germline HGVS (SNPEFF_MANE lists MANE Select first).
"""
import argparse
import datetime
import gzip
import hashlib
import os
import re
import sys

HGVS_MISSENSE = re.compile(r'^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$')


def stars(revstat):
    r = revstat.lower()
    if 'practice_guideline' in r:
        return 4
    if 'expert_panel' in r:
        return 3
    if 'multiple_submitters' in r and 'no_conflicts' in r:
        return 2
    if 'single_submitter' in r or 'conflicting' in r:
        return 1
    return 0


def plp_class(sig):
    s = sig.lower()
    if 'conflicting' in s:
        return None
    if s.startswith('likely_pathogenic'):
        return 'LP'
    if s.startswith('pathogenic'):
        return 'LP' if 'likely_pathogenic' in s else 'P'
    return None


def load_plus_clinical(path):
    with open(path) as fh:
        return {l.split('\t')[0].strip() for l in fh if l.strip() and not l.startswith('#')}


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True, help='ClinVar panel subset annotated with snpEff MANE')
    ap.add_argument('--mane-plus-clinical', required=True,
                    help='MANE Plus Clinical transcripts of the snpEff MANE release (assets/mane)')
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    plus = load_plus_clinical(a.mane_plus_clinical)
    date, snpeff_cmd, rows = None, '', []
    opener = gzip.open if a.vcf.endswith('.gz') else open
    with opener(a.vcf, 'rt') as fh:
        for line in fh:
            if line.startswith('##fileDate='):
                date = line.strip().split('=', 1)[1].replace('-', '')
            elif line.startswith('##SnpEffCmd='):
                snpeff_cmd = line.strip().split('=', 1)[1].strip('"')
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            info = dict(x.split('=', 1) for x in f[7].split(';') if '=' in x)
            cls = plp_class(info.get('CLNSIG', ''))
            if not cls:
                continue
            anns = [x.split('|') for x in info.get('ANN', '').split(',')]
            ann = next((x for x in anns if len(x) > 6 and x[6] not in plus), [])
            if len(ann) < 11 or 'missense_variant' not in ann[1]:
                continue
            m = HGVS_MISSENSE.match(ann[10])
            if not m:
                continue
            rows.append((ann[3], int(m.group(2)), m.group(1), m.group(3), ann[10], cls,
                         stars(info.get('CLNREVSTAT', '')), f[2], f[0], f[1], f[3], f[4],
                         1 if 'splice_region' in ann[1] else 0))
    if not date:
        sys.exit('ERROR: no ##fileDate in ' + a.vcf)

    os.makedirs(a.outdir, exist_ok=True)
    prefix = os.path.join(a.outdir, f'clinvar_{date}_protein_index')
    rows.sort(key=lambda r: (r[0], r[1], r[3], r[8], int(r[9])))
    with gzip.open(prefix + '.tsv.gz', 'wt') as out:
        out.write(f'## clinvar_version={date} mane={snpeff_cmd.split()[-2] if snpeff_cmd else "unknown"}\n')
        out.write('#GENE\tAA_POS\tAA_REF\tAA_ALT\tHGVS_P\tCLASS\tSTARS\tVCV\tCHROM\tPOS\tREF\tALT\tSPLICE_REGION\n')
        for r in rows:
            out.write('\t'.join(map(str, r)) + '\n')
    with open(prefix + '.manifest.txt', 'w') as m:
        m.write('resource:   ClinVar P/LP missense, MANE Select protein coordinates (ACMG PS1 / PM5)\n')
        m.write(f'version:    {date}\n')
        m.write(f'built:      {datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")}\n')
        m.write(f'input:      {os.path.basename(a.vcf)} (md5 {md5(a.vcf)})\n')
        m.write(f'snpeff:     {snpeff_cmd}\n')
        m.write(f'plus_clin:  {os.path.basename(a.mane_plus_clinical)} (md5 {md5(a.mane_plus_clinical)}, '
                f'{len(plus)} MANE Plus Clinical transcripts not used)\n')
        m.write(f'variants:   {len(rows)} (>=2 stars {sum(r[6] >= 2 for r in rows)}, '
                f'splice_region {sum(r[12] for r in rows)})\n')
        m.write(f'md5_tsv:    {md5(prefix + ".tsv.gz")}\n')
    sys.stdout.write(open(prefix + '.manifest.txt').read())


if __name__ == '__main__':
    main()
