#!/usr/bin/env python3
"""
build_mane_exons.py — exon table of the MANE Select transcripts from the NCBI
MANE GTF, in the exons_mane.txt format read by scripts/build_gnomad_panel.sh
(chr, start, end 1-based inclusive, gene_symbol, exon, transcript_id).

Use the MANE release of the snpEff database (GRCh38.mane.1.2.refseq):
  curl -O https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.2/MANE.GRCh38.v1.2.refseq_genomic.gtf.gz
  python3 scripts/build_mane_exons.py --gtf MANE.GRCh38.v1.2.refseq_genomic.gtf.gz \\
      --out /data/refs/mane/exons_mane_v1.2.txt
"""
import argparse
import gzip
import re


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gtf', required=True)
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    rows = []
    with gzip.open(a.gtf, 'rt') as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            f = line.rstrip('\n').split('\t')
            if len(f) < 9 or f[2] != 'exon' or 'tag "MANE Select"' not in f[8]:
                continue
            attr = dict(re.findall(r'(\w+) "([^"]*)"', f[8]))
            rows.append((f[0], int(f[3]), int(f[4]), attr['gene'], attr.get('exon_number', '.'),
                         attr['transcript_id']))
    rows.sort(key=lambda r: (r[0], r[1]))
    with open(a.out, 'w') as out:
        out.write('chr\tstart\tend\tgene_symbol\texon\ttranscript_id\n')
        for r in rows:
            out.write('\t'.join(map(str, r)) + '\n')
    print(f'{len(rows)} exons, {len({r[5] for r in rows})} MANE Select transcripts -> {a.out}')


if __name__ == '__main__':
    main()
