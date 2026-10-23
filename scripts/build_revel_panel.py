#!/usr/bin/env python3
"""
build_revel_panel.py — build the local REVEL resource used by the germline
annotation step (REVEL_ANNOTATE).

Output (in --outdir):
  revel_v1.3_panel.tsv.gz(.tbi)   one row per missense allele inside the panel
                                  regions: CHROM POS REF ALT REVEL_SCORE REVEL_TX REVEL_AA
  revel_v1.3_panel.manifest.txt   version, source, regions, MANE release, md5

REVEL v1.3 gives one score per amino-acid change, with the Ensembl transcripts
that share it. For each allele the score of the MANE Select transcript is kept
(REVEL_TX = MANE); when no row lists a MANE Select transcript (REVEL v1.3 uses
older Ensembl transcript IDs, e.g. ATM, CASR, FANCD2), the highest score is kept
(REVEL_TX = other_tx) and REVEL_AA lists the score of every amino-acid change
("R>Q:0.123,R>W:0.456"), so the pipeline can pick the one matching the MANE HGVS_p.

Input:
  --revel-dir  directory with the REVEL v1.3 zips (revel-v1.3_segments_chrom_*.zip)
               https://zenodo.org/records/7072866 (ODbL)
  --regions    panel regions BED (bgzipped), e.g. gnomad_v4.1_panel_regions.bed.gz
               written by scripts/build_gnomad_panel.sh
  --mane       MANE summary (MANE.GRCh38.v1.x.summary.txt.gz, NCBI)
  --outdir     output directory

Needs python3, bgzip and tabix:
  python3 scripts/build_revel_panel.py --revel-dir /path/to/revel \\
      --regions /path/to/refs/gnomad/gnomad_v4.1_panel_regions.bed.gz \\
      --mane /path/to/MANE.GRCh38.v1.5.summary.txt.gz --outdir /path/to/refs/revel

REVEL: Ioannidis et al., Am J Hum Genet 2016.
"""
import argparse
import bisect
import collections
import csv
import datetime
import glob
import gzip
import hashlib
import io
import os
import subprocess
import sys
import zipfile

VERSION = '1.3'


def md5(path):
    h = hashlib.md5()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_regions(path):
    reg = collections.defaultdict(list)
    with gzip.open(path, 'rt') as fh:
        for line in fh:
            f = line.split('\t')
            reg[f[0]].append((int(f[1]), int(f[2])))
    for c in reg:
        reg[c].sort()
    return reg


def in_regions(reg, chrom, pos):
    iv = reg.get(chrom, [])
    i = bisect.bisect_right(iv, (pos, float('inf'))) - 1
    return i >= 0 and iv[i][0] < pos <= iv[i][1]


def load_mane(path):
    tx = set()
    with gzip.open(path, 'rt') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            if row.get('MANE_status') == 'MANE Select':
                tx.add(row['Ensembl_nuc'].split('.')[0])
    return tx


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--revel-dir', required=True)
    ap.add_argument('--regions', required=True)
    ap.add_argument('--mane', required=True)
    ap.add_argument('--outdir', required=True)
    a = ap.parse_args()

    zips = sorted(glob.glob(os.path.join(a.revel_dir, 'revel-v1.3_segments_chrom_*.zip')))
    if len(zips) < 24:
        sys.exit(f'ERROR: expected 24+ REVEL v1.3 zips in {a.revel_dir}, found {len(zips)}')
    reg = load_regions(a.regions)
    mane = load_mane(a.mane)
    os.makedirs(a.outdir, exist_ok=True)
    prefix = os.path.join(a.outdir, f'revel_v{VERSION}_panel')

    best = {}        # (chrom, pos, ref, alt) -> (is_mane, score)
    by_aa = collections.defaultdict(dict)   # key -> {"R>Q": score}
    n_rows = n_nopos = 0
    for z in zips:
        with zipfile.ZipFile(z) as zf:
            for name in zf.namelist():
                if not name.endswith('.csv'):
                    continue
                with zf.open(name) as raw:
                    for row in csv.DictReader(io.TextIOWrapper(raw)):
                        n_rows += 1
                        pos = row['grch38_pos']
                        if pos in ('.', ''):
                            n_nopos += 1
                            continue
                        chrom = 'chr' + row['chr']
                        pos = int(pos)
                        if not in_regions(reg, chrom, pos):
                            continue
                        try:
                            score = float(row['REVEL'])
                        except ValueError:
                            continue
                        is_mane = any(t.split('.')[0] in mane
                                      for t in row['Ensembl_transcriptid'].split(';'))
                        key = (chrom, pos, row['ref'], row['alt'])
                        aa = row['aaref'] + '>' + row['aaalt']
                        by_aa[key][aa] = max(score, by_aa[key].get(aa, score))
                        cur = best.get(key)
                        # MANE row wins; otherwise keep the highest score
                        if (cur is None or (is_mane and not cur[0])
                                or (is_mane == cur[0] and score > cur[1])):
                            best[key] = (is_mane, score)
        sys.stderr.write(f'[revel] {os.path.basename(z)}: {len(best)} alleles so far\n')

    def chrom_key(c):
        c = c[3:]
        return (0, int(c)) if c.isdigit() else (1, c)

    tmp = prefix + '.tsv'
    counts = collections.Counter()
    with open(tmp, 'w') as out:
        out.write(f'## revel_version={VERSION} source=zenodo.org/records/7072866 '
                  f'mane={os.path.basename(a.mane)}\n')
        out.write('#CHROM\tPOS\tREF\tALT\tREVEL_SCORE\tREVEL_TX\tREVEL_AA\n')
        for key in sorted(best, key=lambda k: (chrom_key(k[0]), k[1], k[2], k[3])):
            is_mane, score = best[key]
            tx = 'MANE' if is_mane else 'other_tx'
            counts[tx] += 1
            aas = ','.join(f'{k}:{v:.3f}' for k, v in sorted(by_aa[key].items()))
            out.write(f'{key[0]}\t{key[1]}\t{key[2]}\t{key[3]}\t{score:.3f}\t{tx}\t{aas}\n')
    subprocess.run(['bgzip', '-f', tmp], check=True)
    subprocess.run(['tabix', '-f', '-s1', '-b2', '-e2', '-c#', prefix + '.tsv.gz'], check=True)

    with open(prefix + '.manifest.txt', 'w') as m:
        m.write('resource:   REVEL v1.3 missense scores, panel subset (MANE Select score per allele)\n')
        m.write(f'version:    {VERSION}\n')
        m.write('source:     https://zenodo.org/records/7072866 (ODbL)\n')
        m.write(f'built:      {datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")}\n')
        m.write(f'regions:    {os.path.basename(a.regions)} (md5 {md5(a.regions)})\n')
        m.write(f'mane:       {os.path.basename(a.mane)} ({len(mane)} MANE Select transcripts)\n')
        m.write(f'rows_read:  {n_rows} ({n_nopos} without GRCh38 position)\n')
        m.write(f'alleles:    {len(best)} (MANE {counts["MANE"]}, other_tx {counts["other_tx"]})\n')
        m.write(f'md5_tsv:    {md5(prefix + ".tsv.gz")}\n')
    sys.stdout.write(open(prefix + '.manifest.txt').read())


if __name__ == '__main__':
    main()
