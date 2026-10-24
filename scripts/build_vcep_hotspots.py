#!/usr/bin/env python3
"""
build_vcep_hotspots.py — print the cancerhotspots.org counts embedded in
bin/acmg_shared.py (VCEP_HOTSPOTS) for the ClinGen VCEP PM1 rules that use them
(TP53 GN009, VHL GN078).

Input: hotspots_v2.xls from https://www.cancerhotspots.org/files/hotspots_v2.xls
(Chang et al., Cancer Discov 2018). Missense changes only (sheet SNV-hotspots;
nonsense, synonymous and splice entries are skipped). Needs pandas + xlrd.

  python3 scripts/build_vcep_hotspots.py hotspots_v2.xls TP53 VHL
"""
import sys

import pandas as pd

AA3 = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu',
       'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe',
       'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'}


def main():
    xls, genes = sys.argv[1], sys.argv[2:]
    s = pd.read_excel(xls, sheet_name='SNV-hotspots')
    print('VCEP_HOTSPOTS = {   # cancerhotspots.org v2: "p.Arg175His" -> somatic count')
    for g in genes:
        out = {}
        for _, r in s[s.Hugo_Symbol == g].iterrows():
            pos = str(r.Amino_Acid_Position)
            ref = str(r.Reference_Amino_Acid).split(':')[0]
            alt, n = str(r.Variant_Amino_Acid).split(':')
            if not pos.isdigit() or ref not in AA3 or alt not in AA3 or alt == ref:
                continue
            out['p.{}{}{}'.format(AA3[ref], pos, AA3[alt])] = int(n)
        items = sorted(out.items(), key=lambda kv: (int(''.join(c for c in kv[0] if c.isdigit())), kv[0]))
        print("    '{}': {{".format(g))
        line = '       '
        for k, v in items:
            tok = " '{}': {},".format(k, v)
            if len(line) + len(tok) > 92:
                print(line)
                line = '       '
            line += tok
        print(line)
        print('    },')
    print('}')


if __name__ == '__main__':
    main()
