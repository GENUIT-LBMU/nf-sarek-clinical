"""Somatic CNV harmonisation and ASCAT quality alerts (bin/cnv_shared.py)."""
import os
from conftest import ASSETS
import cnv_shared as c

SOM = os.path.join(ASSETS, 'gene_panels', 'somatic_38_oncogenes.txt')
GERM = os.path.join(ASSETS, 'gene_panels', 'germline_cancer_predisposition.txt')
HDR = 'sample\tevent_type\tgene\tsource\tchrom\tstart\tend\tcn\tlog2\tminor\tloh_type\n'


def report(tmp_path, rows):
    p = tmp_path / 'pair_somatic_cnv_report.tsv'
    p.write_text(HDR + ''.join('P\t' + '\t'.join(r) + '\n' for r in rows))
    return str(p)


def test_amp_needs_threshold_and_somatic_panel(tmp_path):
    p = report(tmp_path, [
        ('AMP', 'ERBB2', 'ASCAT', '17', '1', '2', '8', '1.9', '', ''),
        ('AMP', 'MYC', 'ASCAT', '8', '1', '2', '4', '1.0', '', ''),        # CN 4 < 5
    ])
    res = c.harmonize_cnv(p, '', SOM, GERM)
    assert [r['gene'] for r in res['amp']] == ['ERBB2']


def test_loh_only_in_germline_panel(tmp_path):
    p = report(tmp_path, [('2nd_hit_LOH', 'MSH2', 'ASCAT', '2', '1', '2', '2', '0', '0', 'cnLOH')])
    res = c.harmonize_cnv(p, '', SOM, GERM)
    assert [r['gene'] for r in res['loh']] == ['MSH2']


def test_missing_inputs_give_no_events():
    assert c.harmonize_cnv('', '', SOM, GERM) == {'homdel': [], 'amp': [], 'loh': []}


def test_ascat_alerts():
    assert c.ascat_alerts({'purity': '15', 'gof': '95'})         # purity < 20 %
    assert c.ascat_alerts({'purity': '60', 'gof': '85'})         # goodness of fit < 90 %
    assert c.ascat_alerts({'purity': '60', 'gof': '95'}) == []
    assert c.ascat_alerts({'purity': 'NA', 'gof': 'NA'}) == []


def test_ploidy_note_only_for_wgd_or_high_ploidy():
    assert 'whole-genome doubling' in c.ploidy_note({'ploidy': '3.1', 'wgd': '1'})
    assert c.ploidy_note({'ploidy': '2.0', 'wgd': '0'}) == ''
