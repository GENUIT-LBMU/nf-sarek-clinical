"""ACMG classification, MSI status and signature calls (bin/acmg_shared.py)."""
from conftest import getter
import acmg_shared as a

BASE = {
    'ANN[*].GENE': 'BRCA2', 'ANN[*].EFFECT': 'frameshift_variant', 'ANN[*].IMPACT': 'HIGH',
    'ANN[*].HGVS_C': 'c.100delA', 'ANN[*].HGVS_P': 'p.Lys34fs', 'ANN[*].RANK': '3/27',
    'ANN[*].AA_POS': '34', 'ANN[*].AA_LEN': '3418', 'ANN[*].FEATUREID': 'NM_000059.4',
    'gnomAD_AF': '0', 'gnomAD_FAF95': '0', 'gnomAD_AF_source': 'gnomAD_v4.1',
}
MISSENSE = dict(BASE, **{'ANN[*].EFFECT': 'missense_variant', 'ANN[*].IMPACT': 'MODERATE',
                         'ANN[*].HGVS_P': 'p.Ala100Val'})


def acmg(**over):
    return a.acmg_criteria_list(getter(dict(BASE, **over)))


def test_rare_frameshift_is_likely_pathogenic():
    cls, crit = acmg()
    assert cls == 'LP' and 'PVS1' in crit


def test_clinvar_pathogenic_two_stars_wins():
    assert acmg(ClinVar_sig='Pathogenic', ClinVar_stars='2')[0] == 'P'


def test_clinvar_benign_two_stars_wins():
    assert acmg(ClinVar_sig='Benign', ClinVar_stars='2')[0] == 'B'


def test_common_variant_is_benign_by_ba1():
    cls, crit = acmg(gnomAD_AF='0.12', gnomAD_FAF95='0.10')
    assert cls == 'B' and 'BA1' in crit


def test_clinvar_conflicting_is_not_pathogenic():
    d = dict(MISSENSE, ClinVar_sig='Conflicting_classifications_of_pathogenicity', ClinVar_stars='1')
    assert a.acmg_criteria_list(getter(d))[0] == 'VUS'


def test_rare_missense_without_evidence_is_vus():
    assert a.acmg_criteria_list(getter(MISSENSE))[0] == 'VUS'


def test_msi_status_thresholds():
    assert a.msi_status('3.5') == 'MSI-H'
    assert a.msi_status('3.49') == 'MSI borderline'
    assert a.msi_status('1.5') == 'MSI borderline'
    assert a.msi_status('1.4') == 'MSS'
    assert a.msi_status('NA') == 'NA' and a.msi_status('') == 'NA'


def test_variant_label_has_no_snpsift_placeholder():
    lbl = a.variant_label(getter({'ANN[*].GENE': 'MSH2', 'ANN[*].HGVS_C': 'c.1386+1G>A', 'ANN[*].HGVS_P': '.'}))
    assert lbl == 'c.1386+1G>A'


def test_mmr_signature_sets():
    assert {'ID1', 'ID2', 'ID7'} <= a.MMR_ID_SIGNATURES
    assert 'SBS6' in a.MMR_SBS_SIGNATURES and 'SBS3' in a.HRD_SBS_SIGNATURES


def _activities(tmp_path, rows):
    p = tmp_path / 'act.txt'
    p.write_text('Samples\tSBS1\tSBS5\tSBS6\n' + ''.join('\t'.join(map(str, r)) + '\n' for r in rows))
    return str(p)


def test_signature_needs_ten_percent(tmp_path):
    p = _activities(tmp_path, [('P1', 700, 205, 95), ('P2', 700, 200, 100)])
    assert 'SBS6' not in a.signature_calls(p, 'P1')['detected']      # 9.5 %
    assert 'SBS6' in a.signature_calls(p, 'P2')['detected']          # 10 %


def test_signature_not_assessable_below_100_sbs(tmp_path):
    p = _activities(tmp_path, [('P1', 40, 10, 30)])
    res = a.signature_calls(p, 'P1')
    assert res['assessable'] is False and res['detected'] == set()


def test_signature_row_matched_by_sample(tmp_path):
    p = _activities(tmp_path, [('OTHER', 0, 0, 900), ('P1', 900, 100, 0)])
    assert a.signature_calls(p, 'P1')['detected'] == {'SBS1', 'SBS5'}
