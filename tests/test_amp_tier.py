"""AMP/ASCO/CAP tiers from CIViC evidence (bin/amp_tier.py)."""
import amp_tier as t

TERMS = {'colorectal cancer'}


def rec(**k):
    r = {'gene': 'KRAS', 'profile': 'KRAS G13D', 'entity_type': 'evidence', 'status': 'accepted',
         'direction': 'SUPPORTS', 'significance': 'RESISTANCE', 'disease': 'Colorectal Cancer',
         'therapies': 'Cetuximab', 'level': 'B', 'amp_category': ''}
    r.update(k)
    return r


def tier(*records):
    return t.classify(list(records), TERMS)[0]


def test_level_a_b_same_tumor():
    assert tier(rec(level='A')) == 'IA'
    assert tier(rec(level='B')) == 'IB'


def test_level_a_b_other_tumor_is_iic():
    assert tier(rec(level='A', disease='Lung Cancer')) == 'IIC'
    assert tier(rec(level='B', disease='Lung Cancer')) == 'IIC'


def test_level_c_d_e():
    assert tier(rec(level='C')) == 'IIC'
    assert tier(rec(level='D')) == 'IID'
    assert tier(rec(level='E')) == 'IID'
    assert tier(rec(level='C', disease='Lung Cancer')) == 'III'


def test_submitted_and_not_supporting_never_set_the_tier():
    assert tier(rec(level='A', status='submitted')) == 'III'
    assert tier(rec(level='A', direction='DOES_NOT_SUPPORT')) == 'III'


def test_best_record_wins():
    assert tier(rec(level='D'), rec(level='B')) == 'IB'


def test_pan_cancer_counts_as_same_tumor():
    assert tier(rec(level='B', disease='Solid Tumor')) == 'IB'


def test_population_polymorphism_is_tier_iv():
    get = lambda k: {'POP_COMMON': '1'}.get(k, '')
    assert t.variant_tier([rec(level='A')], TERMS, get)[0] == 'IV'


def test_no_tumor_type_caps_at_iic():
    assert t.classify([rec(level='A')], set())[0] == 'IIC'


def test_tert_note_when_promoter_outside_capture(tmp_path):
    bed = tmp_path / 'capture.bed'
    bed.write_text('chr5\t1253000\t1295000\n')
    assert 'not assessable' in t.tert_promoter_note(str(bed))
    bed.write_text('chr5\t1295000\t1295200\n')
    assert t.tert_promoter_note(str(bed)) == ''
