"""Shared report helpers (bin/report_shared.py)."""
import report_shared as r


def test_fmt_hrd_prob():
    assert r.fmt_hrd_prob('1.33e-06') == '<0.01'
    assert r.fmt_hrd_prob(0.4567) == '0.46'
    assert r.fmt_hrd_prob('NA') == 'NA' and r.fmt_hrd_prob('') == 'NA'


def test_fmt_num_keeps_integers_and_text():
    assert r.fmt_num('0.029822426928536254') == '0.02982'
    assert r.fmt_num('12') == '12' and r.fmt_num('HET') == 'HET'


def test_unavailable_text(tmp_path, monkeypatch):
    monkeypatch.setenv('REPORT_NOT_RUN', 'germline_cnv,coverage')
    assert 'not analysed in this run' in r.unavailable_text('germline_cnv', '')
    assert 'no result for this sample' in r.unavailable_text('pharmaco', '')
    f = tmp_path / 'x.tsv'
    f.write_text('a\n')
    assert r.unavailable_text('pharmaco', str(f)) == ''
    assert r.kpi_value('coverage', '30') == 'Not run' and r.kpi_value('tmb', '5') == '5'


def test_trial_queries_use_tumor_type_and_biomarkers():
    corr = {'results': [{'gene': 'MLH1', 'acmg': 'LP'}, {'gene': 'ATM', 'acmg': 'VUS'}],
            'biallelic_somatic': [{'gene': 'MSH2'}]}
    q = r.trial_queries(corr, ['KRAS'], 67.9, 'MSI-H', 'HRD-', 'Colorectal Cancer|Colon Cancer')
    terms = [x[2] for x in q]
    assert terms[:3] == ['KRAS', 'MLH1', 'MSH2']            # VUS not queried
    assert all(x[1] in ('Colorectal Cancer', 'solid tumor') for x in q)
    assert any('tumor mutational burden' in t for t in terms)
    assert not any('homologous recombination' in t for t in terms)


def test_trials_for_the_opposite_msi_status_are_dropped():
    mss_trial = {'title': 'Pembrolizumab in MSS colorectal cancer', 'conds': 'Colorectal Cancer'}
    assert r._wrong_msi(mss_trial, 'MSI-H')
    assert not r._wrong_msi(mss_trial, 'MSS')


def test_versions_from_sarek_yaml(tmp_path, monkeypatch):
    y = tmp_path / 'versions.yml'
    y.write_text('ASCAT:\n  ascat: 3.1.1\nMSISENSORPRO_MSISOMATIC:\n  msi-sensor-pro: 1.2.0\n'
                 'STRELKA_SOMATIC:\n  strelka: 2.9.10\n')
    monkeypatch.setenv('REPORT_SAREK_VERSIONS', str(y))
    assert r.tool_label('ASCAT', 'ascat') == 'ASCAT 3.1.1'
    assert r.tool_label('CNVkit', 'cnvkit') == 'CNVkit'
