"""Report model (bin/report_model.py): the single source of both reports."""
from types import SimpleNamespace
import report_model as rm


def write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


def test_metrics_and_statuses(tmp_path):
    a = SimpleNamespace(
        tmb=write(tmp_path, 'tmb.tsv', 'sample\tTMB_mut_per_mb\nP\t12.5\n'),
        msi=write(tmp_path, 'msi', 'Total_Number_of_Sites\tNumber_of_Somatic_Sites\t%\n1000\t40\t4.00\n'),
        ascat=write(tmp_path, 'met.txt', 'purity\tploidy\tgoodness_of_fit\tLOH\tGI\tWGD\n0.52\t2.5\t87\t0.3\t0.9\t1\n'),
        hrd=write(tmp_path, 'hrd.tsv', 'samples\thrd.prob\tstatus\nP\t0.6\tok\n'),
        cov_n=write(tmp_path, 'n.summary.txt', 'chrom\tlength\tbases\tmean\ntotal_region\t1\t1\t110.5\n'),
        cov_t=write(tmp_path, 't.summary.txt', 'chrom\tlength\tbases\tmean\ntotal_region\t1\t1\t46.1\n'))
    m = rm.build_metrics(a)
    assert m['tmb'] == '12.5' and m['tmb_status'] == 'TMB-H'
    assert m['msi'] == '4.00' and m['msi_status'] == 'MSI-H'
    assert m['purity'] == '52.0'
    assert m['hrd_status'] == 'HRD+' and m['hrd_prob_display'] == '0.60'
    assert m['cov_n'] == '110.5' and m['cov_t'] == '46.1'
    assert m['low_tumor_coverage'] is True          # 46X < 50X


def test_hrd_without_probability_keeps_the_reason(tmp_path):
    a = SimpleNamespace(tmb='', msi='', ascat='', cov_n='', cov_t='',
                        hrd=write(tmp_path, 'hrd.tsv', 'samples\thrd.prob\tstatus\nP\tNA\ttoo few SNVs\n'))
    assert rm.build_metrics(a)['hrd_status'] == 'NA (too few SNVs)'


def test_germline_cnv_segments_panel_qa_and_lowmap(tmp_path):
    hdr = 'CHROM\tPOS\tEND\tALT\tCN\tNP\tQA\tQS\tLOWMAP_SEG\tgen\n'
    rows = [
        'chr2\t100\t200\t<DEL>\t1\t3\t50\t60\t0\tMSH2',     # reported
        'chr2\t100\t200\t<DEL>\t1\t3\t50\t60\t0\tEPCAM',    # same segment, second gene
        'chr22\t5\t9\t<DEL>\t0\t2\t40\t200\t1.0\tCHEK2',    # low mappability: apart
        'chr17\t1\t2\t<DUP>\t3\t2\t2\t5\t0\tBRCA1',         # QA below 4
        'chr7\t1\t2\t<DUP>\t3\t9\t90\t90\t0\tEGFR',         # not in the germline panel
        'chr3\t1\t2\t.\t2\t9\t90\t90\t0\tMLH1',             # not an event
    ]
    p = write(tmp_path, 'seg.tsv', hdr + '\n'.join(rows) + '\n')
    res = rm.germline_cnv(p, {'MSH2', 'EPCAM', 'CHEK2', 'BRCA1', 'MLH1'})
    assert [s['genes'] for s in res['segments']] == [['MSH2', 'EPCAM']]
    assert [s['genes'] for s in res['lowmap']] == [['CHEK2']]


def test_signatures_aetiology_and_pair_row(tmp_path):
    p = write(tmp_path, 'act.txt', 'Samples\tSBS1\tSBS6\tSBS999\nOTHER\t5\t5\t5\nP\t300\t100\t0\n')
    s = rm.signatures(p, 'P', 'SBS')
    assert [x['signature'] for x in s['assigned']] == ['SBS1', 'SBS6']
    assert s['assigned'][1]['category'] == 'mmr' and s['assigned'][0]['percent'] == 75.0
    assert s['calls']['detected'] == ['SBS1', 'SBS6']


PGX_HDR = 'CHROM\tPOS\tID\tREF\tALT\tGEN[0].GT\n'


def pgx(tmp_path, *rows):
    return rm.pharmaco(write(tmp_path, 'pgx.tsv', PGX_HDR + ''.join(r + '\n' for r in rows)))


def gene_of(rows, gene):
    return next(r for r in rows if r['gene'] == gene)


def test_dpyd_heterozygous_no_function_is_intermediate_50pct(tmp_path):
    r = gene_of(pgx(tmp_path, 'chr1\t1\trs3918290\tC\tT\t0/1'), 'DPYD')
    assert r['activity'] == '1' and r['phenotype'] == 'Intermediate metabolizer'
    assert '50%' in r['recommendation'] and 'CPIC' in r['guideline']


def test_dpyd_two_variants_unknown_phase_taken_in_trans(tmp_path):
    r = gene_of(pgx(tmp_path, 'chr1\t1\trs3918290\tC\tT\t0/1',
                    'chr1\t2\trs67376798\tT\tA\t0/1'), 'DPYD')
    assert r['activity'] == '0.5' and r['phenotype'] == 'Poor metabolizer' and r['alert'] == 'HIGH'


def test_dpyd_homozygous_no_function_avoid(tmp_path):
    r = gene_of(pgx(tmp_path, 'chr1\t1\trs3918290\tC\tT\t1/1'), 'DPYD')
    assert r['activity'] == '0' and r['recommendation'].startswith('Avoid')


def test_dpyd_decreased_function_heterozygous(tmp_path):
    r = gene_of(pgx(tmp_path, 'chr1\t1\trs56038477\tC\tT\t0|1'), 'DPYD')
    assert r['activity'] == '1.5' and r['phenotype'] == 'Intermediate metabolizer'


def test_ugt1a1_alleles_and_repeat_length(tmp_path):
    two = gene_of(pgx(tmp_path, 'chr2\t1\trs4148323\tG\tA\t0/1',
                      'chr2\t2\trs3064744\tCAT\tCATAT\t0/1'), 'UGT1A1')
    assert two['phenotype'] == 'Poor metabolizer' and '70%' in two['recommendation']
    assert 'DPWG' in two['guideline'] and two['allele'].startswith('*6')
    one = gene_of(pgx(tmp_path, 'chr2\t2\trs3064744\tCAT\tCATAT\t0/1'), 'UGT1A1')
    assert one['phenotype'] == 'Intermediate metabolizer'
    assert pgx(tmp_path, 'chr2\t2\trs3064744\tCATAT\tCAT\t0/1') == []     # TA deletion (*36)
    assert pgx(tmp_path, 'chr2\t1\trs4148323\tG\tA\t0/0') == []           # reference call



SHEET = ('patient,sex,status,sample,lane,fastq_1,fastq_2,patient_name,birth_date,diagnosis,'
         'accession_id,specimen_type,tumor_content\n'
         'P1,XX,0,N1,L1,a,b,Ana Test,1970-01-01,Colon adenocarcinoma,A-2,Blood,\n'
         'P1,XX,1,T1,L1,a,b,Ana Test,1970-01-01,Colon adenocarcinoma,A-1,FFPE,60\n'
         'P1,XX,1,T1,L2,a,b,Ana Test,1970-01-01,Colon adenocarcinoma,A-1,FFPE,60\n')


def test_identification_from_samplesheet(tmp_path):
    idn = rm.identification(write(tmp_path, 's.csv', SHEET), 'T1', 'N1')
    assert idn['patient']['patient_name'] == 'Ana Test' and idn['patient']['sex'] == 'XX'
    assert idn['patient']['patient'] == 'P1'
    assert idn['tumor']['tumor_content'] == '60' and idn['tumor']['specimen_type'] == 'FFPE'
    assert idn['normal']['accession_id'] == 'A-2' and idn['normal']['tumor_content'] == ''
    assert idn['warnings'] == []


def test_identification_conflicts_and_absence(tmp_path):
    bad = SHEET.replace('P1,XX,1,T1,L2,a,b,Ana Test', 'P1,XX,1,T1,L2,a,b,Ana Other')
    idn = rm.identification(write(tmp_path, 's.csv', bad), 'T1', 'N1')
    assert any('patient_name' in w for w in idn['warnings'])
    none = rm.identification('', 'T1', 'N1')
    assert none['source'] == '' and none['patient']['patient_name'] == ''


def test_plain_sarek_samplesheet_has_no_identification(tmp_path):
    plain = 'patient,sex,status,sample,lane,fastq_1,fastq_2\nP1,XX,1,T1,L1,a,b\n'
    assert rm.identification(write(tmp_path, 's.csv', plain), 'T1', 'N1')['source'] == ''


def test_pgx_variant_outside_capture_is_not_assessable(tmp_path):
    bed = write(tmp_path, 'capture.bed', 'chr1\t97000000\t98000000\nchr2\t233760285\t233761153\n')
    assert rm.pharmaco_not_assessable(bed) == ['UGT1A1 *28 (rs3064744)']
    assert rm.pharmaco_not_assessable('') == []          # no capture BED: not checked
