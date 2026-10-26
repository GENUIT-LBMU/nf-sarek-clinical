"""
report_shared.py — shared by the HTML report (FINAL_REPORT) and the PDF report.

Analyses switched off in a run arrive in REPORT_NOT_RUN (comma-separated keys,
workflows/sarek_clinical.nf). A section whose input is missing says why:
"not analysed in this run" (switched off) or "no result for this sample"
(the analysis ran but produced nothing for it), never an empty negative
result such as "no variants detected".
"""

import os

# key -> (section title, parameter that switches it off)
ANALYSES = {
    'germline_snv':     ('Germline SNVs',                      'skip_germline_snv'),
    'germline_cnv':     ('Germline CNVs (GATK gCNV)',          'skip_germline_cnv'),
    'pharmaco':         ('Pharmacogenomics',                   'skip_pharmaco'),
    'somatic_snv':      ('Somatic SNVs',                       'skip_somatic_snv'),
    'somatic_cnv':      ('Somatic CNVs (ASCAT / CNVkit)',      'skip_somatic_cnv'),
    'somatic_cnv_gatk': ('Somatic CNVs (GATK)',                'skip_somatic_cnv_gatk'),
    'tmb':              ('TMB',                                'skip_tmb'),
    'signatures':       ('Mutational signatures',              'skip_signatures'),
    'hrd':              ('HRD',                                'skip_hrd'),
    'coverage':         ('Per-gene coverage',                  'skip_coverage'),
}


def not_run():
    """Set of analysis keys switched off in this run."""
    return {k.strip() for k in os.environ.get('REPORT_NOT_RUN', '').split(',') if k.strip()}


def is_not_run(key):
    return key in not_run()


def _present(path):
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def unavailable_text(key, *paths):
    """Why a section has no data, or '' when at least one of its inputs exists."""
    title, flag = ANALYSES[key]
    if key in not_run():
        return f'{title}: not analysed in this run (--{flag} true).'
    if any(_present(p) for p in paths):
        return ''
    return f'{title}: no result for this sample (the analysis produced no output file).'


def kpi_value(key, value):
    """Dashboard value: 'Not run' for an analysis switched off, else the value."""
    return 'Not run' if key in not_run() else value


# ── Versions (R3) ────────────────────────────────────────────────────────────
# REPORT_INFO: JSON written by the report task from workflows/sarek_clinical.nf
# (pipeline version, git commit, sarek release, resource files).
# REPORT_SAREK_VERSIONS: sarek's pipeline_info/nf_core_sarek_software_mqc_versions.yml.

import json as _json


def report_info():
    path = os.environ.get('REPORT_INFO', '')
    if path and os.path.isfile(path):
        with open(path) as fh:
            return _json.load(fh)
    return {}


def sarek_tool_versions():
    """{tool: version} from sarek's software versions YAML (tool lines only)."""
    path = os.environ.get('REPORT_SAREK_VERSIONS', '')
    out = {}
    if not path or not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            if line.startswith((' ', '\t')) and ':' in line:
                tool, ver = line.strip().split(':', 1)
                tool, ver = tool.strip(), ver.strip().strip('"\'')
                if tool and ver and tool not in out:
                    out[tool] = ver
    return out


def tool_label(name, key):
    """'Strelka2 2.9.10' when sarek recorded the version, else the bare name."""
    ver = sarek_tool_versions().get(key, '')
    return f'{name} {ver}' if ver else name


def pipeline_label():
    """'nf-sarek-clinical 1.0.0 · commit 9c180be' (with '+ local changes' when dirty)."""
    info = report_info()
    if not info:
        return 'nf-sarek-clinical'
    return f"{info.get('pipeline', 'nf-sarek-clinical')} · commit {info.get('commit', 'unknown')}"


def sarek_label():
    info = report_info()
    return f"nf-core/sarek {info.get('sarek') or 'release not recorded'}"


# sarek tools shown in the versions table: (display name, key in the YAML)
SAREK_TOOLS = [
    ('BWA-MEM', 'bwa'), ('fastp', 'fastp'), ('GATK', 'gatk4'), ('Strelka2', 'strelka'),
    ('ASCAT', 'ascat'), ('CNVkit', 'cnvkit'), ('MSIsensor-pro', 'msisensor-pro'),
    ('mosdepth', 'mosdepth'), ('snpEff (sarek)', 'snpeff'),
]


def version_rows():
    """[(label, value)] for the versions table of both reports."""
    info = report_info()
    tools = sarek_tool_versions()
    rows = [('Pipeline', pipeline_label()), ('Upstream pipeline', sarek_label()),
            ('Nextflow', info.get('nextflow', 'NA'))]
    rows += [(name, tools[key]) for name, key in SAREK_TOOLS if key in tools]
    rows += [(k, v) for k, v in (info.get('resources') or {}).items() if v]
    return rows


# ── Number formats (R6): the same value prints the same way everywhere ───────

def fmt_hrd_prob(val):
    """HRD probability: '<0.01' below 0.01, else 2 decimals; 'NA' when missing."""
    try:
        f = float(str(val).strip())
    except (ValueError, TypeError):
        return 'NA'
    if f != f:
        return 'NA'
    return '<0.01' if f < 0.01 else f'{f:.2f}'


def fmt_num(val, sig=4):
    """Numeric table cell with `sig` significant digits; text unchanged."""
    s = str(val).strip()
    try:
        f = float(s)
    except (ValueError, TypeError):
        return s
    if f != f or s.isdigit() or (s.startswith('-') and s[1:].isdigit()):
        return s
    return f'{f:.{sig}g}'


# ── Clinical trials (R7): one query plan for both reports ───────────────────
# Genes with an actionable somatic variant (AMP Tier I/II), the genes of the
# germline-somatic findings and the biomarkers, each searched with the patient's
# tumor type as condition; MSI-H and TMB-H also in any solid tumor. Trials
# restricted to the opposite MSI status are dropped. Live query: the report
# prints its date, and says so when ClinicalTrials.gov cannot be reached.

ACTIONABLE_TIERS = ('IA', 'IB', 'IIC', 'IID')
_MSI_TERM = 'microsatellite instability OR mismatch repair deficient OR dMMR'
_MSS_WORDS = ('mss', 'microsatellite stable', 'pmmr', 'mismatch repair proficient')
_MSI_WORDS = ('msi-h', 'msi high', 'dmmr', 'mismatch repair deficient', 'microsatellite instability-high')


def _condition(tumor_label):
    terms = [t.strip() for t in str(tumor_label or '').split('|') if t.strip()]
    return terms[0] if terms else 'cancer'


def trial_queries(corr_data, actionable_genes, tmb, msi_st, hrd_status, tumor_label, max_genes=4):
    """[(label, condition, term)] in order of priority."""
    cond = _condition(tumor_label)
    genes = list(actionable_genes or [])
    for r in (corr_data or {}).get('results', []):
        if r.get('acmg') in ('P', 'LP'):
            genes.append(r.get('gene', ''))
    for r in (corr_data or {}).get('biallelic_somatic', []):
        genes.append(r.get('gene', ''))
    genes = [g for g in dict.fromkeys(genes) if g][:max_genes]
    q = [(f'{g} · {cond}', cond, g) for g in genes]
    if msi_st == 'MSI-H':
        q += [(f'MSI-H / dMMR · {cond}', cond, _MSI_TERM),
              ('MSI-H / dMMR · any solid tumor', 'solid tumor', _MSI_TERM)]
    try:
        if float(tmb) >= 10:
            q.append(('TMB-H · any solid tumor', 'solid tumor', 'tumor mutational burden'))
    except (ValueError, TypeError):
        pass
    if str(hrd_status) == 'HRD+':
        q.append((f'HRD · {cond}', cond, 'homologous recombination deficiency OR PARP inhibitor'))
    return q


def _wrong_msi(trial, msi_st):
    text = (trial['title'] + ' ' + trial['conds']).lower()
    if msi_st == 'MSI-H':
        return any(w in text for w in _MSS_WORDS) and not any(w in text for w in _MSI_WORDS)
    return any(w in text for w in _MSI_WORDS) and not any(w in text for w in _MSS_WORDS)


def clinical_trials(queries, msi_st='', per_query=3, timeout=8):
    """{'date', 'trials': [...], 'error': str|None} from ClinicalTrials.gov API v2."""
    import datetime
    import urllib.parse
    import urllib.request
    out = {'date': datetime.date.today().isoformat(), 'trials': [], 'error': None}
    seen = set()
    for label, cond, term in queries:
        params = urllib.parse.urlencode({
            'query.cond': cond, 'query.term': term,
            'filter.overallStatus': 'RECRUITING',
            'fields': 'NCTId,BriefTitle,Phase,OverallStatus,Condition,StartDate',
            'pageSize': str(per_query * 3), 'format': 'json',
        })
        try:
            req = urllib.request.Request('https://clinicaltrials.gov/api/v2/studies?' + params,
                                         headers={'User-Agent': 'nf-sarek-clinical'})
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                data = _json.loads(resp.read().decode('utf-8'))
        except Exception as exc:   # network / API down: the report says so
            out['error'] = str(exc)
            break
        kept = 0
        for s in data.get('studies', []):
            proto = s.get('protocolSection', {})
            ident = proto.get('identificationModule', {})
            stat = proto.get('statusModule', {})
            t = {
                'nct':    ident.get('nctId', ''),
                'title':  ident.get('briefTitle', ''),
                'phase':  ', '.join(proto.get('designModule', {}).get('phases', [])) or 'N/A',
                'status': stat.get('overallStatus', ''),
                'conds':  ', '.join(proto.get('conditionsModule', {}).get('conditions', [])[:2]),
                'start':  (stat.get('startDateStruct') or {}).get('date', ''),
                'query':  label,
            }
            if not t['nct'] or t['nct'] in seen or _wrong_msi(t, msi_st):
                continue
            seen.add(t['nct'])
            out['trials'].append(t)
            kept += 1
            if kept >= per_query:
                break
    return out
