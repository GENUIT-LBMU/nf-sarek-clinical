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
