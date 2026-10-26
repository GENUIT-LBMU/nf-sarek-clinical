"""The report and correlation scripts load (a broken import or a name moved to a
shared module fails here, not at the end of a run)."""
import importlib

import pytest


@pytest.mark.parametrize('module', [
    'clinical_pdf_report', 'germline_somatic_correlate', 'final_report_sections',
    'acmg_shared', 'amp_tier', 'cnv_shared', 'report_shared', 'gene_sets',
    'report_model', 'render_html',
])
def test_module_imports(module):
    importlib.import_module(module)
