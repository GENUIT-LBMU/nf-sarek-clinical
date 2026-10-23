"""
ACMG germline classification helpers for the nf-sarek-clinical final report.
Imported by modules/local/final_report/main.nf at run time via sys.path.

The actual classification logic lives in `bin/acmg_shared.py` (single source of
truth for ACMG, in silico panel, and gene/signature sets). This module exposes
the (row, hdr) signature expected by callers in the final-report scripts.
"""
from acmg_shared import (
    acmg_criteria_list as _acmg_criteria_list,
    count_insil as _count_insil_shared,
    getter_for_row_hdr,
    INSIL_DAMAGING as _INSIL_DAMAGING,  # re-exported for callers that introspect
    INSIL_BENIGN as _INSIL_BENIGN,
)

_BENIGN_CRIT = {'BA1', 'BS1', 'BP4_Sup', 'BP4_Mod', 'BP4_St'}


def _count_insil(col_fn):
    """Count damaging vs benign votes given a column-getter callable."""
    return _count_insil_shared(col_fn)


def acmg_criteria_list(row, hdr):
    """Returns (cls, criteria_list) for a (row:list, hdr:dict[str,int]) pair."""
    return _acmg_criteria_list(getter_for_row_hdr(row, hdr))


def acmg_class(row, hdr):
    """Thin wrapper — returns only the classification string."""
    cls, _ = acmg_criteria_list(row, hdr)
    return cls
