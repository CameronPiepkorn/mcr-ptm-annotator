"""
mcr_ptm_annotator
~~~~~~~~~~~~~~~~~
Annotate McrA protein sequences with known post-translational modification
sites from the methyl-coenzyme M reductase active site.

Based on:
  Ermler et al. (1997) Science 278:1457  -- PDB 1MRO, 5 PTMs identified
  Nayak et al. (2017) eLife 6:e29218     -- thioamidation genetics
  Nayak et al. (2020) PLoS Biol 18:e3000507 -- PTM functional interactions

Quick start::

    from mcr_ptm_annotator import McrAPTMAnnotator, report

    annotator = McrAPTMAnnotator()
    hits = annotator.annotate_sequence(my_mcra_sequence, seq_id="MA0528")

    print(report.summary(hits, seq_id="MA0528"))
    report.to_tsv(hits, seq_id="MA0528", path="ptm_hits.tsv")
"""

from .annotator import McrAPTMAnnotator, PTMHit
from .ptm_database import KNOWN_PTMS, KnownPTM
from . import report

__version__ = "0.1.0"
__all__ = ["McrAPTMAnnotator", "PTMHit", "KNOWN_PTMS", "KnownPTM", "report"]
