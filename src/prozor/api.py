"""The one public prozor module for other anndata_bridge packages."""

from __future__ import annotations

from prozor.inference.greedy import greedy_parsimony
from prozor.inference.ties import TieCandidate
from prozor.matching.annotation import (
    ProteinSequenceRecord,
    annotate_peptides,
    annotate_peptides_streaming,
)
from prozor.matching.automaton import resolve_backend

__all__ = [
    "ProteinSequenceRecord",
    "TieCandidate",
    "annotate_peptides",
    "annotate_peptides_streaming",
    "greedy_parsimony",
    "resolve_backend",
]
