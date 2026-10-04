"""The one public prozor module for other anndata_bridge packages."""

from __future__ import annotations

from prozor.matching.annotation import annotate_peptides
from prozor.matching.automaton import resolve_backend

__all__ = ["annotate_peptides", "resolve_backend"]
