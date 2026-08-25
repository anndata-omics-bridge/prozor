from __future__ import annotations

import pytest

from prozor.inference.ties import TieCandidate, resolve_current_tie


def test_current_tie_prefers_larger_identical_evidence_group() -> None:
    candidates = (
        TieCandidate(("A", "B"), frozenset({"p1", "p2"})),
        TieCandidate(("C",), frozenset({"p1", "p3"})),
    )
    assert resolve_current_tie(candidates) == candidates[0]


def test_current_tie_uses_accession_order_after_group_size() -> None:
    candidates = (
        TieCandidate(("B",), frozenset({"p1", "p2"})),
        TieCandidate(("A",), frozenset({"p1", "p3"})),
    )
    assert resolve_current_tie(candidates) == candidates[1]


def test_current_tie_requires_candidates() -> None:
    with pytest.raises(ValueError, match="at least one"):
        resolve_current_tie(())
