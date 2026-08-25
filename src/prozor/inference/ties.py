"""Tie candidates and resolution for greedy protein inference."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class TieCandidate:
    """One indistinguishable protein group competing in a consequential tie."""

    proteins: tuple[str, ...]
    unexplained_peptides: frozenset[str]


type TieResolver = Callable[[Sequence[TieCandidate]], TieCandidate]


def resolve_current_tie(candidates: Sequence[TieCandidate]) -> TieCandidate:
    """Preserve Prozor's deterministic group-size and accession tie rule."""
    if not candidates:
        raise ValueError("at least one tie candidate is required")
    return min(candidates, key=lambda candidate: (-len(candidate.proteins), candidate.proteins))
