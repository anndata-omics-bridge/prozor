"""Explore tie scores and quantitative-profile coherence without changing production."""

from __future__ import annotations

import json
import math
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

from prozor.inference.ties import TieCandidate, resolve_current_tie

type Profile = tuple[float | None, ...]
type StudyResolver = Callable[[TieStudyCase], TieCandidate]


@dataclass(frozen=True, slots=True)
class TieStudyCase:
    """Synthetic tie with independently known truth and score evidence."""

    name: str
    candidates: tuple[TieCandidate, ...]
    expected: tuple[str, ...] | None
    peptide_degeneracy: Mapping[str, int]
    sequence_coverage: Mapping[str, float]
    external_scores: Mapping[str, float]


@dataclass(frozen=True, slots=True)
class ProfileScenario:
    """Peptide profiles and anchors for a quantitative-coherence experiment."""

    name: str
    profiles: Mapping[str, Profile]
    protein_peptides: Mapping[str, frozenset[str]]
    interpretation: str


def main() -> None:
    """Write the deterministic exploratory report."""
    resolvers: tuple[tuple[str, StudyResolver], ...] = (
        ("current", lambda case: resolve_current_tie(case.candidates)),
        ("inverse_degeneracy", _by_inverse_degeneracy),
        ("sequence_coverage", _by_sequence_coverage),
        ("external_score", _by_external_score),
    )
    cases = _tie_cases()
    tie_results = [_tie_result(case, resolvers) for case in cases]
    profile_results = [_profile_result(scenario) for scenario in _profile_scenarios()]
    output = Path(__file__).parent / "results" / "tie_scoring_after.json"
    output.write_text(
        json.dumps(
            {
                "status": "exploratory; no production policy changed",
                "tie_scoring": tie_results,
                "accuracy_on_known_truth": _accuracy(cases, resolvers),
                "quantitative_coherence": profile_results,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )


def _tie_result(
    case: TieStudyCase,
    resolvers: Sequence[tuple[str, StudyResolver]],
) -> dict[str, object]:
    return {
        "case": case.name,
        "expected": case.expected,
        "selections": {name: resolver(case).proteins for name, resolver in resolvers},
        "single_peptide_dropout_stability": {
            name: _dropout_stability(case, resolver) for name, resolver in resolvers
        },
    }


def _accuracy(
    cases: Sequence[TieStudyCase],
    resolvers: Sequence[tuple[str, StudyResolver]],
) -> dict[str, dict[str, int]]:
    known = tuple(case for case in cases if case.expected is not None)
    return {
        name: {
            "correct": sum(resolver(case).proteins == case.expected for case in known),
            "cases": len(known),
        }
        for name, resolver in resolvers
    }


def _dropout_stability(case: TieStudyCase, resolver: StudyResolver) -> dict[str, int]:
    baseline = resolver(case).proteins
    stable = 0
    trials = 0
    for candidate_index, candidate in enumerate(case.candidates):
        for peptide in candidate.unexplained_peptides:
            candidates = list(case.candidates)
            candidates[candidate_index] = TieCandidate(
                proteins=candidate.proteins,
                unexplained_peptides=candidate.unexplained_peptides - {peptide},
            )
            perturbed = TieStudyCase(
                name=case.name,
                candidates=tuple(candidates),
                expected=case.expected,
                peptide_degeneracy=case.peptide_degeneracy,
                sequence_coverage=case.sequence_coverage,
                external_scores=case.external_scores,
            )
            stable += resolver(perturbed).proteins == baseline
            trials += 1
    return {"stable": stable, "trials": trials}


def _by_inverse_degeneracy(case: TieStudyCase) -> TieCandidate:
    return _highest(
        case.candidates,
        lambda candidate: sum(
            1 / case.peptide_degeneracy[peptide] for peptide in candidate.unexplained_peptides
        ),
    )


def _by_sequence_coverage(case: TieStudyCase) -> TieCandidate:
    return _highest(
        case.candidates,
        lambda candidate: max(case.sequence_coverage[protein] for protein in candidate.proteins),
    )


def _by_external_score(case: TieStudyCase) -> TieCandidate:
    return _highest(
        case.candidates,
        lambda candidate: max(case.external_scores[protein] for protein in candidate.proteins),
    )


def _highest(
    candidates: Sequence[TieCandidate],
    score: Callable[[TieCandidate], float],
) -> TieCandidate:
    return min(candidates, key=lambda candidate: (-score(candidate), candidate.proteins))


def _tie_cases() -> tuple[TieStudyCase, ...]:
    a = TieCandidate(("A",), frozenset({"shared", "a_evidence"}))
    b = TieCandidate(("B",), frozenset({"shared", "b_evidence"}))
    return (
        TieStudyCase(
            name="informative_unique_evidence",
            candidates=(a, b),
            expected=("A",),
            peptide_degeneracy={"shared": 2, "a_evidence": 1, "b_evidence": 4},
            sequence_coverage={"A": 0.35, "B": 0.30},
            external_scores={"A": 0.9, "B": 0.7},
        ),
        TieStudyCase(
            name="coverage_favors_b",
            candidates=(a, b),
            expected=("B",),
            peptide_degeneracy={"shared": 2, "a_evidence": 2, "b_evidence": 2},
            sequence_coverage={"A": 0.20, "B": 0.55},
            external_scores={"A": 0.6, "B": 0.8},
        ),
        TieStudyCase(
            name="deliberately_ambiguous",
            candidates=(a, b),
            expected=None,
            peptide_degeneracy={"shared": 2, "a_evidence": 2, "b_evidence": 2},
            sequence_coverage={"A": 0.30, "B": 0.30},
            external_scores={"A": 0.5, "B": 0.5},
        ),
    )


def _profile_result(scenario: ProfileScenario) -> dict[str, object]:
    clusters = _correlation_clusters(scenario.profiles, threshold=0.8)
    anchors = {
        protein: sorted(
            peptide
            for peptide in peptides
            if sum(peptide in other for other in scenario.protein_peptides.values()) == 1
        )
        for protein, peptides in scenario.protein_peptides.items()
    }
    return {
        "case": scenario.name,
        "clusters": clusters,
        "unique_anchors": anchors,
        "interpretation": scenario.interpretation,
    }


def _correlation_clusters(
    profiles: Mapping[str, Profile],
    *,
    threshold: float,
) -> list[list[str]]:
    adjacency: dict[str, set[str]] = {peptide: set() for peptide in profiles}
    names = sorted(profiles)
    for left_index, left in enumerate(names):
        for right in names[left_index + 1 :]:
            correlation = _pearson(profiles[left], profiles[right])
            if correlation is not None and correlation >= threshold:
                adjacency[left].add(right)
                adjacency[right].add(left)
    clusters: list[list[str]] = []
    unseen = set(names)
    while unseen:
        pending = [min(unseen)]
        cluster: set[str] = set()
        while pending:
            peptide = pending.pop()
            if peptide not in unseen:
                continue
            unseen.remove(peptide)
            cluster.add(peptide)
            pending.extend(adjacency[peptide])
        clusters.append(sorted(cluster))
    return clusters


def _pearson(left: Profile, right: Profile) -> float | None:
    pairs = [(x, y) for x, y in zip(left, right, strict=True) if x is not None and y is not None]
    if len(pairs) < 3:
        return None
    left_mean = sum(x for x, _y in pairs) / len(pairs)
    right_mean = sum(y for _x, y in pairs) / len(pairs)
    numerator = sum((x - left_mean) * (y - right_mean) for x, y in pairs)
    left_scale = math.sqrt(sum((x - left_mean) ** 2 for x, _y in pairs))
    right_scale = math.sqrt(sum((y - right_mean) ** 2 for _x, y in pairs))
    scale = left_scale * right_scale
    return numerator / scale if scale else None


def _profile_scenarios() -> tuple[ProfileScenario, ...]:
    return (
        ProfileScenario(
            name="homogeneous_protein",
            profiles={
                "a1": (1.0, 2.0, 4.0, 8.0),
                "a2": (1.1, 2.1, 4.2, 8.1),
                "a3": (0.9, 1.9, 3.8, 7.9),
            },
            protein_peptides={"A": frozenset({"a1", "a2", "a3"})},
            interpretation="One coherent cluster is consistent with one quantitative entity.",
        ),
        ProfileScenario(
            name="two_proteins_with_unique_support",
            profiles={
                "a_unique": (1.0, 2.0, 4.0, 8.0),
                "shared": (1.1, 2.1, 4.2, 8.1),
                "b_unique": (8.0, 4.0, 2.0, 1.0),
            },
            protein_peptides={
                "A": frozenset({"a_unique", "shared"}),
                "B": frozenset({"shared", "b_unique"}),
            },
            interpretation=(
                "Unique anchors support two proteins; the shared peptide profile may reveal "
                "which signal dominates but is not needed to establish both."
            ),
        ),
        ProfileScenario(
            name="shared_cluster_with_b_anchor",
            profiles={
                "p1": (1.0, 2.0, 4.0, 8.0),
                "p2": (1.1, 2.2, 4.1, 8.2),
                "p3": (0.9, 1.8, 3.9, 7.8),
                "p4": (8.0, 4.0, 2.0, 1.0),
                "p5": (7.8, 3.9, 2.1, 1.1),
                "b_unique": (8.2, 4.1, 1.9, 0.9),
            },
            protein_peptides={
                "A": frozenset({"p1", "p2", "p3", "p4", "p5"}),
                "B": frozenset({"p4", "p5", "b_unique"}),
            },
            interpretation=(
                "The second cluster has a B-unique anchor, so it is evidence for B rather than "
                "merely heterogeneity within A."
            ),
        ),
        ProfileScenario(
            name="shared_cluster_without_b_anchor",
            profiles={
                "p1": (1.0, 2.0, 4.0, 8.0),
                "p2": (1.1, 2.2, 4.1, 8.2),
                "p3": (0.9, 1.8, 3.9, 7.8),
                "p4": (8.0, 4.0, 2.0, 1.0),
                "p5": (7.8, 3.9, 2.1, 1.1),
            },
            protein_peptides={
                "A": frozenset({"p1", "p2", "p3", "p4", "p5"}),
                "B": frozenset({"p4", "p5"}),
            },
            interpretation=(
                "Two clusters diagnose disagreement, but without a B-unique anchor they do not "
                "identify B; proteoforms or interference remain alternatives."
            ),
        ),
        ProfileScenario(
            name="one_protein_two_proteoforms",
            profiles={
                "form_1a": (1.0, 2.0, 4.0, 8.0),
                "form_1b": (1.1, 2.1, 4.2, 8.2),
                "form_2a": (8.0, 4.0, 2.0, 1.0),
                "form_2b": (8.2, 4.1, 2.1, 1.1),
            },
            protein_peptides={"A": frozenset({"form_1a", "form_1b", "form_2a", "form_2b"})},
            interpretation=(
                "Two clusters can occur within one protein and therefore cannot alone prove a "
                "second protein assignment."
            ),
        ),
        ProfileScenario(
            name="technical_interference",
            profiles={
                "a1": (1.0, 2.0, 4.0, 8.0),
                "a2": (1.1, 2.1, 4.1, 8.2),
                "interfered": (8.0, 1.0, 7.0, 2.0),
            },
            protein_peptides={"A": frozenset({"a1", "a2", "interfered"})},
            interpretation=(
                "An isolated incoherent peptide is a quality warning; it supplies no protein "
                "identity evidence."
            ),
        ),
        ProfileScenario(
            name="structured_missingness",
            profiles={
                "a1": (1.0, None, 4.0, 8.0, 16.0),
                "a2": (1.1, 2.0, None, 8.2, 16.1),
                "a3": (None, 2.1, 4.2, 8.1, 15.9),
            },
            protein_peptides={"A": frozenset({"a1", "a2", "a3"})},
            interpretation="Pairwise-complete correlation remains diagnostic, not assignment proof.",
        ),
    )


if __name__ == "__main__":
    main()
