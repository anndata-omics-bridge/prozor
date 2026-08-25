"""Benchmark-local equivalent implementations of greedy parsimony."""

from __future__ import annotations

import heapq
from collections.abc import Callable, Iterable, Mapping, Sequence

from prozor.inference.greedy import GreedyResult, PeptideProteinEdge, ProteinGroup, greedy_parsimony
from prozor.inference.ties import TieCandidate, TieResolver

type GreedyImplementation = Callable[[Iterable[PeptideProteinEdge], TieResolver], GreedyResult]


def production_incremental(
    edges: Iterable[PeptideProteinEdge],
    resolve_tie: TieResolver,
) -> GreedyResult:
    """Run the production integer-indexed incremental implementation."""
    return greedy_parsimony(edges, resolve_tie=resolve_tie)


def direct_name_sets(
    edges: Iterable[PeptideProteinEdge],
    resolve_tie: TieResolver,
) -> GreedyResult:
    """Run a readable oracle that recomputes every active protein score."""
    protein_peptides, peptide_proteins = _name_incidence(edges)
    active_peptides = set(peptide_proteins)
    active_proteins = set(protein_peptides)
    groups: list[ProteinGroup] = []
    while active_peptides and active_proteins:
        evidence = {
            protein: protein_peptides[protein] & active_peptides for protein in active_proteins
        }
        max_count = max(map(len, evidence.values()), default=0)
        if max_count == 0:
            break
        selections = _choose_name_groups(
            {
                protein: peptides
                for protein, peptides in evidence.items()
                if len(peptides) == max_count
            },
            resolve_tie,
        )
        for winners, covered in selections:
            active_evidence = {
                protein: protein_peptides[protein] & active_peptides for protein in active_proteins
            }
            subsumed = _name_subsumed(
                winners,
                covered,
                active_evidence,
                peptide_proteins,
                active_proteins,
            )
            selected = tuple(sorted((*winners, *subsumed)))
            groups.append(ProteinGroup(proteins=list(selected), peptides=sorted(covered)))
            active_peptides.difference_update(covered)
            active_proteins.difference_update(selected)
    return _ordered_result(groups)


def lazy_priority_queue(
    edges: Iterable[PeptideProteinEdge],
    resolve_tie: TieResolver,
) -> GreedyResult:
    """Run a name-set implementation whose maximum scores come from a lazy heap."""
    protein_peptides, peptide_proteins = _name_incidence(edges)
    active_peptides = set(peptide_proteins)
    active_proteins = set(protein_peptides)
    counts = {protein: len(peptides) for protein, peptides in protein_peptides.items()}
    heap = [(-count, protein) for protein, count in counts.items()]
    heapq.heapify(heap)
    groups: list[ProteinGroup] = []
    while active_peptides and active_proteins:
        max_proteins = _pop_current_maximum(heap, counts, active_proteins)
        if not max_proteins:
            break
        evidence = {
            protein: protein_peptides[protein] & active_peptides for protein in max_proteins
        }
        for winners, covered in _choose_name_groups(evidence, resolve_tie):
            active_evidence = {
                protein: protein_peptides[protein] & active_peptides for protein in active_proteins
            }
            subsumed = _name_subsumed(
                winners,
                covered,
                active_evidence,
                peptide_proteins,
                active_proteins,
            )
            selected = tuple(sorted((*winners, *subsumed)))
            groups.append(ProteinGroup(proteins=list(selected), peptides=sorted(covered)))
            active_peptides.difference_update(covered)
            active_proteins.difference_update(selected)
            affected = {
                protein
                for peptide in covered
                for protein in peptide_proteins[peptide] & active_proteins
            }
            for protein in affected:
                counts[protein] -= len(protein_peptides[protein] & covered)
                heapq.heappush(heap, (-counts[protein], protein))
    return _ordered_result(groups)


def _name_incidence(
    edges: Iterable[PeptideProteinEdge],
) -> tuple[dict[str, frozenset[str]], dict[str, frozenset[str]]]:
    protein_peptides: dict[str, set[str]] = {}
    peptide_proteins: dict[str, set[str]] = {}
    for peptide, protein in set(edges):
        protein_peptides.setdefault(protein, set()).add(peptide)
        peptide_proteins.setdefault(peptide, set()).add(protein)
    return (
        {protein: frozenset(peptides) for protein, peptides in protein_peptides.items()},
        {peptide: frozenset(proteins) for peptide, proteins in peptide_proteins.items()},
    )


def _choose_name_groups(
    evidence: Mapping[str, frozenset[str]],
    resolve_tie: TieResolver,
) -> tuple[tuple[tuple[str, ...], frozenset[str]], ...]:
    signature_groups: dict[frozenset[str], list[str]] = {}
    for protein, peptides in evidence.items():
        signature_groups.setdefault(peptides, []).append(protein)
    groups = tuple(
        (tuple(sorted(proteins)), peptides)
        for peptides, proteins in sorted(
            signature_groups.items(),
            key=lambda item: tuple(sorted(item[1])),
        )
    )
    signatures = tuple(peptides for _proteins, peptides in groups)
    selections = [
        _choose_name_component(component, groups, resolve_tie)
        for component in _overlap_components(signatures)
    ]
    return tuple(sorted(selections, key=lambda selection: (-len(selection[0]), selection[0])))


def _overlap_components(signatures: Sequence[frozenset[str]]) -> tuple[tuple[int, ...], ...]:
    peptide_groups: dict[str, list[int]] = {}
    for group_index, signature in enumerate(signatures):
        for peptide in signature:
            peptide_groups.setdefault(peptide, []).append(group_index)
    unseen = set(range(len(signatures)))
    components: list[tuple[int, ...]] = []
    while unseen:
        pending = [min(unseen)]
        component: set[int] = set()
        while pending:
            group_index = pending.pop()
            if group_index not in unseen:
                continue
            unseen.remove(group_index)
            component.add(group_index)
            for peptide in signatures[group_index]:
                pending.extend(peptide_groups[peptide])
        components.append(tuple(sorted(component)))
    return tuple(components)


def _choose_name_component(
    component: tuple[int, ...],
    groups: Sequence[tuple[tuple[str, ...], frozenset[str]]],
    resolve_tie: TieResolver,
) -> tuple[tuple[str, ...], frozenset[str]]:
    choices = tuple(groups[index] for index in component)
    if len(choices) == 1:
        return choices[0]
    candidates = tuple(TieCandidate(proteins, peptides) for proteins, peptides in choices)
    selected = resolve_tie(candidates)
    try:
        return choices[candidates.index(selected)]
    except ValueError as error:
        raise ValueError("tie resolver must return one of the candidates it received") from error


def _ordered_result(groups: list[ProteinGroup]) -> GreedyResult:
    groups.sort(key=lambda group: (-group.n_peptides, -group.n_proteins, tuple(group.proteins)))
    return GreedyResult(groups)


def _name_subsumed(
    winners: tuple[str, ...],
    covered: frozenset[str],
    evidence: Mapping[str, frozenset[str]],
    peptide_proteins: Mapping[str, frozenset[str]],
    active_proteins: set[str],
) -> tuple[str, ...]:
    candidates = {
        protein for peptide in covered for protein in peptide_proteins[peptide]
    } & active_proteins
    candidates.difference_update(winners)
    return tuple(protein for protein in sorted(candidates) if evidence[protein] <= covered)


def _pop_current_maximum(
    heap: list[tuple[int, str]],
    counts: Mapping[str, int],
    active_proteins: set[str],
) -> tuple[str, ...]:
    while heap:
        negative_count, protein = heap[0]
        if protein in active_proteins and -negative_count == counts[protein]:
            break
        heapq.heappop(heap)
    if not heap or heap[0][0] == 0:
        return ()
    maximum = heap[0][0]
    proteins: set[str] = set()
    retained: list[tuple[int, str]] = []
    while heap and heap[0][0] == maximum:
        entry = heapq.heappop(heap)
        negative_count, protein = entry
        if protein in active_proteins and -negative_count == counts[protein]:
            proteins.add(protein)
            retained.append(entry)
    for entry in retained:
        heapq.heappush(heap, entry)
    return tuple(sorted(proteins))
