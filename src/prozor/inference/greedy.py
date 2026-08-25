"""Deterministic greedy-parsimony protein inference over unique edges."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass

from prozor.inference.ties import (
    TieCandidate,
    TieResolver,
    resolve_current_tie,
)

type PeptideProteinEdge = tuple[str, str]


@dataclass(slots=True)
class ProteinGroup:
    """Proteins selected together and the peptides assigned to them."""

    proteins: list[str]
    peptides: list[str]

    @property
    def protein_id(self) -> str:
        """Return the semicolon-joined protein group identifier."""
        return ";".join(self.proteins)

    @property
    def n_peptides(self) -> int:
        """Return the number of assigned peptides."""
        return len(self.peptides)

    @property
    def n_proteins(self) -> int:
        """Return the number of grouped proteins."""
        return len(self.proteins)


@dataclass(slots=True)
class GreedyResult:
    """Protein groups selected by greedy parsimony."""

    groups: list[ProteinGroup]

    def __len__(self) -> int:
        return len(self.groups)

    def __iter__(self) -> Iterator[ProteinGroup]:
        return iter(self.groups)

    @property
    def n_proteins(self) -> int:
        """Return the total number of selected and subsumed proteins."""
        return sum(group.n_proteins for group in self.groups)

    @property
    def n_groups(self) -> int:
        """Return the number of inferred protein groups."""
        return len(self.groups)

    @property
    def n_peptides(self) -> int:
        """Return the number of distinct assigned peptides."""
        return len({peptide for group in self.groups for peptide in group.peptides})

    def to_dict(self) -> dict[str, str]:
        """Return a peptide-to-protein-group mapping."""
        return {peptide: group.protein_id for group in self.groups for peptide in group.peptides}


@dataclass(frozen=True, slots=True)
class _Selection:
    winners: tuple[int, ...]
    covered_peptides: frozenset[int]


def greedy_parsimony(
    edges: Iterable[PeptideProteinEdge],
    *,
    resolve_tie: TieResolver = resolve_current_tie,
    subsume: bool = True,
) -> GreedyResult:
    """Select a deterministic parsimonious set of protein groups.

    Repeated ``(peptide, protein)`` edges collapse before inference. Proteins
    with identical remaining peptide evidence are grouped. Disjoint equal
    candidates are ordered deterministically, while overlapping non-identical
    candidates are passed to ``resolve_tie``.

    Args:
        edges: Peptide and protein identifier pairs.
        resolve_tie: Operation selecting one consequential tie candidate.
        subsume: Whether proteins supported by a subset of selected evidence
            remain in the selected group.

    Returns:
        Deterministically ordered inferred protein groups.
    """
    peptides, proteins, peptide_proteins, protein_peptides = _build_incidence(edges)
    active_peptides = set(range(len(peptides)))
    active_proteins = set(range(len(proteins)))
    counts = [len(protein_evidence) for protein_evidence in protein_peptides]
    groups: list[ProteinGroup] = []

    while active_peptides and active_proteins:
        selections = _select_winners(
            active_peptides,
            active_proteins,
            protein_peptides,
            counts,
            peptides,
            proteins,
            resolve_tie,
        )
        if not selections:
            break
        for selection in selections:
            subsumed = (
                _find_subsumed(
                    selection,
                    active_peptides,
                    active_proteins,
                    peptide_proteins,
                    protein_peptides,
                )
                if subsume
                else ()
            )
            group_indices = tuple(sorted((*selection.winners, *subsumed)))
            groups.append(
                ProteinGroup(
                    proteins=sorted(proteins[index] for index in group_indices),
                    peptides=sorted(peptides[index] for index in selection.covered_peptides),
                )
            )
            active_peptides.difference_update(selection.covered_peptides)
            active_proteins.difference_update(group_indices)
            _decrement_counts(
                selection.covered_peptides,
                active_proteins,
                peptide_proteins,
                counts,
            )

    groups.sort(
        key=lambda group: (
            -group.n_peptides,
            -group.n_proteins,
            tuple(group.proteins),
        )
    )
    return GreedyResult(groups=groups)


def _build_incidence(
    edges: Iterable[PeptideProteinEdge],
) -> tuple[list[str], list[str], list[frozenset[int]], list[frozenset[int]]]:
    unique_edges = set(edges)
    peptides = sorted({peptide for peptide, _protein in unique_edges})
    proteins = sorted({protein for _peptide, protein in unique_edges})
    peptide_indices = {peptide: index for index, peptide in enumerate(peptides)}
    protein_indices = {protein: index for index, protein in enumerate(proteins)}
    peptide_proteins: list[set[int]] = [set() for _peptide in peptides]
    protein_peptides: list[set[int]] = [set() for _protein in proteins]
    for peptide, protein in unique_edges:
        peptide_index = peptide_indices[peptide]
        protein_index = protein_indices[protein]
        peptide_proteins[peptide_index].add(protein_index)
        protein_peptides[protein_index].add(peptide_index)
    return (
        peptides,
        proteins,
        [frozenset(indices) for indices in peptide_proteins],
        [frozenset(indices) for indices in protein_peptides],
    )


def _select_winners(
    active_peptides: set[int],
    active_proteins: set[int],
    protein_peptides: Sequence[frozenset[int]],
    counts: Sequence[int],
    peptide_names: Sequence[str],
    protein_names: Sequence[str],
    resolve_tie: TieResolver,
) -> tuple[_Selection, ...]:
    ordered_active = sorted(active_proteins)
    if not ordered_active:
        return ()
    max_count = max(counts[index] for index in ordered_active)
    if max_count == 0:
        return ()
    max_proteins = [index for index in ordered_active if counts[index] == max_count]
    signature_groups: dict[frozenset[int], list[int]] = {}
    for index in max_proteins:
        signature = frozenset(protein_peptides[index] & active_peptides)
        signature_groups.setdefault(signature, []).append(index)

    grouped = tuple(
        tuple(indices)
        for _signature, indices in sorted(
            signature_groups.items(),
            key=lambda item: tuple(protein_names[index] for index in item[1]),
        )
    )
    signatures = tuple(
        frozenset(protein_peptides[indices[0]] & active_peptides) for indices in grouped
    )
    components = _overlap_components(signatures)
    selections = [
        _selection_for_component(
            component,
            grouped,
            signatures,
            peptide_names,
            protein_names,
            resolve_tie,
        )
        for component in components
    ]
    return tuple(sorted(selections, key=lambda selection: _selection_key(selection, protein_names)))


def _selection_key(
    selection: _Selection,
    protein_names: Sequence[str],
) -> tuple[int, tuple[str, ...]]:
    return (
        -len(selection.winners),
        tuple(protein_names[index] for index in selection.winners),
    )


def _overlap_components(signatures: Sequence[frozenset[int]]) -> tuple[tuple[int, ...], ...]:
    peptide_groups: dict[int, list[int]] = {}
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


def _selection_for_component(
    component: tuple[int, ...],
    groups: Sequence[tuple[int, ...]],
    signatures: Sequence[frozenset[int]],
    peptide_names: Sequence[str],
    protein_names: Sequence[str],
    resolve_tie: TieResolver,
) -> _Selection:
    component_groups = tuple(groups[index] for index in component)
    component_signatures = tuple(signatures[index] for index in component)
    winners = (
        component_groups[0]
        if len(component_groups) == 1
        else _resolve_overlapping_tie(
            component_groups,
            component_signatures,
            peptide_names,
            protein_names,
            resolve_tie,
        )
    )
    winner_index = component_groups.index(winners)
    return _Selection(winners=winners, covered_peptides=component_signatures[winner_index])


def _resolve_overlapping_tie(
    groups: Sequence[tuple[int, ...]],
    signatures: Sequence[frozenset[int]],
    peptide_names: Sequence[str],
    protein_names: Sequence[str],
    resolve_tie: TieResolver,
) -> tuple[int, ...]:
    candidates = tuple(
        TieCandidate(
            proteins=tuple(protein_names[index] for index in group),
            unexplained_peptides=frozenset(peptide_names[index] for index in signature),
        )
        for group, signature in zip(groups, signatures, strict=True)
    )
    selected = resolve_tie(candidates)
    try:
        selected_index = candidates.index(selected)
    except ValueError as error:
        raise ValueError("tie resolver must return one of the candidates it received") from error
    return groups[selected_index]


def _find_subsumed(
    selection: _Selection,
    active_peptides: set[int],
    active_proteins: set[int],
    peptide_proteins: Sequence[frozenset[int]],
    protein_peptides: Sequence[frozenset[int]],
) -> tuple[int, ...]:
    candidates = {
        protein for peptide in selection.covered_peptides for protein in peptide_proteins[peptide]
    }
    candidates.difference_update(selection.winners)
    candidates.intersection_update(active_proteins)
    return tuple(
        protein
        for protein in sorted(candidates)
        if protein_peptides[protein] & active_peptides <= selection.covered_peptides
    )


def _decrement_counts(
    removed_peptides: frozenset[int],
    active_proteins: set[int],
    peptide_proteins: Sequence[frozenset[int]],
    counts: list[int],
) -> None:
    for peptide in removed_peptides:
        for protein in peptide_proteins[peptide] & active_proteins:
            counts[protein] -= 1
