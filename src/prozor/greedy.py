"""Deterministic greedy-parsimony protein inference."""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
from scipy import sparse

from prozor.sparse_matrix import PeptideProteinMatrix


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
    peptide_protein: PeptideProteinMatrix,
    subsume: bool = True,
) -> GreedyResult:
    """Select a deterministic parsimonious set of protein groups.

    Proteins with identical remaining peptide evidence are grouped. When
    ``subsume`` is true, proteins whose remaining peptides are a subset of the
    selected evidence are retained in that group instead of being discarded.
    """
    peptide_proteins, protein_peptides = _build_incidence(peptide_protein.matrix)
    active_peptides = set(range(peptide_protein.n_peptides))
    active_proteins = set(range(peptide_protein.n_proteins))
    counts = np.asarray([len(peptides) for peptides in protein_peptides], dtype=np.int64)
    groups: list[ProteinGroup] = []

    while active_peptides and active_proteins:
        selection = _select_winners(
            active_peptides,
            active_proteins,
            protein_peptides,
            counts,
            peptide_protein.proteins,
        )
        if selection is None:
            break
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
                proteins=sorted(peptide_protein.proteins[index] for index in group_indices),
                peptides=sorted(
                    peptide_protein.peptides[index] for index in selection.covered_peptides
                ),
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

    return GreedyResult(groups=groups)


def _build_incidence(
    matrix: sparse.csr_matrix,
) -> tuple[list[frozenset[int]], list[frozenset[int]]]:
    topology = matrix.astype(bool).astype(np.int8).tocsr()
    topology.sum_duplicates()
    peptide_proteins = [
        frozenset(
            int(index)
            for index in topology.indices[topology.indptr[row] : topology.indptr[row + 1]]
        )
        for row in range(topology.shape[0])
    ]
    transposed = topology.transpose().tocsr()
    protein_peptides = [
        frozenset(
            int(index)
            for index in transposed.indices[
                transposed.indptr[column] : transposed.indptr[column + 1]
            ]
        )
        for column in range(transposed.shape[0])
    ]
    return peptide_proteins, protein_peptides


def _select_winners(
    active_peptides: set[int],
    active_proteins: set[int],
    protein_peptides: Sequence[frozenset[int]],
    counts: npt.NDArray[np.int64],
    protein_names: Sequence[str],
) -> _Selection | None:
    ordered_active = sorted(active_proteins)
    if not ordered_active:
        return None
    max_count = max(int(counts[index]) for index in ordered_active)
    if max_count == 0:
        return None
    candidates = [index for index in ordered_active if counts[index] == max_count]
    signature_groups: dict[frozenset[int], list[int]] = {}
    for index in candidates:
        signature = frozenset(protein_peptides[index] & active_peptides)
        signature_groups.setdefault(signature, []).append(index)
    winner_groups = sorted(
        signature_groups.values(),
        key=lambda group: (
            -len(group),
            tuple(protein_names[index] for index in group),
        ),
    )
    winners = tuple(winner_groups[0])
    covered = frozenset(protein_peptides[winners[0]] & active_peptides)
    return _Selection(winners=winners, covered_peptides=covered)


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
    counts: npt.NDArray[np.int64],
) -> None:
    for peptide in removed_peptides:
        for protein in peptide_proteins[peptide] & active_proteins:
            counts[protein] -= 1
