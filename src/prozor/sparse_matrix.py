"""Sparse peptide--protein topology used by protein inference."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Literal, Self, cast

import numpy as np
import numpy.typing as npt
from scipy import sparse

type Weighting = Literal["binary", "inverse"]
type IndexSelection = Sequence[int] | npt.NDArray[np.int_]

_VALID_WEIGHTINGS = frozenset({"binary", "inverse"})


@dataclass(slots=True)
class PeptideProteinMatrix:
    """Sparse peptide-by-protein matrix with explicit row and column labels."""

    matrix: sparse.csr_matrix
    peptides: list[str]
    proteins: list[str]

    def __post_init__(self) -> None:
        expected = (len(self.peptides), len(self.proteins))
        if self.matrix.shape != expected:
            raise ValueError(f"matrix shape {self.matrix.shape} does not match labels {expected}")

    @property
    def shape(self) -> tuple[int, int]:
        """Return ``(n_peptides, n_proteins)``."""
        rows, columns = self.matrix.shape
        return int(rows), int(columns)

    @property
    def n_peptides(self) -> int:
        """Return the number of peptide rows."""
        return len(self.peptides)

    @property
    def n_proteins(self) -> int:
        """Return the number of protein columns."""
        return len(self.proteins)

    @property
    def density(self) -> float:
        """Return the fraction of populated peptide--protein edges."""
        total = self.matrix.shape[0] * self.matrix.shape[1]
        return float(self.matrix.nnz / total) if total else 0.0

    def peptides_per_protein(self) -> npt.NDArray[np.int64]:
        """Count peptide edges incident to every protein."""
        return np.asarray(self.matrix.getnnz(axis=0), dtype=np.int64).ravel()

    def proteins_per_peptide(self) -> npt.NDArray[np.int64]:
        """Count protein edges incident to every peptide."""
        return np.asarray(self.matrix.getnnz(axis=1), dtype=np.int64).ravel()

    def proteotypic_peptides(self) -> list[str]:
        """Return peptides connected to exactly one protein."""
        counts = self.proteins_per_peptide()
        return [peptide for peptide, count in zip(self.peptides, counts, strict=True) if count == 1]

    def proteotypic_fraction(self) -> float:
        """Return the fraction of peptides connected to exactly one protein."""
        counts = self.proteins_per_peptide()
        if not counts.size:
            return 0.0
        proteotypic_count = sum(count == 1 for count in counts)
        return proteotypic_count / len(counts)

    @classmethod
    def from_edges(
        cls,
        edges: Iterable[tuple[str, str]],
        weighting: str = "binary",
    ) -> Self:
        """Build a matrix from ``(peptide, protein_id)`` edges."""
        validated_weighting = _validate_weighting(weighting)
        edge_list = list(edges)
        peptides = sorted({peptide for peptide, _protein in edge_list})
        proteins = sorted({protein for _peptide, protein in edge_list})
        peptide_indices = {peptide: index for index, peptide in enumerate(peptides)}
        protein_indices = {protein: index for index, protein in enumerate(proteins)}
        rows = [peptide_indices[peptide] for peptide, _protein in edge_list]
        columns = [protein_indices[protein] for _peptide, protein in edge_list]
        matrix = _matrix_from_coordinates(
            rows,
            columns,
            shape=(len(peptides), len(proteins)),
            weighting=validated_weighting,
        )
        return cls(matrix=matrix, peptides=peptides, proteins=proteins)

    def to_dense(self) -> npt.NDArray[np.float64]:
        """Return a dense floating-point representation."""
        return np.asarray(self.matrix.toarray(), dtype=np.float64)

    def subset_peptides(self, peptide_indices: IndexSelection) -> PeptideProteinMatrix:
        """Return a matrix containing selected peptide rows."""
        indices = np.asarray(peptide_indices, dtype=np.intp)
        new_matrix = sparse.csr_matrix(self.matrix[indices, :])
        new_peptides = [self.peptides[int(index)] for index in indices]
        return PeptideProteinMatrix(
            matrix=new_matrix,
            peptides=new_peptides,
            proteins=self.proteins.copy(),
        )

    def subset_proteins(self, protein_indices: IndexSelection) -> PeptideProteinMatrix:
        """Return a matrix containing selected protein columns."""
        indices = np.asarray(protein_indices, dtype=np.intp)
        new_matrix = sparse.csr_matrix(self.matrix[:, indices])
        new_proteins = [self.proteins[int(index)] for index in indices]
        return PeptideProteinMatrix(
            matrix=new_matrix,
            peptides=self.peptides.copy(),
            proteins=new_proteins,
        )

    def remove_zero_rows(self) -> PeptideProteinMatrix:
        """Remove peptide rows without protein edges."""
        indices = np.flatnonzero(self.proteins_per_peptide()).astype(np.int_)
        return self.subset_peptides(indices)

    def remove_zero_cols(self) -> PeptideProteinMatrix:
        """Remove protein columns without peptide edges."""
        indices = np.flatnonzero(self.peptides_per_protein()).astype(np.int_)
        return self.subset_proteins(indices)


def _matrix_from_coordinates(
    rows: Sequence[int],
    columns: Sequence[int],
    *,
    shape: tuple[int, int],
    weighting: Weighting,
) -> sparse.csr_matrix:
    matrix = sparse.csr_matrix(
        (np.ones(len(rows), dtype=np.float64), (rows, columns)),
        shape=shape,
    )
    matrix.sum_duplicates()
    if matrix.nnz:
        matrix.data.fill(1.0)
    if weighting == "binary":
        return matrix

    counts = np.asarray(matrix.getnnz(axis=1), dtype=np.float64).ravel()
    scales = np.ones_like(counts)
    nonzero = counts > 0
    scales[nonzero] = 1.0 / counts[nonzero]
    return sparse.csr_matrix(sparse.diags(scales, format="csr") @ matrix)


def _validate_weighting(weighting: str) -> Weighting:
    if weighting not in _VALID_WEIGHTINGS:
        raise ValueError(f"weighting must be one of {sorted(_VALID_WEIGHTINGS)}, got {weighting!r}")
    return cast(Weighting, weighting)
