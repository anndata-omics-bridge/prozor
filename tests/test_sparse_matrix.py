from __future__ import annotations

import numpy as np
import pytest

from prozor.annotate import AnnotationResult, PeptideAnnotation
from prozor.sparse_matrix import PeptideProteinMatrix


def _annotations() -> AnnotationResult:
    return AnnotationResult(
        annotations=[
            PeptideAnnotation("PEP1", "PROT1", 0, 4),
            PeptideAnnotation("PEP1", "PROT1", 5, 9),
            PeptideAnnotation("PEP1", "PROT2", 0, 4),
            PeptideAnnotation("PEP2", "PROT1", 10, 14),
            PeptideAnnotation("PEP3", "PROT3", 0, 4),
        ]
    )


def test_matrix_is_binary_and_deduplicates_repeated_sites() -> None:
    matrix = _annotations().to_sparse_matrix()
    assert matrix.shape == (3, 3)
    assert matrix.matrix.nnz == 4
    assert set(matrix.to_dense().ravel()) <= {0.0, 1.0}


def test_topology_counts_do_not_depend_on_inverse_weights() -> None:
    matrix = _annotations().to_sparse_matrix(weighting="inverse")
    np.testing.assert_array_equal(matrix.proteins_per_peptide(), [2, 1, 1])
    np.testing.assert_array_equal(matrix.peptides_per_protein(), [2, 1, 1])
    assert matrix.proteotypic_peptides() == ["PEP2", "PEP3"]
    assert matrix.proteotypic_fraction() == pytest.approx(2 / 3)
    np.testing.assert_allclose(matrix.to_dense().sum(axis=1), [1.0, 1.0, 1.0])


def test_empty_matrix_has_defined_density() -> None:
    matrix = AnnotationResult(annotations=[]).to_sparse_matrix()
    assert matrix.shape == (0, 0)
    assert matrix.density == 0.0
    assert matrix.proteotypic_fraction() == 0.0


def test_invalid_weighting_fails() -> None:
    with pytest.raises(ValueError, match="weighting must be one of"):
        _annotations().to_sparse_matrix(weighting="unsupported")


def test_subsetting_and_zero_removal_preserve_labels() -> None:
    matrix = _annotations().to_sparse_matrix()
    peptide_subset = matrix.subset_peptides(np.array([0, 2], dtype=np.int_))
    protein_subset = matrix.subset_proteins(np.array([0, 2], dtype=np.int_))
    assert peptide_subset.peptides == ["PEP1", "PEP3"]
    assert protein_subset.proteins == ["PROT1", "PROT3"]
    assert matrix.remove_zero_rows().peptides == matrix.peptides
    assert matrix.remove_zero_cols().proteins == matrix.proteins


def test_shape_must_match_labels() -> None:
    with pytest.raises(ValueError, match="does not match labels"):
        PeptideProteinMatrix(
            matrix=_annotations().to_sparse_matrix().matrix,
            peptides=["only-one"],
            proteins=["PROT1", "PROT2", "PROT3"],
        )
