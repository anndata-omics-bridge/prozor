from __future__ import annotations

from prozor.annotate import AnnotationResult, PeptideAnnotation, annotate_peptides
from prozor.greedy import ProteinGroup, greedy_parsimony
from prozor.sparse_matrix import PeptideProteinMatrix


def _matrix(edges: list[tuple[str, str]]) -> PeptideProteinMatrix:
    return PeptideProteinMatrix.from_edges(edges)


def test_indistinguishable_proteins_are_grouped() -> None:
    inference = greedy_parsimony(
        _matrix(
            [
                ("PEP1", "PROT1"),
                ("PEP1", "PROT2"),
                ("PEP2", "PROT1"),
                ("PEP2", "PROT2"),
            ]
        )
    )
    assert inference.n_groups == 1
    assert inference.groups[0].proteins == ["PROT1", "PROT2"]


def test_greedy_selection_is_deterministic() -> None:
    matrix = _matrix(
        [
            ("PEP1", "PROT_B"),
            ("PEP2", "PROT_B"),
            ("PEP3", "PROT_A"),
            ("PEP4", "PROT_A"),
        ]
    )
    first = greedy_parsimony(matrix)
    second = greedy_parsimony(matrix)
    assert first.to_dict() == second.to_dict()
    assert first.groups[0].protein_id == "PROT_A"


def test_subsumed_proteins_are_explicitly_configurable() -> None:
    matrix = _matrix(
        [
            ("PEP1", "PROT1"),
            ("PEP2", "PROT1"),
            ("PEP1", "PROT2"),
        ]
    )
    assert greedy_parsimony(matrix, subsume=True).groups[0].proteins == ["PROT1", "PROT2"]
    assert greedy_parsimony(matrix, subsume=False).groups[0].proteins == ["PROT1"]


def test_inverse_weighting_does_not_change_inference_topology() -> None:
    annotations = AnnotationResult(
        annotations=[
            PeptideAnnotation("SHARED", "PROT1", 0, 6),
            PeptideAnnotation("SHARED", "PROT2", 0, 6),
        ]
    )
    binary = greedy_parsimony(annotations.to_sparse_matrix())
    inverse = greedy_parsimony(annotations.to_sparse_matrix(weighting="inverse"))
    assert binary.to_dict() == inverse.to_dict()
    assert inverse.groups[0].proteins == ["PROT1", "PROT2"]


def test_empty_matrix_produces_empty_result() -> None:
    matrix = AnnotationResult(annotations=[]).to_sparse_matrix()
    result = greedy_parsimony(matrix)
    assert result.n_groups == 0
    assert result.n_peptides == 0
    assert result.n_proteins == 0


def test_full_matching_and_inference_workflow() -> None:
    annotations = annotate_peptides(
        ["PEPTIDE", "SEQUENCER", "UNIQUE"],
        {
            "PROT_A": "MKPEPTIDESEQUENCER",
            "PROT_B": "MSEQUENCEROTHER",
            "PROT_C": "MXUNIQUEPEPTIDEX",
        },
    )
    inference = greedy_parsimony(annotations.to_sparse_matrix())
    assert inference.n_peptides == 3
    assert inference.groups[0].protein_id == "PROT_A;PROT_B"


def test_current_r_prozor_fixture_after_group_canonicalization() -> None:
    """Match a fixture generated from the local R prozor 0.3.4 source.

    R returned ``PROT_B;PROT_C;PROT_A`` for the equivalent group. Python
    deliberately canonicalizes the same membership lexicographically.
    """
    matrix = _matrix(
        [
            ("PEP1", "PROT_A"),
            ("PEP1", "PROT_B"),
            ("PEP1", "PROT_C"),
            ("PEP2", "PROT_A"),
            ("PEP2", "PROT_B"),
            ("PEP2", "PROT_C"),
            ("PEP3", "PROT_D"),
        ]
    )
    assert greedy_parsimony(matrix, subsume=False).to_dict() == {
        "PEP1": "PROT_A;PROT_B;PROT_C",
        "PEP2": "PROT_A;PROT_B;PROT_C",
        "PEP3": "PROT_D",
    }


def test_protein_group_properties() -> None:
    group = ProteinGroup(proteins=["P1", "P2"], peptides=["PEP1", "PEP2"])
    assert group.protein_id == "P1;P2"
    assert group.n_proteins == 2
    assert group.n_peptides == 2
