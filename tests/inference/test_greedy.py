from __future__ import annotations

from collections.abc import Sequence

import pytest

from prozor.inference.greedy import ProteinGroup, greedy_parsimony
from prozor.inference.ties import TieCandidate
from prozor.matching.annotation import annotate_peptides


class _UnexpectedResolver:
    def __call__(self, candidates: Sequence[TieCandidate]) -> TieCandidate:
        raise AssertionError(f"resolver should not receive {candidates!r}")


def test_indistinguishable_proteins_are_grouped_without_tie_resolution() -> None:
    inference = greedy_parsimony(
        [
            ("PEP1", "PROT1"),
            ("PEP1", "PROT2"),
            ("PEP2", "PROT1"),
            ("PEP2", "PROT2"),
        ],
        resolve_tie=_UnexpectedResolver(),
    )
    assert inference.groups[0].proteins == ["PROT1", "PROT2"]


def test_disjoint_equal_candidates_are_ordered_without_tie_resolution() -> None:
    result = greedy_parsimony(
        [("PEP1", "PROT_B"), ("PEP2", "PROT_A")],
        resolve_tie=_UnexpectedResolver(),
    )
    assert [group.protein_id for group in result] == ["PROT_A", "PROT_B"]


def test_overlapping_nonidentical_tie_is_injected() -> None:
    received: list[tuple[TieCandidate, ...]] = []

    def choose_b(candidates: Sequence[TieCandidate]) -> TieCandidate:
        received.append(tuple(candidates))
        return next(candidate for candidate in candidates if candidate.proteins == ("B",))

    result = greedy_parsimony(
        [
            ("p1", "A"),
            ("p2", "A"),
            ("p1", "B"),
            ("p3", "B"),
            ("p2", "C"),
            ("p4", "C"),
        ],
        resolve_tie=choose_b,
        subsume=False,
    )
    assert received
    assert result.groups[0].protein_id == "B"
    assert result.n_groups == 2


def test_resolver_must_return_a_received_candidate() -> None:
    def invent_candidate(_candidates: Sequence[TieCandidate]) -> TieCandidate:
        return TieCandidate(("INVENTED",), frozenset({"p1"}))

    with pytest.raises(ValueError, match="must return one of"):
        greedy_parsimony(
            [("p1", "A"), ("p2", "A"), ("p1", "B"), ("p3", "B")],
            resolve_tie=invent_candidate,
        )


def test_duplicate_and_input_order_do_not_change_result() -> None:
    edges = [("p1", "B"), ("p2", "A"), ("p1", "B"), ("p3", "A")]
    assert greedy_parsimony(edges).to_dict() == greedy_parsimony(reversed(edges)).to_dict()


def test_subsumed_proteins_are_explicitly_configurable() -> None:
    edges = [("PEP1", "PROT1"), ("PEP2", "PROT1"), ("PEP1", "PROT2")]
    assert greedy_parsimony(edges, subsume=True).groups[0].proteins == ["PROT1", "PROT2"]
    assert greedy_parsimony(edges, subsume=False).groups[0].proteins == ["PROT1"]


def test_empty_edges_produce_empty_result() -> None:
    result = greedy_parsimony([])
    assert result.n_groups == 0
    assert result.n_peptides == 0
    assert result.n_proteins == 0


def test_full_matching_and_inference_workflow_uses_edges() -> None:
    annotations = annotate_peptides(
        ["PEPTIDE", "SEQUENCER", "UNIQUE"],
        {
            "PROT_A": "MKPEPTIDESEQUENCER",
            "PROT_B": "MSEQUENCEROTHER",
            "PROT_C": "MXUNIQUEPEPTIDEX",
        },
    )
    edges = {(match.peptide, match.protein_id) for match in annotations}
    inference = greedy_parsimony(edges)
    assert inference.n_peptides == 3
    assert inference.groups[0].protein_id == "PROT_A;PROT_B"


def test_current_r_prozor_fixture_after_group_canonicalization() -> None:
    edges = [
        ("PEP1", "PROT_A"),
        ("PEP1", "PROT_B"),
        ("PEP1", "PROT_C"),
        ("PEP2", "PROT_A"),
        ("PEP2", "PROT_B"),
        ("PEP2", "PROT_C"),
        ("PEP3", "PROT_D"),
    ]
    assert greedy_parsimony(edges, subsume=False).to_dict() == {
        "PEP1": "PROT_A;PROT_B;PROT_C",
        "PEP2": "PROT_A;PROT_B;PROT_C",
        "PEP3": "PROT_D",
    }


def test_protein_group_properties() -> None:
    group = ProteinGroup(proteins=["P1", "P2"], peptides=["PEP1", "PEP2"])
    assert group.protein_id == "P1;P2"
    assert group.n_proteins == 2
    assert group.n_peptides == 2
