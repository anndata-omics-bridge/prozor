from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass

import pytest

from prozor.matching.annotation import (
    AnnotationResult,
    ProteinSequenceRecord,
    annotate_peptides,
    annotate_peptides_streaming,
)
from prozor.matching.automaton import create_automaton, get_available_backends

PROTEINS = {
    "sp|P12345|PROT1": "MKWVTFISLLFSSAYSRGVFRRDTHK",
    "sp|P67890|PROT2": "MRGVFRRDTHKSEQ",
    "sp|Q11111|PROT3": "MXXUNIQUESEQXXX",
}


@dataclass(frozen=True, slots=True)
class _ProteinRecord:
    id: str
    sequence: str


def _protein_records() -> Iterator[ProteinSequenceRecord]:
    for protein_id, sequence in PROTEINS.items():
        yield _ProteinRecord(id=protein_id, sequence=sequence)


def _records(result: AnnotationResult) -> set[tuple[str, str, int, int]]:
    return {
        (annotation.peptide, annotation.protein_id, annotation.start, annotation.end)
        for annotation in result
    }


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_annotation_matches_mapping_and_streaming_input(backend: str) -> None:
    peptides = ["GVFRR", "DTHK", "UNIQUE"]
    mapped = annotate_peptides(peptides, PROTEINS, backend=backend)
    streamed = annotate_peptides_streaming(peptides, _protein_records(), backend=backend)
    assert _records(mapped) == _records(streamed)


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_annotation_reports_positions_within_each_sequence(backend: str) -> None:
    # THKMR spans the end of PROT1 and the start of PROT2 and must not match.
    peptides = ["GVFRR", "DTHK", "MK", "XXX", "THKMR"]
    expected = {
        (match.keyword, protein_id, match.start, match.end)
        for protein_id, sequence in PROTEINS.items()
        for match in create_automaton(peptides, backend=backend).find_all(sequence)
    }
    assert _records(annotate_peptides(peptides, PROTEINS, backend=backend)) == expected
    assert "THKMR" not in {peptide for peptide, *_ in expected}


def test_annotation_rejects_peptides_with_line_breaks() -> None:
    with pytest.raises(ValueError, match="line breaks"):
        annotate_peptides(["PEP\nTIDE"], PROTEINS)


def test_annotation_deduplicates_patterns_but_keeps_sites() -> None:
    result = annotate_peptides(["GVFRR", "GVFRR"], PROTEINS)
    assert result.peptides == {"GVFRR"}
    assert len(result) == 2


def test_annotation_handles_empty_and_unmatched_patterns() -> None:
    empty = annotate_peptides([], PROTEINS)
    unmatched = annotate_peptides(["ZZZZZ"], PROTEINS)
    assert len(empty) == 0
    assert len(unmatched) == 0
    assert empty.requested_backend == "auto"
    assert empty.resolved_backend in get_available_backends()


def test_tryptic_filtering() -> None:
    proteins = {"P1": "MKPEPTIDEARK"}
    result = annotate_peptides(["PEPTIDE", "MK"], proteins)
    assert len(result.filter_tryptic(proteins)) == 2
    assert len(result.filter_tryptic(proteins, allow_n_term=False)) == 1
