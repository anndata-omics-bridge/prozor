from __future__ import annotations

import pytest

from prozor.matching.annotation import (
    AnnotationResult,
    annotate_peptides,
    annotate_peptides_streaming,
)
from prozor.matching.automaton import get_available_backends

PROTEINS = {
    "sp|P12345|PROT1": "MKWVTFISLLFSSAYSRGVFRRDTHK",
    "sp|P67890|PROT2": "MRGVFRRDTHKSEQ",
    "sp|Q11111|PROT3": "MXXUNIQUESEQXXX",
}


def _records(result: AnnotationResult) -> set[tuple[str, str, int, int]]:
    return {
        (annotation.peptide, annotation.protein_id, annotation.start, annotation.end)
        for annotation in result
    }


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_annotation_matches_mapping_and_streaming_input(backend: str) -> None:
    peptides = ["GVFRR", "DTHK", "UNIQUE"]
    mapped = annotate_peptides(peptides, PROTEINS, backend=backend)
    streamed = annotate_peptides_streaming(peptides, iter(PROTEINS.items()), backend=backend)
    assert _records(mapped) == _records(streamed)


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
