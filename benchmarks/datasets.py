"""Deterministic inputs shared by Prozor benchmarks."""

from __future__ import annotations

import random
from dataclasses import dataclass

_AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"

type PeptideProteinEdge = tuple[str, str]


@dataclass(frozen=True, slots=True)
class MatchingCase:
    """One Aho-Corasick benchmark input."""

    name: str
    keywords: tuple[str, ...]
    protein_records: tuple[tuple[str, str], ...]

    @property
    def residues(self) -> int:
        """Return the total number of protein-sequence residues."""
        return sum(len(sequence) for _protein, sequence in self.protein_records)


@dataclass(frozen=True, slots=True)
class InferenceCase:
    """One greedy-inference benchmark input."""

    name: str
    edges: tuple[PeptideProteinEdge, ...]

    @property
    def peptide_count(self) -> int:
        """Return the number of distinct peptides."""
        return len({peptide for peptide, _protein in self.edges})

    @property
    def protein_count(self) -> int:
        """Return the number of distinct proteins."""
        return len({protein for _peptide, protein in self.edges})


def matching_cases() -> tuple[MatchingCase, ...]:
    """Return deterministic sparse-match, dense-match, and overlap cases."""
    return (
        _random_matching_case(
            name="sparse_128x100k",
            seed=31,
            keyword_count=128,
            record_count=100,
            record_length=1_000,
            insertions_per_record=2,
        ),
        _random_matching_case(
            name="sparse_1024x500k",
            seed=37,
            keyword_count=1_024,
            record_count=250,
            record_length=2_000,
            insertions_per_record=3,
        ),
        MatchingCase(
            name="nested_overlap",
            keywords=("AA", "AAA", "AAAA", "AAAAA", "CAA", "AAC"),
            protein_records=tuple(
                (f"OVERLAP_{index:03d}", "CA" + "A" * 196 + "C") for index in range(40)
            ),
        ),
    )


def inference_cases() -> tuple[InferenceCase, ...]:
    """Return graph shapes that exercise different greedy costs."""
    return (
        _unique_inference_case(protein_count=600, peptides_per_protein=5),
        _component_inference_case(
            name="components_30x20",
            seed=41,
            component_count=30,
            proteins_per_component=20,
            peptides_per_component=80,
            proteins_per_peptide=3,
        ),
        _component_inference_case(
            name="dense_component",
            seed=43,
            component_count=1,
            proteins_per_component=300,
            peptides_per_component=2_000,
            proteins_per_peptide=8,
        ),
        _overlapping_tie_case(component_count=500),
    )


def _random_matching_case(
    *,
    name: str,
    seed: int,
    keyword_count: int,
    record_count: int,
    record_length: int,
    insertions_per_record: int,
) -> MatchingCase:
    generator = random.Random(seed)
    keywords = _unique_random_peptides(generator, keyword_count)
    records: list[tuple[str, str]] = []
    for record_index in range(record_count):
        sequence = [generator.choice(_AMINO_ACIDS) for _ in range(record_length)]
        for insertion_index in range(insertions_per_record):
            keyword = keywords[
                (record_index * insertions_per_record + insertion_index) % keyword_count
            ]
            available = record_length - len(keyword)
            start = generator.randrange(available + 1)
            sequence[start : start + len(keyword)] = keyword
        records.append((f"PROTEIN_{record_index:05d}", "".join(sequence)))
    return MatchingCase(name=name, keywords=keywords, protein_records=tuple(records))


def _unique_random_peptides(generator: random.Random, count: int) -> tuple[str, ...]:
    peptides: dict[str, None] = {}
    while len(peptides) < count:
        length = generator.randrange(7, 19)
        peptide = "".join(generator.choice(_AMINO_ACIDS) for _ in range(length))
        peptides[peptide] = None
    return tuple(peptides)


def _unique_inference_case(
    *,
    protein_count: int,
    peptides_per_protein: int,
) -> InferenceCase:
    edges = tuple(
        (f"PEPTIDE_{protein:05d}_{peptide:02d}", f"PROTEIN_{protein:05d}")
        for protein in range(protein_count)
        for peptide in range(peptides_per_protein)
    )
    return InferenceCase(name="mostly_unique", edges=edges)


def _component_inference_case(
    *,
    name: str,
    seed: int,
    component_count: int,
    proteins_per_component: int,
    peptides_per_component: int,
    proteins_per_peptide: int,
) -> InferenceCase:
    generator = random.Random(seed)
    edges: list[PeptideProteinEdge] = []
    for component in range(component_count):
        proteins = [
            f"PROTEIN_{component:04d}_{protein:04d}" for protein in range(proteins_per_component)
        ]
        for protein in proteins:
            edges.append((f"ANCHOR_{protein}", protein))
        for peptide in range(peptides_per_component):
            peptide_name = f"PEPTIDE_{component:04d}_{peptide:05d}"
            for protein in generator.sample(proteins, proteins_per_peptide):
                edges.append((peptide_name, protein))
    return InferenceCase(name=name, edges=tuple(edges))


def _overlapping_tie_case(*, component_count: int) -> InferenceCase:
    edges: list[PeptideProteinEdge] = []
    for component in range(component_count):
        protein_a = f"TIE_A_{component:05d}"
        protein_b = f"TIE_B_{component:05d}"
        edges.extend(
            (
                (f"LEFT_{component:05d}", protein_a),
                (f"SHARED_{component:05d}", protein_a),
                (f"SHARED_{component:05d}", protein_b),
                (f"RIGHT_{component:05d}", protein_b),
            )
        )
    return InferenceCase(name="overlapping_ties", edges=tuple(edges))
