"""Peptide-to-protein annotation over mappings or streaming protein records."""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass
from typing import TYPE_CHECKING

from prozor.ahocorasick import BackendName, BackendRequest, create_automaton, resolve_backend

if TYPE_CHECKING:
    from prozor.sparse_matrix import PeptideProteinMatrix


@dataclass(frozen=True, slots=True)
class PeptideAnnotation:
    """One peptide occurrence in a protein sequence."""

    peptide: str
    protein_id: str
    start: int
    end: int

    @property
    def length(self) -> int:
        """Return the matched peptide length."""
        return len(self.peptide)


@dataclass(slots=True)
class AnnotationResult:
    """Collection of peptide matches and matching-backend provenance."""

    annotations: list[PeptideAnnotation]
    requested_backend: BackendRequest = "ahocorapy"
    resolved_backend: BackendName = "ahocorapy"

    def __len__(self) -> int:
        return len(self.annotations)

    def __iter__(self) -> Iterator[PeptideAnnotation]:
        return iter(self.annotations)

    @property
    def peptides(self) -> set[str]:
        """Return the distinct matched peptide sequences."""
        return {annotation.peptide for annotation in self.annotations}

    @property
    def proteins(self) -> set[str]:
        """Return the distinct matched protein identifiers."""
        return {annotation.protein_id for annotation in self.annotations}

    def filter_tryptic(
        self,
        proteins: Mapping[str, str],
        prefix_residues: str = "RK",
        allow_n_term: bool = True,
        allow_after_init_met: bool = True,
    ) -> AnnotationResult:
        """Return matches with a supported tryptic N-terminal context."""
        filtered: list[PeptideAnnotation] = []
        for annotation in self.annotations:
            sequence = proteins.get(annotation.protein_id, "")
            if not sequence:
                continue
            valid_prefix = (
                (annotation.start == 0 and allow_n_term)
                or (annotation.start == 1 and allow_after_init_met and sequence.startswith("M"))
                or (annotation.start > 0 and sequence[annotation.start - 1] in prefix_residues)
            )
            if valid_prefix:
                filtered.append(annotation)
        return AnnotationResult(
            annotations=filtered,
            requested_backend=self.requested_backend,
            resolved_backend=self.resolved_backend,
        )

    def to_sparse_matrix(self, weighting: str = "binary") -> PeptideProteinMatrix:
        """Build the sparse peptide--protein topology for these matches."""
        from prozor.sparse_matrix import PeptideProteinMatrix

        edges = ((annotation.peptide, annotation.protein_id) for annotation in self.annotations)
        return PeptideProteinMatrix.from_edges(edges, weighting=weighting)


def annotate_peptides(
    peptides: Iterable[str],
    proteins: Mapping[str, str],
    backend: str = "auto",
    filter_tryptic: bool = False,
) -> AnnotationResult:
    """Annotate peptides against an in-memory protein mapping."""
    result = annotate_peptides_streaming(peptides, proteins.items(), backend=backend)
    return result.filter_tryptic(proteins) if filter_tryptic else result


def annotate_peptides_streaming(
    peptides: Iterable[str],
    protein_records: Iterable[tuple[str, str]],
    backend: str = "auto",
) -> AnnotationResult:
    """Annotate peptides against one-pass ``(protein_id, sequence)`` records."""
    peptide_list = list(dict.fromkeys(peptides))
    if not peptide_list:
        return AnnotationResult(
            annotations=[],
            requested_backend=_requested_backend(backend),
            resolved_backend=resolve_backend(backend),
        )

    automaton = create_automaton(peptide_list, backend=backend)
    annotations: list[PeptideAnnotation] = []
    for protein_id, sequence in protein_records:
        annotations.extend(
            PeptideAnnotation(
                peptide=match.keyword,
                protein_id=protein_id,
                start=match.start,
                end=match.end,
            )
            for match in automaton.find_all(sequence)
        )
    return AnnotationResult(
        annotations=annotations,
        requested_backend=automaton.requested_backend,
        resolved_backend=automaton.resolved_backend,
    )


def _requested_backend(backend: str) -> BackendRequest:
    # resolve_backend provides the public runtime validation. The three branches
    # below narrow the already-validated value for the type checker.
    resolve_backend(backend)
    if backend == "ahocorapy":
        return "ahocorapy"
    if backend == "ahocorasick_rs":
        return "ahocorasick_rs"
    return "auto"
