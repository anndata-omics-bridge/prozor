# Getting started

## Requirements

Prozor supports Python 3.12 and newer.

## Installation

Until the first package-index release, install directly from GitHub:

```bash
python -m pip install "prozor @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

Both the accelerated `ahocorasick_rs` implementation and the portable
`ahocorapy` fallback are installed. Public matching operations default to
`backend="auto"`: Rust is selected when importable, otherwise matching falls
back to pure Python. Once released on a package index, install with
`pip install prozor`.

## Match peptides against proteins

Use a mapping when all proteins are already in memory:

```python
from prozor.matching.annotation import annotate_peptides

proteins = {
    "P1": "MYPEPTIDESEQUENCE",
    "P2": "XXSEQUENCEXX",
}
result = annotate_peptides(["PEPTIDE", "SEQUENCE"], proteins)

for match in result:
    print(match.peptide, match.protein_id, match.start, match.end)
```

The call above uses Rust by default. Pass `backend="ahocorapy"` to select the
portable implementation explicitly.

Use a one-pass iterable when a consumer already streams FASTA or database
records:

```python
from collections.abc import Iterator

from prozor.matching.annotation import annotate_peptides_streaming


def protein_records() -> Iterator[tuple[str, str]]:
    yield "P1", "MYPEPTIDESEQUENCE"
    yield "P2", "XXSEQUENCEXX"


result = annotate_peptides_streaming(
    ["PEPTIDE", "SEQUENCE"],
    protein_records(),
)
```

!!! important

    Prozor matches the strings it receives. It does not uppercase sequences,
    remove modifications, apply I/L equivalence, classify decoys, or parse
    FASTA headers. Those policies belong to the consuming application.

## Infer protein groups

Convert occurrence-level matches to unique peptide--protein edges and run
greedy parsimony. Repeated sites collapse naturally in the set:

```python
from prozor.inference.greedy import greedy_parsimony

edges = {(match.peptide, match.protein_id) for match in result}
inference = greedy_parsimony(edges)

for group in inference:
    print(group.protein_id, group.peptides)
```

See [Protein inference](protein-inference.md) for grouping, tie resolution, and
subsumption semantics.
