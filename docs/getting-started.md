# Getting started

## Requirements

Prozor supports Python 3.12 and newer.

## Installation

Until the first package-index release, install directly from GitHub:

```bash
python -m pip install "prozor @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

The portable `ahocorapy` implementation is always installed. To also install
the accelerated Rust backend:

```bash
python -m pip install "prozor[fast] @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

Once released on a package index, the corresponding commands are
`pip install prozor` and `pip install "prozor[fast]"`.

## Match peptides against proteins

Use a mapping when all proteins are already in memory:

```python
from prozor.annotate import annotate_peptides

proteins = {
    "P1": "MYPEPTIDESEQUENCE",
    "P2": "XXSEQUENCEXX",
}
result = annotate_peptides(["PEPTIDE", "SEQUENCE"], proteins)

for match in result:
    print(match.peptide, match.protein_id, match.start, match.end)
```

Use a one-pass iterable when a consumer already streams FASTA or database
records:

```python
from collections.abc import Iterator

from prozor.annotate import annotate_peptides_streaming


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

Convert occurrence-level matches to a deduplicated sparse topology and run
greedy parsimony:

```python
from prozor.greedy import greedy_parsimony

matrix = result.to_sparse_matrix()
inference = greedy_parsimony(matrix)

for group in inference:
    print(group.protein_id, group.peptides)
```

See [Protein inference](protein-inference.md) for weighting, grouping, and
subsumption semantics.
