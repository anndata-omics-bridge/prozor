# Prozor

Prozor is a typed Python library for matching peptide sequences to proteins and
inferring deterministic parsimonious protein groups.

It provides two connected building blocks:

- backend-neutral Aho--Corasick matching over in-memory or streaming protein
  records; and
- sparse peptide--protein topology with deterministic greedy-parsimony protein
  inference.

[Get started](getting-started.md){ .md-button .md-button--primary }
[API reference](api.md){ .md-button }

## Why Prozor?

- **Complete matching:** nested, repeated, and overlapping peptide occurrences
  are retained with half-open coordinates.
- **Streaming input:** protein records can be consumed from any one-pass
  iterable; Prozor does not prescribe a FASTA parser.
- **Auditable backends:** results record both the requested backend and the
  concrete implementation selected at runtime.
- **Deterministic inference:** equivalent protein groups and tie-breaking use a
  stable ordering.
- **Small boundary:** the core does not depend on pandas, AnnData, MuData,
  MuLink, workflow engines, or consumer CLIs.

## Minimal workflow

```python
from prozor.annotate import annotate_peptides_streaming
from prozor.greedy import greedy_parsimony

matches = annotate_peptides_streaming(
    ["PEPTIDE", "SEQUENCE"],
    [
        ("P1", "MYPEPTIDESEQUENCE"),
        ("P2", "XXSEQUENCEXX"),
    ],
)

protein_groups = greedy_parsimony(matches.to_sparse_matrix())
print(protein_groups.to_dict())
```

## Project status

The Python package is derived from the algorithmic behavior of the R
[`prozor`](https://github.com/wolski/prozor) reference package. The matching and
inference core is covered by pure-Python/Rust backend-equivalence tests, strict
Pyright, Ruff, branch coverage, dependency validation, and distribution builds.
