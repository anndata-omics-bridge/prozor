# Design boundaries

Prozor is intentionally an algorithm package rather than a proteomics workflow.
This keeps peptide matching and protein inference reusable by APB,
`diann_runner`, and other consumers without coupling the core to their data
models.

## Dependency direction

```text
consumer application
  ├── FASTA parsing and identifier policy
  ├── sequence normalization and validation policy
  ├── decoy/contaminant classification
  ├── persistence and provenance
  └── prozor
       ├── Aho--Corasick matching
       ├── occurrence annotations
       ├── sparse peptide--protein topology
       └── deterministic greedy parsimony
```

Prozor must not import consumer packages. In particular, its runtime core does
not depend on pandas, AnnData, MuData, MuLink, FASTA parsers, or workflow
engines.

## Public module boundaries

| Module | Responsibility |
| --- | --- |
| `prozor.ahocorasick` | Backend selection and exact multi-pattern matching. |
| `prozor.annotate` | Peptide occurrences over mappings or streaming records. |
| `prozor.sparse_matrix` | Label-aware sparse peptide--protein topology. |
| `prozor.greedy` | Deterministic greedy-parsimony protein inference. |

The package `__init__.py` remains empty. Public objects are imported from their
defining modules so ownership stays explicit.

## Reproducibility guarantees

- backend requests and resolved implementations are both exposed;
- match coordinates are half-open and occurrence-level;
- sparse labels and inferred group members are sorted;
- duplicate occurrence sites do not change matrix topology;
- greedy tie-breaking is deterministic; and
- both matching backends are exercised by the same behavioral fixtures.

## Provenance

The Python implementation is derived from the behavior of Witold Wolski's R
[`prozor`](https://github.com/wolski/prozor) package. A fixed fixture generated
from the R 0.3.4 source guards equivalent protein-group membership; the Python
implementation canonicalizes group identifier ordering.
