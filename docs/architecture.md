# Architecture

Prozor is an algorithm package, not a proteomics workflow. Peptide matching and
protein inference remain reusable without coupling the core to a file format,
dataframe, analysis container, or consuming application.

## Dependency direction

```text
consumer application
  ├── FASTA parsing and identifier policy
  ├── sequence normalization and validation policy
  ├── decoy/contaminant classification
  ├── persistence and provenance
  └── prozor
       ├── matching/
       │    ├── automaton.py
       │    └── annotation.py
       └── inference/
            ├── ties.py
            └── greedy.py
```

`matching/` and `inference/` are independent algorithm families. Neither
imports the other; the consumer turns annotations into edges and composes them.
This boundary is enforced by import-linter. Both package markers are empty.

Production Prozor imports none of NumPy, SciPy, Pydantic, pandas, Polars,
AnnData, MuData, FASTA parsers, or workflow engines. FASTA parsing belongs in a
separate package because header schemas and validation are input policy, not
matching or inference algorithms.

## Public module boundaries

| Module | Responsibility |
| --- | --- |
| `prozor.matching.automaton` | Backend selection and exact multi-pattern matching. |
| `prozor.matching.annotation` | Peptide occurrences over mappings or streaming records. |
| `prozor.inference.ties` | Consequential tie value and default resolution operation. |
| `prozor.inference.greedy` | Edge-based deterministic greedy-parsimony inference. |

The root `__init__.py` remains empty. Public objects are imported from their
defining modules so ownership stays explicit.

## Reproducibility guarantees

- backend requests and resolved implementations are both exposed;
- `auto` selects Rust by default and retains a tested pure-Python fallback;
- match coordinates are half-open and occurrence-level;
- duplicate occurrences collapse to unique edges before inference;
- identical-evidence proteins are grouped rather than arbitrarily selected;
- disjoint candidates do not invoke the tie resolver;
- overlapping, non-identical ties receive complete unexplained evidence;
- inferred members and group output order are stable; and
- both matching backends are exercised by the same behavioral fixtures.

## Measured implementation choice

The benchmark suite retains readable name-set and lazy-priority-queue variants
outside production. The integer-indexed incidence implementation is the
production choice because it returns the same canonical groups while winning
all measured graph shapes. Independent overlap components are processed
together; this is both the correct tie boundary and the largest speedup.

Alternative tie scores and quantitative-profile clustering remain benchmark
experiments. Correlation can reveal incoherent peptide subsets, but a second
cluster without protein-unique evidence cannot identify which protein,
proteoform, or interference generated it.

## Executable enforcement

This document is not generated from Python. Its claims are kept executable in
the smallest appropriate places:

- `.importlinter` enforces independence between `matching/` and `inference/`
  and forbids data/storage framework imports;
- `tests/test_package.py` keeps package markers empty;
- tests under `tests/matching/` run the same behavioral contracts against both
  matching implementations; and
- tests under `tests/inference/` define grouping, tie, subsumption, and ordering
  behavior.

## Provenance

The Python implementation is derived from the behavior of Witold Wolski's R
[`prozor`](https://github.com/wolski/prozor) package. A fixed fixture generated
from the R 0.3.4 source guards equivalent protein-group membership; the Python
implementation canonicalizes group identifier ordering.
