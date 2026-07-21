# Prozor

Typed peptide-to-protein matching and deterministic greedy-parsimony protein
inference. This Python implementation is derived from the algorithmic behavior
of Witold Wolski's R `prozor` package and consolidates the implementations that
were previously copied into APB and `diann_runner`.

The core accepts protein records; it deliberately does not parse FASTA files or
depend on pandas, AnnData, MuData, MuLink, workflow engines, or consumer CLIs.

## Usage

```python
from prozor.annotate import annotate_peptides_streaming
from prozor.greedy import greedy_parsimony

matches = annotate_peptides_streaming(
    ["PEPTIDE", "SEQUENCE"],
    [("P1", "MYPEPTIDESEQUENCE")],
)
protein_groups = greedy_parsimony(matches.to_sparse_matrix())
```

`backend="auto"` prefers the optional Rust implementation while recording both
the requested and concrete backend. Install `prozor[fast]` to require the Rust
backend; `ahocorapy` remains the portable implementation.

## Development

```bash
uv sync --group dev
make check
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

All Python commands run from the synchronized project `.venv`.

## Provenance and license

The R reference implementation is maintained at
<https://github.com/wolski/prozor>. This Python distribution is licensed under
GPL-3.0-only, matching the declared license of that reference package.
