# Prozor

[![Quality](https://github.com/anndata-omics-bridge/prozor/actions/workflows/quality.yml/badge.svg)](https://github.com/anndata-omics-bridge/prozor/actions/workflows/quality.yml)
[![Documentation](https://github.com/anndata-omics-bridge/prozor/actions/workflows/docs.yml/badge.svg)](https://anndata-omics-bridge.github.io/prozor/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.12-3776AB.svg)](https://www.python.org/)
[![License: GPL-3.0-only](https://img.shields.io/badge/License-GPL--3.0--only-blue.svg)](https://spdx.org/licenses/GPL-3.0-only.html)

Typed peptide-to-protein matching and deterministic greedy-parsimony protein
inference.

Prozor offers backend-neutral Aho--Corasick matching over mappings or streaming
protein records, occurrence-level annotations, sparse peptide--protein
topology, and deterministic protein inference. The core deliberately does not
parse FASTA files or depend on pandas, AnnData, MuData, MuLink, workflow engines,
or consumer CLIs.

**[Documentation](https://anndata-omics-bridge.github.io/prozor/)** ·
**[Getting started](https://anndata-omics-bridge.github.io/prozor/getting-started/)** ·
**[API reference](https://anndata-omics-bridge.github.io/prozor/api/)**

## Installation

Until the first package-index release, install from GitHub:

```bash
python -m pip install "prozor @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

Install the optional Rust matcher with:

```bash
python -m pip install "prozor[fast] @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

The portable `ahocorapy` backend is always available. With `backend="auto"`,
Prozor prefers `ahocorasick_rs` when installed and records both the requested
and concrete backend in the result.

## Quick start

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

for match in matches:
    print(match.peptide, match.protein_id, match.start, match.end)

protein_groups = greedy_parsimony(matches.to_sparse_matrix())
print(protein_groups.to_dict())
```

Matches include nested, overlapping, and repeated sites using half-open
coordinates. Prozor matches the exact strings supplied by the consumer; FASTA
header interpretation, normalization, decoy classification, and persistence
remain application policy.

## Development

```bash
uv sync --group dev --group docs
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
make check
```

`make check` runs Ruff, strict Pyright, dependency validation, branch-coverage
tests, a strict documentation build, and wheel/sdist validation. All Python
commands run from the synchronized project `.venv`.

## Provenance and license

This implementation is derived from the algorithmic behavior of Witold
Wolski's R [`prozor`](https://github.com/wolski/prozor) package and consolidates
the matching and inference core for Python consumers such as APB and
`diann_runner`.

The Python distribution is licensed under GPL-3.0-only, matching the declared
license of the R reference package.
