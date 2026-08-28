# Prozor

[![Quality](https://github.com/anndata-omics-bridge/prozor/actions/workflows/quality.yml/badge.svg)](https://github.com/anndata-omics-bridge/prozor/actions/workflows/quality.yml)
[![Documentation](https://github.com/anndata-omics-bridge/prozor/actions/workflows/docs.yml/badge.svg)](https://anndata-omics-bridge.github.io/prozor/)
[![Python](https://img.shields.io/badge/Python-%E2%89%A53.12-3776AB.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://spdx.org/licenses/MIT.html)

Typed peptide-to-protein matching and deterministic greedy-parsimony protein
inference.

Prozor offers backend-neutral Aho--Corasick matching over mappings or streaming
protein records, occurrence-level annotations, unique peptide--protein edges,
and deterministic protein inference. The core deliberately does not
parse FASTA files or depend on pandas, AnnData, MuData, MuLink, workflow engines,
or consumer CLIs.

**[Documentation](https://anndata-omics-bridge.github.io/prozor/)** ·
**[Getting started](https://anndata-omics-bridge.github.io/prozor/getting-started/)** ·
**[API reference](https://anndata-omics-bridge.github.io/prozor/api/)**

Documentation: <https://anndata-omics-bridge.github.io/prozor/>

## Installation

Until the first package-index release, install from GitHub:

```bash
python -m pip install "prozor @ git+https://github.com/anndata-omics-bridge/prozor.git"
```

Both matching implementations are installed. The public matching operations
default to `backend="auto"`, which selects `ahocorasick_rs`. The portable
`ahocorapy` implementation remains directly selectable and is the automatic
runtime fallback if Rust cannot be imported. Results record both the requested
and concrete backend.

## Quick start

```python
from prozor.inference.greedy import greedy_parsimony
from prozor.matching.annotation import annotate_peptides_streaming

matches = annotate_peptides_streaming(
    ["PEPTIDE", "SEQUENCE"],
    [
        ("P1", "MYPEPTIDESEQUENCE"),
        ("P2", "XXSEQUENCEXX"),
    ],
)

for match in matches:
    print(match.peptide, match.protein_id, match.start, match.end)

edges = {(match.peptide, match.protein_id) for match in matches}
protein_groups = greedy_parsimony(edges)
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

The Python distribution is an independent implementation of that algorithm
rather than a translation of the R source, and it is licensed under MIT. The R
reference package remains GPL-3, and its license does not extend here: the
copyright holder of both packages is the same author, and algorithms are not
themselves subject to copyright.

MIT keeps the whole anndata-omics-bridge spine under one permissive license, so
that consumers such as APB can be redistributed without inheriting copyleft
obligations from this package.
