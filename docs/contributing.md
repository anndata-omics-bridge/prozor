# Contributing

## Set up the environment

Prozor uses uv and a project-local `.venv`:

```bash
uv sync --group dev --group docs
.venv/bin/pre-commit install --hook-type pre-commit --hook-type pre-push
```

## Quality gates

```bash
make format       # apply Ruff formatting and safe lint fixes
make lint         # Ruff checks
make typecheck    # standard Pyright strict
make imports      # directed matching/inference boundary
make test         # pytest with branch coverage
make docs         # strict MkDocs build
make build        # wheel/sdist plus metadata validation
make check        # every merge-blocking gate
```

All Python commands use `.venv/bin`. Runtime dependencies, development tools,
and documentation tools must remain in their corresponding `pyproject.toml`
sections, and `uv.lock` must be updated with dependency changes.

## Code expectations

- Fully annotate source and tests; strict Pyright must remain at zero errors.
- Keep public APIs minimal and import them from their defining modules.
- Preserve nested and overlapping matches across both backends.
- Keep ordering and protein inference deterministic.
- Keep `matching/` and `inference/` independent; compose them in consumers.
- Add focused tests for behavioral changes and run `make check` before review.

See [`AGENTS.md`](https://github.com/anndata-omics-bridge/prozor/blob/main/AGENTS.md)
for the complete project conventions and dependency boundaries.
