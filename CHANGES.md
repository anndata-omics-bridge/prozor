# Changes

- 2026-08-28: Relicensed from GPL-3.0-only to MIT. The Python package is an independent
  implementation of the R `prozor` algorithm rather than a translation of its source, and the same
  author holds copyright on both; the R package stays GPL-3. This removes the copyleft conflict with
  MIT-licensed APB, which imports `prozor.matching`. Added the `LICENSE` file and declared
  `license-files`, without which uv_build shipped wheels carrying no licence text.

- 2026-08-25: Split the algorithm core into independent `matching/` and `inference/`
  packages. Protein inference now consumes unique string edges directly; NumPy, SciPy, and
  `scipy-stubs` are removed. Consequential overlapping ties use an injected resolver, while
  identical evidence is grouped and disjoint components are processed independently. Reproducible
  benchmarks compare both Aho--Corasick backends and three greedy implementations; the selected
  integer-incidence implementation preserves all baseline result hashes and is faster than the
  former SciPy path on every measured graph shape. Rust matching is now installed and selected by
  default, while the pure-Python backend remains installed, explicitly selectable, and tested as the
  fallback.

- 2026-08-11: Add `make carpets`, running the sibling `carpet_scan` package over `src` and `tests`.
  Report-only and never blocking, matching this repo's Makefile idiom (`$(VENV_BIN)/carpet-scan`)
  rather than introducing a `manual` pre-commit stage this repo does not use. First run is **clean
  on every check**: 17 modules, no public function without an external caller, nothing at call depth
  >= 3, no dead code, no cross-module private access. `carpet_scan` needed its `requires-python`
  lowered from 3.13 to 3.12 to be installable here.
