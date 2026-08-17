# Changes

- 2026-08-11: Add `make carpets`, running the sibling `carpet_scan` package over `src` and `tests`.
  Report-only and never blocking, matching this repo's Makefile idiom (`$(VENV_BIN)/carpet-scan`)
  rather than introducing a `manual` pre-commit stage this repo does not use. First run is **clean
  on every check**: 17 modules, no public function without an external caller, nothing at call depth
  >= 3, no dead code, no cross-module private access. `carpet_scan` needed its `requires-python`
  lowered from 3.13 to 3.12 to be installable here.
