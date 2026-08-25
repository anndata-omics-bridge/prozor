# Prozor algorithm benchmark — 2026-08-25

The before and after JSON files were produced in the same Prozor environment
with three measured repetitions after a warm-up. Every after implementation
was rejected before timing unless its canonical result hash matched the
readable reference. The production hashes also match the former SciPy path for
all four graph shapes.

The human-readable interpretation lives here with the rest of the maintained
documentation. Reproducibility code remains in `benchmarks/`, and its
machine-readable outputs remain in `benchmarks/results/`; they are evidence,
not production package modules.

## Matching

The source move did not change matching behavior or materially change timings.
On the after run, `ahocorasick_rs` scanned the two sparse fixtures 10.1x and
13.1x faster than `ahocorapy`. On the deliberately dense nested-overlap fixture
it was 1.18x faster because materializing 31,200 Python match records dominates
the automaton search.

## Greedy inference

Median total seconds, including construction of the inference representation:

| Graph shape | Former SciPy CSR | Integer incidence | Speedup |
| --- | ---: | ---: | ---: |
| Mostly unique, 3,000 edges | 0.083891 | 0.005081 | 16.51x |
| 30 components, 7,800 edges | 0.040341 | 0.009132 | 4.42x |
| One dense component, 16,300 edges | 0.040615 | 0.024529 | 1.66x |
| 500 overlapping tie components, 2,000 edges | 0.162191 | 0.007770 | 20.87x |

The integer-indexed implementation won every measured graph shape and is the
production implementation. The direct-name-set oracle and lazy priority queue
remain benchmark-only. Processing independent overlap components in one pass
was decisive for the mostly-unique and many-tie fixtures; it also ensures an
injected resolver sees one consequential connected tie rather than unrelated
candidates.

`tracemalloc` reports Python allocations only and cannot see SciPy's native
buffers, so its before/after peak numbers are retained in JSON but are not used
to claim a memory improvement.

## Tie-score and quantitative-profile exploration

`benchmarks/results/tie_scoring_after.json` compares the current resolver with inverse peptide
degeneracy, sequence coverage, and external-score examples. These results are
exploratory and do not change production policy.

The profile simulations support a narrower conclusion: correlation can expose
two coherent peptide clusters, but cannot by itself identify their biological
source. When the second cluster contains a peptide unique to protein B, it is
evidence for B. When it contains only peptides shared with B, the same pattern
is also compatible with proteoforms, interference, or other heterogeneity and
must remain unresolved.

## Reproduction

```bash
PROZOR_BENCHMARK_LABEL=after .venv/bin/python -m benchmarks.benchmark_matching
PROZOR_BENCHMARK_LABEL=after .venv/bin/python -m benchmarks.benchmark_greedy
.venv/bin/python -m benchmarks.explore_tie_scoring
```
