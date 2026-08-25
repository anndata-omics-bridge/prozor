"""Benchmark equivalent pure-Python greedy-inference implementations."""

from __future__ import annotations

import gc
import hashlib
import importlib.metadata
import json
import os
import platform
import random
import statistics
import subprocess
import time
import tracemalloc
from collections.abc import Callable, Sequence
from pathlib import Path
from typing import TypedDict

from benchmarks.datasets import InferenceCase, inference_cases
from benchmarks.greedy_implementations import (
    GreedyImplementation,
    direct_name_sets,
    lazy_priority_queue,
    production_incremental,
)
from prozor.inference.greedy import GreedyResult
from prozor.inference.ties import TieCandidate, resolve_current_tie

_REPETITIONS = 3
_IMPLEMENTATIONS: tuple[tuple[str, GreedyImplementation], ...] = (
    ("direct_name_sets", direct_name_sets),
    ("integer_incremental", production_incremental),
    ("lazy_priority_queue", lazy_priority_queue),
)


class Timing(TypedDict):
    """Repeated timing summary in seconds."""

    median_seconds: float
    seconds: list[float]


class RunResult(TypedDict):
    """Inference result plus the count of consequential ties."""

    result: GreedyResult
    genuine_ties: int


def main() -> None:
    """Run inference benchmarks and write one machine-readable result."""
    results = [
        _benchmark_implementation(case, name, implementation)
        for case in inference_cases()
        for name, implementation in _IMPLEMENTATIONS
    ]
    _write_result(results)


def _benchmark_implementation(
    case: InferenceCase,
    name: str,
    implementation: GreedyImplementation,
) -> dict[str, object]:
    expected = _run(direct_name_sets, case)["result"]
    observed = _run(implementation, case)
    result = observed["result"]
    if _result_digest(result) != _result_digest(expected):
        raise AssertionError(f"{name} disagrees with direct_name_sets for {case.name}")
    shuffled_edges = list(case.edges)
    random.Random(71).shuffle(shuffled_edges)
    shuffled_case = InferenceCase(name=case.name, edges=tuple(shuffled_edges))
    if _result_digest(_run(implementation, shuffled_case)["result"]) != _result_digest(result):
        raise AssertionError(f"{name} depends on edge order for {case.name}")
    total = _measure(lambda: _run(implementation, case))
    return {
        "case": case.name,
        "implementation": name,
        "edge_count": len(case.edges),
        "peptide_count": case.peptide_count,
        "protein_count": case.protein_count,
        "group_count": result.n_groups,
        "round_count": result.n_groups,
        "genuine_tie_count": observed["genuine_ties"],
        "result_sha256": _result_digest(result),
        "total": total,
        "python_peak_bytes": _python_peak_bytes(lambda: _run(implementation, case)),
    }


def _run(implementation: GreedyImplementation, case: InferenceCase) -> RunResult:
    tie_count = 0

    def counted_resolver(candidates: Sequence[TieCandidate]) -> TieCandidate:
        nonlocal tie_count
        tie_count += 1
        return resolve_current_tie(candidates)

    return {
        "result": implementation(case.edges, counted_resolver),
        "genuine_ties": tie_count,
    }


def _result_digest(result: GreedyResult) -> str:
    groups = [{"proteins": group.proteins, "peptides": group.peptides} for group in result.groups]
    encoded = json.dumps(groups, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def _measure(operation: Callable[[], object]) -> Timing:
    measurements: list[float] = []
    operation()
    for _index in range(_REPETITIONS):
        gc.collect()
        started = time.perf_counter()
        operation()
        measurements.append(time.perf_counter() - started)
    return {"median_seconds": statistics.median(measurements), "seconds": measurements}


def _python_peak_bytes(operation: Callable[[], object]) -> int:
    gc.collect()
    tracemalloc.start()
    operation()
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def _write_result(results: list[dict[str, object]]) -> None:
    label = os.environ.get("PROZOR_BENCHMARK_LABEL", "manual")
    output = Path(__file__).parent / "results" / f"greedy_{label}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "label": label,
        "commit": _commit(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "prozor": importlib.metadata.version("prozor"),
        "repetitions": _REPETITIONS,
        "results": results,
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def _commit() -> str:
    completed = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


if __name__ == "__main__":
    main()
