"""Benchmark Prozor's Aho-Corasick implementations over identical inputs."""

from __future__ import annotations

import gc
import importlib.metadata
import json
import os
import platform
import statistics
import subprocess
import time
import tracemalloc
from collections.abc import Callable, Iterable
from pathlib import Path
from typing import TypedDict

from benchmarks.datasets import MatchingCase, matching_cases
from prozor.matching.automaton import AhoCorasickBase, create_automaton, get_available_backends

_REPETITIONS = 3

type MatchRecord = tuple[str, str, int, int]


class Timing(TypedDict):
    """Repeated timing summary in seconds."""

    median_seconds: float
    seconds: list[float]


def main() -> None:
    """Run matching benchmarks and write one machine-readable result."""
    cases = matching_cases()
    backends = get_available_backends()
    if set(backends) != {"ahocorapy", "ahocorasick_rs"}:
        raise RuntimeError(f"benchmark requires both backends, found {backends!r}")

    results: list[dict[str, object]] = []
    for case in cases:
        reference = _canonical_matches(case, "ahocorapy")
        for backend in backends:
            results.append(_benchmark_backend(case, backend, reference))
    _write_result("matching", results)


def _benchmark_backend(
    case: MatchingCase,
    backend: str,
    reference: frozenset[MatchRecord],
) -> dict[str, object]:
    observed = _canonical_matches(case, backend)
    if observed != reference:
        raise AssertionError(f"{backend} disagrees with ahocorapy for {case.name}")
    automaton = create_automaton(case.keywords, backend=backend)
    build = _measure(lambda: create_automaton(case.keywords, backend=backend))
    scan = _measure(lambda: _scan(automaton, case.protein_records))
    python_peak_bytes = _python_peak_bytes(
        lambda: _scan(
            create_automaton(case.keywords, backend=backend),
            case.protein_records,
        )
    )
    return {
        "case": case.name,
        "backend": backend,
        "keyword_count": len(case.keywords),
        "protein_count": len(case.protein_records),
        "residues": case.residues,
        "match_count": len(observed),
        "build": build,
        "scan": scan,
        "scan_residues_per_second": case.residues / scan["median_seconds"],
        "python_peak_bytes": python_peak_bytes,
    }


def _canonical_matches(case: MatchingCase, backend: str) -> frozenset[MatchRecord]:
    automaton = create_automaton(case.keywords, backend=backend)
    return frozenset(
        (protein, match.keyword, match.start, match.end)
        for protein, sequence in case.protein_records
        for match in automaton.find_all(sequence)
    )


def _scan(
    automaton: AhoCorasickBase,
    protein_records: Iterable[tuple[str, str]],
) -> int:
    return sum(
        1 for _protein, sequence in protein_records for _match in automaton.find_all(sequence)
    )


def _measure(operation: Callable[[], object]) -> Timing:
    measurements: list[float] = []
    operation()
    for _index in range(_REPETITIONS):
        gc.collect()
        started = time.perf_counter()
        operation()
        measurements.append(time.perf_counter() - started)
    return {
        "median_seconds": statistics.median(measurements),
        "seconds": measurements,
    }


def _python_peak_bytes(operation: Callable[[], object]) -> int:
    gc.collect()
    tracemalloc.start()
    operation()
    _current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    return peak


def _write_result(kind: str, results: list[dict[str, object]]) -> None:
    label = os.environ.get("PROZOR_BENCHMARK_LABEL", "manual")
    output = Path(__file__).parent / "results" / f"{kind}_{label}.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "label": label,
        "commit": _commit(),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "prozor": importlib.metadata.version("prozor"),
        "ahocorapy": importlib.metadata.version("ahocorapy"),
        "ahocorasick_rs": importlib.metadata.version("ahocorasick-rs"),
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
