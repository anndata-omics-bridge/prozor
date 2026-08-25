from __future__ import annotations

from collections.abc import Iterator

import pytest

import prozor.matching.automaton as aho
from prozor.matching.automaton import Match, create_automaton, get_available_backends


def _match_tuples(keywords: list[str], text: str, backend: str) -> set[tuple[str, int, int]]:
    return {
        (match.keyword, match.start, match.end)
        for match in create_automaton(keywords, backend=backend).find_all(text)
    }


def test_development_environment_exercises_both_backends() -> None:
    assert get_available_backends() == ["ahocorapy", "ahocorasick_rs"]


@pytest.mark.parametrize("backend", ["ahocorapy", "ahocorasick_rs"])
def test_nested_overlapping_and_repeated_match_contract(backend: str) -> None:
    matches = _match_tuples(
        ["SAMPLER", "SAMPLERPEPTIDER", "AA", "ABSENT"],
        "SAMPLERPEPTIDERXAAA",
        backend,
    )
    assert matches == {
        ("SAMPLER", 0, 7),
        ("SAMPLERPEPTIDER", 0, 15),
        ("AA", 16, 18),
        ("AA", 17, 19),
    }


def test_cross_backend_results_are_identical() -> None:
    keywords = ["SAMPLER", "SAMPLERPEPTIDER", "GVFRR", "AA"]
    text = "MKSAMPLERPEPTIDERKGVFRRXAAAX"
    assert _match_tuples(keywords, text, "ahocorasick_rs") == _match_tuples(
        keywords,
        text,
        "ahocorapy",
    )


def test_single_match_uses_half_open_coordinates() -> None:
    assert list(create_automaton(["PEPTIDE"]).find_all("MYPEPTIDE")) == [
        Match(keyword="PEPTIDE", start=2, end=9)
    ]


def test_auto_records_requested_and_resolved_backend() -> None:
    automaton = create_automaton(["PEP"], backend="auto")
    assert automaton.requested_backend == "auto"
    assert automaton.resolved_backend == "ahocorasick_rs"


def test_auto_falls_back_when_rust_is_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(aho, "_rust_available", lambda: False)
    automaton = create_automaton(["PEP"], backend="auto")
    assert automaton.requested_backend == "auto"
    assert automaton.resolved_backend == "ahocorapy"


def test_explicit_missing_rust_backend_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(aho, "_rust_available", lambda: False)
    with pytest.raises(ImportError, match="working ahocorasick-rs installation"):
        create_automaton(["PEP"], backend="ahocorasick_rs")


def test_invalid_backend_fails() -> None:
    with pytest.raises(ValueError, match="backend must be one of"):
        create_automaton(["PEP"], backend="unsupported")


def test_case_insensitive_matching_preserves_original_keyword() -> None:
    match = next(create_automaton(["PeP"], case_sensitive=False).find_all("xxpep"))
    assert match == Match(keyword="PeP", start=2, end=5)


def test_match_is_immutable() -> None:
    match = Match(keyword="PEP", start=0, end=3)
    attribute = "keyword"
    with pytest.raises(AttributeError):
        setattr(match, attribute, "OTHER")


def test_automaton_iterator_is_lazy() -> None:
    matches = create_automaton(["PEP"]).find_all("PEP")
    assert isinstance(matches, Iterator)
