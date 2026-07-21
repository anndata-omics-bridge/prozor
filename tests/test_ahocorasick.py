from __future__ import annotations

import pytest

import prozor.ahocorasick as aho
from prozor.ahocorasick import Match, create_automaton, get_available_backends


def _match_tuples(keywords: list[str], text: str, backend: str) -> set[tuple[str, int, int]]:
    return {
        (match.keyword, match.start, match.end)
        for match in create_automaton(keywords, backend=backend).find_all(text)
    }


def test_single_and_absent_matches() -> None:
    assert list(create_automaton(["PEPTIDE"]).find_all("MYPEPTIDE")) == [
        Match(keyword="PEPTIDE", start=2, end=9)
    ]
    assert list(create_automaton(["ABSENT"]).find_all("MYPEPTIDE")) == []


def test_nested_overlapping_and_repeated_matches() -> None:
    matches = _match_tuples(
        ["SAMPLER", "SAMPLERPEPTIDER", "AA"],
        "SAMPLERPEPTIDERXAAA",
        "ahocorapy",
    )
    assert ("SAMPLER", 0, 7) in matches
    assert ("SAMPLERPEPTIDER", 0, 15) in matches
    assert {match for match in matches if match[0] == "AA"} == {
        ("AA", 16, 18),
        ("AA", 17, 19),
    }


@pytest.mark.parametrize("backend", get_available_backends())
def test_all_backends_have_identical_exact_matches(backend: str) -> None:
    keywords = ["SAMPLER", "SAMPLERPEPTIDER", "GVFRR", "AA"]
    text = "MKSAMPLERPEPTIDERKGVFRRXAAAX"
    assert _match_tuples(keywords, text, backend) == _match_tuples(
        keywords,
        text,
        "ahocorapy",
    )


def test_auto_records_requested_and_resolved_backend() -> None:
    automaton = create_automaton(["PEP"], backend="auto")
    assert automaton.requested_backend == "auto"
    assert automaton.resolved_backend in get_available_backends()


def test_auto_falls_back_when_rust_is_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(aho, "_rust_available", lambda: False)
    automaton = create_automaton(["PEP"], backend="auto")
    assert automaton.requested_backend == "auto"
    assert automaton.resolved_backend == "ahocorapy"


def test_explicit_missing_rust_backend_fails(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(aho, "_rust_available", lambda: False)
    with pytest.raises(ImportError, match=r"prozor\[fast\]"):
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
