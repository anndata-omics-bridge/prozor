"""Backend-neutral Aho--Corasick peptide matching."""

from __future__ import annotations

from abc import ABC, abstractmethod
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from importlib.util import find_spec
from typing import Literal, cast, override

type BackendName = Literal["ahocorapy", "ahocorasick_rs"]
type BackendRequest = Literal["auto", "ahocorapy", "ahocorasick_rs"]

_VALID_BACKENDS = frozenset({"auto", "ahocorapy", "ahocorasick_rs"})


@dataclass(frozen=True, slots=True)
class Match:
    """One half-open keyword match in a searched sequence."""

    keyword: str
    start: int
    end: int


class AhoCorasickBase(ABC):
    """Common interface implemented by every matching backend."""

    requested_backend: BackendRequest
    resolved_backend: BackendName

    @abstractmethod
    def find_all(self, text: str) -> Iterator[Match]:
        """Yield every match, including nested and overlapping matches."""


class AhoCorasickPure(AhoCorasickBase):
    """Portable matcher implemented with :mod:`ahocorapy`."""

    def __init__(
        self,
        keywords: Iterable[str],
        *,
        case_sensitive: bool = True,
        requested_backend: BackendRequest = "ahocorapy",
    ) -> None:
        from ahocorapy.keywordtree import KeywordTree

        self.requested_backend = requested_backend
        self.resolved_backend: BackendName = "ahocorapy"
        self._tree = KeywordTree(case_insensitive=not case_sensitive)
        for keyword in _unique_keywords(keywords):
            self._tree.add(keyword)
        self._tree.finalize()

    @override
    def find_all(self, text: str) -> Iterator[Match]:
        """Yield all pure-Python matches in ``text``."""
        for keyword, start in self._tree.search_all(text):
            yield Match(keyword=keyword, start=start, end=start + len(keyword))


class AhoCorasickRust(AhoCorasickBase):
    """Accelerated matcher implemented with :mod:`ahocorasick_rs`."""

    def __init__(
        self,
        keywords: Iterable[str],
        *,
        case_sensitive: bool = True,
        requested_backend: BackendRequest = "ahocorasick_rs",
    ) -> None:
        import ahocorasick_rs

        self.requested_backend = requested_backend
        self.resolved_backend: BackendName = "ahocorasick_rs"
        self._keywords = _unique_keywords(keywords)
        self._case_sensitive = case_sensitive
        indexed_keywords = (
            self._keywords
            if case_sensitive
            else tuple(keyword.lower() for keyword in self._keywords)
        )
        self._automaton = ahocorasick_rs.AhoCorasick(indexed_keywords)

    @override
    def find_all(self, text: str) -> Iterator[Match]:
        """Yield all Rust matches in ``text``."""
        search_text = text if self._case_sensitive else text.lower()
        for index, start, end in self._automaton.find_matches_as_indexes(
            search_text,
            overlapping=True,
        ):
            yield Match(keyword=self._keywords[index], start=start, end=end)


def create_automaton(
    keywords: Iterable[str],
    backend: str = "auto",
    case_sensitive: bool = True,
) -> AhoCorasickBase:
    """Create a matcher and expose both requested and resolved backend names.

    Args:
        keywords: Peptide patterns to search for.
        backend: ``auto``, ``ahocorapy``, or ``ahocorasick_rs``.
        case_sensitive: Whether matching distinguishes letter case.

    Returns:
        A backend-neutral matcher.

    Raises:
        ImportError: The Rust backend was requested but is not installed.
        ValueError: ``backend`` is not supported.
    """
    requested_backend = _validate_backend(backend)
    resolved_backend = resolve_backend(requested_backend)
    if resolved_backend == "ahocorasick_rs":
        return AhoCorasickRust(
            keywords,
            case_sensitive=case_sensitive,
            requested_backend=requested_backend,
        )
    return AhoCorasickPure(
        keywords,
        case_sensitive=case_sensitive,
        requested_backend=requested_backend,
    )


def resolve_backend(backend: str = "auto") -> BackendName:
    """Resolve a requested backend to the concrete implementation name."""
    requested_backend = _validate_backend(backend)
    if requested_backend == "ahocorasick_rs":
        if not _rust_available():
            raise ImportError(
                "backend 'ahocorasick_rs' requires a working ahocorasick-rs installation"
            )
        return "ahocorasick_rs"
    if requested_backend == "ahocorapy":
        return "ahocorapy"
    return "ahocorasick_rs" if _rust_available() else "ahocorapy"


def get_available_backends() -> list[BackendName]:
    """Return the concrete matching backends available in this environment."""
    backends: list[BackendName] = ["ahocorapy"]
    if _rust_available():
        backends.append("ahocorasick_rs")
    return backends


def _validate_backend(backend: str) -> BackendRequest:
    if backend not in _VALID_BACKENDS:
        raise ValueError(f"backend must be one of {sorted(_VALID_BACKENDS)}, got {backend!r}")
    return cast(BackendRequest, backend)


def _rust_available() -> bool:
    return find_spec("ahocorasick_rs") is not None


def _unique_keywords(keywords: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(keywords))
