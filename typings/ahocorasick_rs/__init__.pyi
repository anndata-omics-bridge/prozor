from collections.abc import Sequence

class AhoCorasick:
    def __init__(self, patterns: Sequence[str]) -> None: ...
    def find_matches_as_indexes(
        self,
        haystack: str,
        *,
        overlapping: bool = False,
    ) -> list[tuple[int, int, int]]: ...
