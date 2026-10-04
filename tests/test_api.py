"""The public module other anndata_bridge packages import."""

from __future__ import annotations

from prozor import api


def test_api_exports_exactly_the_approved_names() -> None:
    assert sorted(api.__all__) == [
        "annotate_peptides",
        "resolve_backend",
    ]
