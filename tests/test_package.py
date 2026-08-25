from __future__ import annotations

import importlib
from pathlib import Path


def test_package_is_importable() -> None:
    """The generated distribution exposes its configured import package."""
    assert importlib.import_module("prozor") is not None


def test_package_markers_are_empty() -> None:
    package = Path(__file__).parents[1] / "src" / "prozor"
    for marker in (
        package / "__init__.py",
        package / "matching/__init__.py",
        package / "inference/__init__.py",
    ):
        assert not marker.read_text().strip()
