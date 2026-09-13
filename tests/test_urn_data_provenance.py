"""Synthetic research fixtures must not overwrite or impersonate measurements."""
from __future__ import annotations

import ast
import csv
import os
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/universal_relaxation_network/data/real_world/nasa_battery"


class _DataFrame:
    """Exercise generator I/O without importing its download or pandas setup."""

    def __init__(self, columns):
        self.columns = columns

    def to_csv(self, stream, index=False):
        assert index is False
        writer = csv.writer(stream)
        writer.writerow(self.columns)
        writer.writerows(zip(*self.columns.values()))


@pytest.mark.parametrize(("filename", "function", "outputs"), [
    ("extract_nasa_eis.py", "create_representative_eis_csv", ["synthetic_18650_eis.csv"]),
    ("download_nasa_data.py", "create_eis_csv", ["synthetic_18650_eis_fresh.csv", "synthetic_18650_eis_aged.csv"]),
])
def test_synthetic_generators_label_data_and_preserve_measurements(tmp_path, filename, function, outputs):
    protected = ["nasa_18650_eis.csv", "nasa_18650_eis_fresh.csv", "nasa_18650_eis_aged.csv", "README.md"]
    for name in protected:
        (tmp_path / name).write_text("measurement-owned sentinel", encoding="utf-8")
    path = DATA / filename
    tree = ast.parse(path.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == function)
    namespace = {"np": np, "pd": SimpleNamespace(DataFrame=_DataFrame), "os": os}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    namespace[function](str(tmp_path))
    for name in protected:
        assert (tmp_path / name).read_text(encoding="utf-8") == "measurement-owned sentinel"
    assert set(p.name for p in tmp_path.iterdir()) == set(protected + outputs)
    for name in outputs:
        text = (tmp_path / name).read_text(encoding="utf-8")
        assert text.startswith("# SYNTHETIC")
        assert "REAL MEASUREMENT DATA" not in text
        rows = list(csv.reader(line for line in text.splitlines() if line and not line.startswith("#")))
        assert len(rows) == 51
        assert rows[0][0] == "frequency_Hz"


def test_synthetic_documentation_does_not_replace_measured_readme(tmp_path):
    path = DATA / "extract_nasa_eis.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    node = next(n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "create_readme")
    (tmp_path / "README.md").write_text("measurement-owned sentinel", encoding="utf-8")
    namespace = {"os": os}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    namespace["create_readme"](str(tmp_path))
    assert (tmp_path / "README.md").read_text(encoding="utf-8") == "measurement-owned sentinel"
    assert "not extracted from NASA" in (tmp_path / "SYNTHETIC_DATA.md").read_text(encoding="utf-8")
