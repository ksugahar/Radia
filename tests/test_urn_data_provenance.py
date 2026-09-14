"""Synthetic research fixtures must not overwrite or impersonate measurements."""
from __future__ import annotations

import ast
import csv
import os
import io
import zipfile
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/universal_relaxation_network/data/real_world/nasa_battery"


@pytest.mark.parametrize(("filename", "function"), [
    ("generate_paper_figures.py", "load_nasa_battery_data"),
    ("demo_spice_timedomain.py", "load_battery_data"),
])
def test_historical_nasa_consumers_reject_unverified_frequency(tmp_path, filename, function):
    path = DATA.parents[2] / filename
    node = next(n for n in ast.parse(path.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == function)
    namespace = {"Path": Path, "__file__": str(path), "np": np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    fixture = tmp_path / 'unverified.csv'
    fixture.write_text('frequency_Hz,Z_real_Ohm,Z_imag_Ohm\n1,2,3\n2,3,4\n')
    with pytest.raises(ValueError, match="frequency axis is unverified"):
        namespace[function](fixture)


def test_documented_nasa_figure_loader_parses_columns_not_fixed_header_lines(tmp_path):
    path = DATA.parents[2] / "generate_paper_figures.py"
    node = next(n for n in ast.parse(path.read_text(encoding="utf-8")).body
                if isinstance(n, ast.FunctionDef) and n.name == "load_nasa_battery_data")
    destination = tmp_path / "data/real_world/nasa_battery"
    destination.mkdir(parents=True)
    (destination / "nasa_18650_eis.csv").write_text(
        "# documented fixture\n#   Frequency source: instrument.csv:rows2-3\n"
        "frequency_Hz,Z_real_Ohm,Z_imag_Ohm\n5000,1,2\n1,3,4\n", encoding="utf-8")
    namespace = {"Path": Path, "__file__": str(tmp_path / path.name), "np": np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), namespace)
    frequency, impedance = namespace["load_nasa_battery_data"]()
    np.testing.assert_array_equal(frequency, [5000, 1])
    np.testing.assert_array_equal(impedance, [1+2j, 3+4j])


@pytest.mark.parametrize(('filename', 'function'), [
    ('generate_paper_figures.py', 'load_nasa_battery_data'),
    ('demo_spice_timedomain.py', 'load_battery_data'),
])
def test_nasa_consumers_accept_private_path_and_environment(tmp_path, monkeypatch, filename, function):
    path = DATA.parents[2] / filename
    node = next(n for n in ast.parse(path.read_text(encoding='utf-8')).body
                if isinstance(n, ast.FunctionDef) and n.name == function)
    namespace = {'Path': Path, '__file__': str(tmp_path / filename), 'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    private = tmp_path / 'private.csv'
    private.write_text('#   Frequency source: synthetic test only\n'
                       'frequency_Hz,Z_real_Ohm,Z_imag_Ohm\n5,1,2\n1,3,4\n')
    monkeypatch.setenv('RADIA_NASA_EIS_CSV', str(tmp_path / 'missing.csv'))
    np.testing.assert_array_equal(namespace[function](private)[0], [5, 1])
    monkeypatch.setenv('RADIA_NASA_EIS_CSV', str(private))
    np.testing.assert_array_equal(namespace[function]()[1], [1+2j, 3+4j])
    monkeypatch.delenv('RADIA_NASA_EIS_CSV')
    if function == 'load_battery_data':
        with pytest.raises(FileNotFoundError, match='RADIA_NASA_EIS_CSV'):
            namespace[function]()
    else:
        assert namespace[function]() == (None, None)


@pytest.mark.parametrize(('filename', 'function'), [
    ('generate_paper_figures.py', 'load_nasa_battery_data'),
    ('demo_spice_timedomain.py', 'load_battery_data'),
])
@pytest.mark.parametrize('rows', [
    '1,1,2\n', '0,1,2\n2,3,4\n', '1,1,2\n1,3,4\n',
    'nan,1,2\n2,3,4\n', '1,nan,2\n2,3,4\n', '1,1,inf\n2,3,4\n',
])
def test_documented_nasa_consumers_still_validate_numerical_data(tmp_path, filename, function, rows):
    path = DATA.parents[2] / filename
    node = next(n for n in ast.parse(path.read_text(encoding='utf-8')).body
                if isinstance(n, ast.FunctionDef) and n.name == function)
    destination = tmp_path / 'data/real_world/nasa_battery'
    destination.mkdir(parents=True)
    (destination / 'nasa_18650_eis.csv').write_text(
        '#   Frequency source: instrument.csv\nfrequency_Hz,Z_real_Ohm,Z_imag_Ohm\n'+rows,
        encoding='utf-8')
    namespace = {'Path': Path, '__file__': str(tmp_path / filename), 'np': np}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), namespace)
    with pytest.raises(ValueError, match='positive finite'):
        namespace[function]()


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


def _extractor_namespace():
    path = DATA / "extract_real_eis.py"
    tree = ast.parse(path.read_text(encoding="utf-8"))
    nodes = [n for n in tree.body if isinstance(n, ast.FunctionDef)]
    namespace = {"np": np, "os": os, "io": io, "zipfile": zipfile}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


@pytest.mark.parametrize(("values", "source"), [
    (None, "instrument log"), ([1, 2], ""), ([1, 2], "\nforged"),
    ([1], "log"), ([0, 2], "log"), ([1, np.nan], "log"),
    ([1, np.inf], "log"), ([1, 1], "log"), ([[1, 2]], "log"),
])
def test_measured_frequency_rejects_undocumented_or_invalid_axis(values, source):
    with pytest.raises(ValueError):
        _extractor_namespace()["validate_frequency_axis"](values, 2, source)


def test_measured_frequency_preserves_instrument_sample_order():
    actual = _extractor_namespace()["validate_frequency_axis"]([5000, 1, 100], 3, "instrument.csv:rows 2-4")
    np.testing.assert_array_equal(actual, [5000, 1, 100])


def test_real_extraction_never_extracts_archive_members_or_invents_frequency(tmp_path):
    import scipy.io

    data = {"Battery_impedance": np.array([1 + 2j, 3 + 4j])}
    stream = io.BytesIO()
    cycles = np.empty((1, 1), dtype=[("type", "O"), ("data", "O")])
    cycles[0, 0] = ("impedance", data)
    scipy.io.savemat(stream, {"B0005": {"cycle": cycles}})
    with zipfile.ZipFile(tmp_path / "1. BatteryAgingARC-FY08Q4.zip", "w") as archive:
        archive.writestr("nested/B0005.mat", stream.getvalue())
    sentinel = tmp_path / "B0005.mat"
    sentinel.write_bytes(b"caller-owned")
    namespace = _extractor_namespace()
    namespace.update(NASA_DATA_DIR=str(tmp_path), scipy=SimpleNamespace(io=scipy.io))
    extract = namespace["extract_eis_from_nasa"]
    with pytest.raises(ValueError, match="documented frequency"):
        extract()
    result = extract(frequency_hz=[5000, 1], frequency_source="instrument record")
    np.testing.assert_array_equal(result["frequency"], [5000, 1])
    np.testing.assert_array_equal(result["Z"], [1 + 2j, 3 + 4j])
    assert sentinel.read_bytes() == b"caller-owned"
    assert not (tmp_path / "nested").exists()
    with pytest.raises(ValueError, match="cycle_type"):
        extract(cycle_type="typo")
