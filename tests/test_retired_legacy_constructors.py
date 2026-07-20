"""Unsafe legacy extruded constructors must not return to any API surface."""

import ctypes
import importlib
from pathlib import Path
import sys

import pytest
import radia


RETIRED_PYTHON_NAMES = ("ObjMltExtPgn", "ObjMltExtRtg", "ObjMltExtTri")
RETIRED_C_ABI_SYMBOLS = ("RadObjMltExtPgn", "RadObjMltExtRtg", "RadObjMltExtTri")


def test_unsafe_legacy_extruded_constructors_are_not_public():
    for name in RETIRED_PYTHON_NAMES:
        assert not hasattr(radia, name)


def test_unsafe_legacy_extruded_constructors_are_absent_from_c_abi_sources():
    root = Path(__file__).resolve().parents[1]
    for relative_path in (
        Path("src/lib/radentry.h"),
        Path("src/lib/radentry.cpp"),
        Path("src/lib/raddll.def"),
    ):
        source = (root / relative_path).read_text(encoding="utf-8", errors="ignore")
        for symbol in RETIRED_C_ABI_SYMBOLS:
            assert symbol not in source, f"{symbol} leaked into {relative_path}"


@pytest.mark.skipif(sys.platform != "win32", reason="Windows DLL export check")
def test_unsafe_legacy_extruded_constructors_are_not_exported():
    extension = importlib.import_module("radia._radia_pybind")
    root = Path(__file__).resolve().parents[1]
    build_artifacts = list(root.glob("build*/_radia_pybind.cp*-win_amd64.pyd"))
    candidates = [Path(extension.__file__), *build_artifacts]
    binary = max(candidates, key=lambda path: path.stat().st_mtime_ns)
    dll = ctypes.WinDLL(str(binary))
    for symbol in RETIRED_C_ABI_SYMBOLS:
        with pytest.raises(AttributeError):
            getattr(dll, symbol)
