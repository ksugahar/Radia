"""The patched Netgen library reaches a run-local environment only through its receipt.

Guards: the patch changes the pyramid scale and nothing else; a changed library, a receipt from
another patch, a shared interpreter, a package outside the environment or another Netgen version
is refused; install replaces nglib.dll alone; the installed lane refuses a stock library.
Assumes: nothing is compiled here. That the patched library meets the 1e-11 criterion is shown by
test_hdiv_pyramid_material.py in the GPU lane.
"""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


overlay_tool = _load("netgen_overlay", "tools/netgen_overlay.py")


def _overlay(tmp_path, library=b"patched"):
    root = tmp_path / "overlay"
    root.mkdir()
    (root / overlay_tool.LIBRARY).write_bytes(library)
    (root / overlay_tool.RECEIPT).write_text(json.dumps(dict(
        schema=overlay_tool.SCHEMA, upstream_sha=overlay_tool.UPSTREAM_SHA, patch_sha256=overlay_tool.digest(overlay_tool.PATCH),
        library=overlay_tool.LIBRARY, library_sha256=overlay_tool.digest(root / overlay_tool.LIBRARY))), encoding="utf-8")
    return root


def test_patch_changes_the_pyramid_scale_and_nothing_else():
    lines = overlay_tool.PATCH.read_text(encoding="utf-8").splitlines()
    files = [line for line in lines if line.startswith("diff --git")]
    removed = [line[1:].strip() for line in lines if line.startswith("-") and not line.startswith("---")]
    assert files == ["diff --git a/libsrc/meshing/curvedelems.cpp b/libsrc/meshing/curvedelems.cpp"], files
    assert removed == ["z *= (1-1e-10);"], removed
    assert sum(line.startswith("@@") for line in lines) == 1, lines


@pytest.mark.parametrize("field,value,message", [("library_sha256", "0"*64, "Overlay library"), ("patch_sha256", "0"*64, "Overlay patch"),
                                                 ("upstream_sha", "0"*40, "expected radia.netgen-overlay.v1")])
def test_a_changed_overlay_is_refused(tmp_path, field, value, message):
    root = _overlay(tmp_path)
    receipt = json.loads((root / overlay_tool.RECEIPT).read_text(encoding="utf-8"))
    receipt[field] = value
    (root / overlay_tool.RECEIPT).write_text(json.dumps(receipt), encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        overlay_tool.install(root, "unused")


@pytest.mark.parametrize("found,message", [
    (dict(package="/env/netgen", prefix="/base", base="/base", version="6.2.2608"), "not a virtual environment"),
    (dict(package="/base/netgen", prefix="/env", base="/base", version="6.2.2608"), "outside"),
    (dict(package="/env/netgen", prefix="/env", base="/base", version="6.2.2607"), "netgen-mesher 6.2.2607")])
def test_only_the_pinned_netgen_of_a_virtual_environment_is_touched(monkeypatch, found, message):
    monkeypatch.setattr(overlay_tool.subprocess, "check_output", lambda *a, **k: json.dumps(found))
    with pytest.raises(ValueError, match=message):
        overlay_tool.locate("python")


def test_install_replaces_only_the_library(tmp_path, monkeypatch):
    root = _overlay(tmp_path)
    package = tmp_path / "netgen"
    package.mkdir()
    (package / overlay_tool.LIBRARY).write_bytes(b"stock")
    (package / "ngcore.dll").write_bytes(b"other")
    monkeypatch.setattr(overlay_tool, "locate", lambda python: package)
    record = overlay_tool.install(root, "unused")
    assert (package / overlay_tool.LIBRARY).read_bytes() == b"patched"
    assert (package / "ngcore.dll").read_bytes() == b"other"
    assert record["library_sha256"] == overlay_tool.digest(root / overlay_tool.LIBRARY) != record["previous_sha256"], record


def test_installed_lane_refuses_a_stock_library(tmp_path, monkeypatch):
    lane = _load("hdiv_lane_netgen", "validation_test/hdiv_cuda/run_installed.py")
    root = _overlay(tmp_path)
    package = tmp_path / "netgen"
    package.mkdir()
    (package / overlay_tool.LIBRARY).write_bytes(b"stock")
    monkeypatch.setattr("importlib.util.find_spec", lambda name: type("Spec", (), dict(submodule_search_locations=[str(package)])))
    monkeypatch.setattr(lane.sys, "prefix", str(tmp_path))
    with pytest.raises(RuntimeError, match="nglib.dll is"):
        lane.check_netgen(root)
    (package / overlay_tool.LIBRARY).write_bytes(b"patched")
    assert lane.check_netgen(root)["library_sha256"] == overlay_tool.digest(root / overlay_tool.LIBRARY)
