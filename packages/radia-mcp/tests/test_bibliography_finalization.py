"""Compiled citation decisions and cross-process writer exclusion."""
import subprocess
import sys
from types import SimpleNamespace

import pytest

from radia_mcp.bibliography._compiled_aux import read_compiled_aux
from radia_mcp.bibliography._write_lock import target_lock
from radia_mcp.bibliography import _source_edit, _write_lock


def test_aux_resolves_compiled_citation_order_and_nested_inputs(tmp_path):
    (tmp_path / "child.aux").write_text(r"\citation{two,one}\bibstyle{plain}", encoding="utf-8")
    main = tmp_path / "paper.aux"
    main.write_text("% \\citation{fake}\n" + r"\citation{one}\@input{child.aux}\citation{*}", encoding="utf-8")
    keys, style, snapshots = read_compiled_aux(main)
    assert keys == ["one", "two", "*"] and style == "plain"
    assert len(snapshots) == 2 and all(len(s["sha256"]) == 64 for s in snapshots)


@pytest.mark.parametrize("body", [r"\citation{one}\citation{broken", r"\@input{missing.aux}",
    r"\citation{one}\@input{paper.aux}", r"\citation{one}\bibstyle{plain}\bibstyle{IEEEtran}",
    r"\abx@aux@cite{0}{one}"])
def test_invalid_aux_is_not_a_partial_success(tmp_path, body):
    main = tmp_path / "paper.aux"
    main.write_text(body, encoding="utf-8")
    with pytest.raises((ValueError, OSError)):
        read_compiled_aux(main)


@pytest.mark.parametrize("changed", [False, True])
def test_bbl_uses_compiled_decisions_instead_of_static_macros(tmp_path, monkeypatch, changed):
    from radia_mcp.bibliography.plans import T14_canonical as canonical
    bib = tmp_path / "fixture.bib"
    bib.write_text("@misc{selected,title={Selected}}", encoding="utf-8")
    monkeypatch.setattr(canonical, "CANONICAL", bib)
    monkeypatch.setattr(canonical.shutil, "which", lambda _: "bibtex")
    (tmp_path / "plain.bst").write_bytes(b"test style")
    tex = tmp_path / "paper.tex"
    tex.write_text(r"\newcommand{\mycite}{\cite{selected}}\mycite\iffalse\cite{excluded}\fi", encoding="utf-8")
    aux = tmp_path / "paper.aux"
    aux.write_text(r"\citation{selected}\bibstyle{plain}", encoding="utf-8")
    out = tmp_path / "paper.bbl"
    out.write_bytes(b"previous")
    def run(*args, cwd, **kwargs):
        assert r"\citation{selected}" in (cwd / "manuscript.aux").read_text()
        assert "excluded" not in (cwd / "manuscript.aux").read_text()
        (cwd / "manuscript.bbl").write_bytes(br"\bibitem{selected}Selected")
        if changed:
            aux.write_text(r"\citation{other}", encoding="utf-8")
        return SimpleNamespace(returncode=0, stdout=b"", stderr=b"")
    monkeypatch.setattr(canonical.subprocess, "run", run)
    result = canonical.bibliography_make_bbl(str(tex), aux_path=str(aux))
    if changed:
        assert result.startswith("Error:") and out.read_bytes() == b"previous"
    else:
        assert result.startswith("bibliography_make_bbl:") and "compiled aux" in result


def test_lock_excludes_another_process_and_releases_on_exception(tmp_path):
    path = tmp_path / "target.bib"
    code = ("from pathlib import Path; import runpy,sys; "
            "target_lock=runpy.run_path(sys.argv[2])['target_lock']; "
            "ctx=target_lock(Path(__import__('sys').argv[1])); ctx.__enter__(); ctx.__exit__(None,None,None)")
    # Test this exact implementation, not editable-package startup over SMB.
    command = [sys.executable, "-S", "-c", code, str(path), _write_lock.__file__]
    with pytest.raises(RuntimeError):
        with target_lock(path):
            result = subprocess.run(command, capture_output=True, timeout=20)
            assert result.returncode != 0 and b"target busy" in result.stderr
            raise RuntimeError("operation failed")
    result = subprocess.run(command, capture_output=True, timeout=20)
    assert result.returncode == 0, result.stderr
    assert path.with_name(".target.bib.radia-write-lock").exists()


def test_locked_source_is_not_rewritten(tmp_path):
    path = tmp_path / "fixture.bib"
    path.write_bytes(b"@misc{old,title={First}}")
    original, text, entries = _source_edit.read_source(path)
    with target_lock(path):
        with pytest.raises(OSError, match="busy"):
            _source_edit.write_source_edits(path, original, text, [(*entries[0].key_span, "new")])
    assert path.read_bytes() == original
