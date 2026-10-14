import pytest
from radia_mcp.md2html import md_to_html


@pytest.mark.parametrize("optional_args", [0, 1, 2])
def test_legacy_cli_delegates_to_package(tmp_path, optional_args):
    import os
    from pathlib import Path
    import subprocess
    import sys

    root = Path(__file__).resolve().parents[2]
    source = tmp_path / "input file.md"
    source.write_text("# Heading\n\n```math\nx^2\n```\n", encoding="utf-8")
    output = tmp_path / "explicit output.html" if optional_args else source.with_suffix(".html")
    args = [str(source)]
    if optional_args:
        args.append(str(output))
    if optional_args == 2:
        args.append("Custom title")
    env = dict(os.environ, PYTHONPATH=str(root / "packages/radia-mcp/src"))
    run = subprocess.run(
        [sys.executable, str(root / ".agents/skills/md2html/md2html.py"), *args],
        cwd=tmp_path, env=env, capture_output=True, text=True, timeout=30)
    assert run.returncode == 0, run.stderr
    assert str(output) in run.stdout
    expected = tmp_path / "expected.html"
    md_to_html(str(source), str(expected), "Custom title" if optional_args == 2 else None)
    assert output.read_bytes() == expected.read_bytes()
    if optional_args:
        assert not source.with_suffix(".html").exists()
