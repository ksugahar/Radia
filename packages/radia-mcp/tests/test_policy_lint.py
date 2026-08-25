# -*- coding: utf-8 -*-
"""Unit tests for tools/policy_lint.py (publish-boundary guard) on synthetic
fixtures -- independent of the live repo's current state, so CI stays green
while any real violation is remediated separately."""
import importlib.util
import subprocess
import sys
from pathlib import Path

_TOOLS = Path(__file__).resolve().parents[1] / "tools"
_spec = importlib.util.spec_from_file_location("policy_lint", _TOOLS / "policy_lint.py")
policy_lint = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(policy_lint)


def _make_repo(tmp, scripts, exclude=None, token_files=None, dependencies=None,
               optional_dependencies=None):
    src = tmp / "src" / "radia_mcp"
    src.mkdir(parents=True, exist_ok=True)
    py = ['[project]', 'name = "radia-mcp-test"', 'version = "0.0.0"']
    if dependencies is not None:
        py.append(f"dependencies = {list(dependencies)!r}")
    if optional_dependencies is not None:
        py += ['', '[project.optional-dependencies]']
        for key, values in optional_dependencies.items():
            py.append(f"{key} = {list(values)!r}")
    py += ['', '[project.scripts]']
    for k, v in scripts.items():
        py.append(f'{k} = "{v}"')
    if exclude is not None:
        py += ['', '[tool.setuptools.packages.find]', f'exclude = {list(exclude)!r}']
    (tmp / "pyproject.toml").write_text("\n".join(py) + "\n", encoding="utf-8")
    for v in scripts.values():
        sub = v.split(":")[0].split(".")[1]
        d = src / sub
        d.mkdir(parents=True, exist_ok=True)
        (d / "__init__.py").write_text("", encoding="utf-8")
        (d / "server.py").write_text("def main():\n    pass\n", encoding="utf-8")
    for relpath, content in (token_files or {}).items():
        f = src / relpath
        f.parent.mkdir(parents=True, exist_ok=True)
        f.write_text(content, encoding="utf-8")
    return tmp


def _run(root, strict=False):
    argv = ["policy_lint", "--root", str(root)] + (["--strict"] if strict else [])
    old = sys.argv
    sys.argv = argv
    try:
        return policy_lint.main()
    finally:
        sys.argv = old


def test_clean_open_servers_pass(tmp_path):
    _make_repo(tmp_path, {"mcp-server-cubit": "radia_mcp.cubit.server:main"})
    assert _run(tmp_path) == 0


def test_denied_wrapper_in_scripts_fails(tmp_path):
    _make_repo(tmp_path, {
        "mcp-server-cubit": "radia_mcp.cubit.server:main",
        "mcp-server-comsol-converter": "radia_mcp.comsol_converter.server:main",
    })
    assert _run(tmp_path) == 1   # wired public + shipped in wheel -> ERROR


def test_denied_excluded_and_unwired_warns_only(tmp_path):
    _make_repo(tmp_path, {"mcp-server-cubit": "radia_mcp.cubit.server:main"},
               exclude=["radia_mcp.comsol_converter*"])
    cc = tmp_path / "src" / "radia_mcp" / "comsol_converter"
    cc.mkdir(parents=True)
    (cc / "__init__.py").write_text("", encoding="utf-8")
    assert _run(tmp_path) == 0            # source present -> warning only
    assert _run(tmp_path, strict=True) == 1


def test_commercial_named_module_in_public_sub_warns(tmp_path):
    # a module NAMED for a commercial tool in a public sub -> review warning
    _make_repo(tmp_path, {"mcp-server-motor": "radia_mcp.motor.server:main"},
               token_files={"motor/comsol_bridge.py": "# bridge\n"})
    assert _run(tmp_path) == 0
    assert _run(tmp_path, strict=True) == 1


def test_benchmark_mention_not_flagged(tmp_path):
    # a plain mention of .mph / JMAG-Designer in a normally-named module is OK
    # for BOTH guards (structure + provenance): no tool-attributed number, no
    # "vs <tool>", no internal path -- just a bare format mention.
    _make_repo(tmp_path, {"mcp-server-motor": "radia_mcp.motor.server:main"},
               token_files={"motor/notes.py": "# supports the .mph and JMAG-Designer formats\n"})
    assert _run(tmp_path) == 0
    assert _run(tmp_path, strict=True) == 0   # mention alone is not a violation


def test_optuna_server_must_stay_external(tmp_path):
    _make_repo(tmp_path, {
        "mcp-server-optuna": "radia_mcp.optuna.server:main",
    })
    assert _run(tmp_path) == 1


def test_optuna_dependency_must_stay_external(tmp_path):
    _make_repo(tmp_path, {"mcp-server-cubit": "radia_mcp.cubit.server:main"},
               dependencies=["optuna>=4"],
               optional_dependencies={"external": ["optuna-mcp>=0.2"]})
    assert _run(tmp_path) == 1


def test_optuna_import_must_stay_external(tmp_path):
    _make_repo(tmp_path, {"mcp-server-motor": "radia_mcp.motor.server:main"},
               token_files={
                   "motor/objective.py": "import optuna\n",
                   "motor/dashboard.py": "from optuna_dashboard import run_server\n",
               })
    assert _run(tmp_path) == 1


def test_optuna_import_scan_accepts_utf8_bom_python_files(tmp_path):
    _make_repo(tmp_path, {"mcp-server-motor": "radia_mcp.motor.server:main"},
               token_files={"motor/bom_ok.py": "\ufeff# UTF-8 BOM is allowed\nx = 1\n"})
    assert _run(tmp_path) == 0


def test_optuna_matlab_difference_support_is_allowed(tmp_path):
    _make_repo(
        tmp_path,
        {"mcp-server-radia-matlab": "radia_mcp.matlab.server:main"},
        token_files={
            "matlab/optuna_boundary.py": (
                "def matlab_optuna_mcp_route():\n"
                "    return {'owner': 'optuna/optuna-mcp', "
                "'difference': 'MATLAB table/MAT and Simulink'}\n"
            ),
        },
    )
    assert _run(tmp_path) == 0


def test_tracked_only_provenance_scan_for_repo_root(tmp_path):
    examples = tmp_path / "examples"
    examples.mkdir()
    tracked = examples / "public_leak.py"
    tracked.write_text(r'PATH = r"W:\00_CAE\COMSOL\_crossval\case.py"' + "\n",
                       encoding="utf-8")
    untracked = examples / "private_scratch.py"
    untracked.write_text(r'PATH = r"S:\ELF_MAGIC\_crossval\case.meg"' + "\n",
                         encoding="utf-8")

    subprocess.run(["git", "init"], cwd=tmp_path, check=True,
                   stdout=subprocess.DEVNULL)
    subprocess.run(["git", "add", "examples/public_leak.py"], cwd=tmp_path,
                   check=True, stdout=subprocess.DEVNULL)

    tracked_findings = policy_lint.scan_text_tree(tmp_path, tracked_only=True)
    assert [item[0] for item in tracked_findings] == ["examples/public_leak.py"]

    all_findings = policy_lint.scan_text_tree(tmp_path)
    assert {item[0] for item in all_findings} == {
        "examples/public_leak.py",
        "examples/private_scratch.py",
    }


def test_tracked_only_provenance_scan_uses_all_scan_suffixes(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    tracked = src / "public_leak.m"
    tracked.write_text("(* PATH = W:/00_CAE/NGSolve/leak *)\n", encoding="utf-8")

    subprocess.run(["git", "init"], cwd=tmp_path, check=True,
                   stdout=subprocess.DEVNULL)
    subprocess.run(["git", "add", "src/public_leak.m"], cwd=tmp_path,
                   check=True, stdout=subprocess.DEVNULL)

    tracked_findings = policy_lint.scan_text_tree(tmp_path, tracked_only=True)
    assert [item[0] for item in tracked_findings] == ["src/public_leak.m"]
