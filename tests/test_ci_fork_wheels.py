import hashlib
import importlib.util
from io import BytesIO
import json
from pathlib import Path
import sys
from types import SimpleNamespace
import zipfile

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ci_fork_wheels", ROOT / "tools/ci_fork_wheels.py")
fork = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fork)


def test_modified_wheel_cannot_reach_pip(tmp_path, monkeypatch):
    filename = "netgen_mesher-test.whl"
    url = "https://example.invalid/netgen_mesher-test.whl"
    monkeypatch.setattr(fork, "urlopen", lambda url, timeout: BytesIO(b"modified"))
    lock = {"wheels": {"netgen-mesher": {
        "filename": filename, "url": url, "sha256": hashlib.sha256(b"verified").hexdigest()
    }}}
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    monkeypatch.setattr(fork.subprocess, "run", lambda *a, **k: pytest.fail("pip received an unverified wheel"))
    with pytest.raises(RuntimeError, match="Incorrect wheel SHA-256"):
        fork.install(lock)


def test_lock_identifies_both_cpu_wheels_and_native_dlls():
    from packaging.specifiers import SpecifierSet
    import tomllib

    lock = json.loads(fork.LOCK.read_text(encoding="utf-8"))
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    assert set(lock["wheels"]) == {"netgen-mesher", "ngsolve"}
    assert set(lock["dlls"]) == {"nglib.dll", "ngcore.dll", "libngsolve.dll"}
    for package, wheel in lock["wheels"].items():
        pin, = [dep[len(package):] for dep in deps if dep.startswith(package + "==")]
        assert SpecifierSet(pin).contains(wheel["version"])
        expected = lock["sources"]["netgen" if package == "netgen-mesher" else "ngsolve"]
        assert "+radia." + expected[:8] in wheel["filename"]
    assert lock["wheels"]["ngsolve"]["version"].endswith(".cpu")


@pytest.mark.parametrize("outside", ["netgen", "ngsolve"])
def test_verify_rejects_packages_outside_the_venv(tmp_path, monkeypatch, outside):
    from types import SimpleNamespace

    prefix = tmp_path / "venv"
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.setattr(fork.importlib.util, "find_spec",
                        lambda name: SimpleNamespace(submodule_search_locations=[
                            str((tmp_path if name == outside else prefix) / name)
                        ]))
    with pytest.raises(RuntimeError, match=outside + " is outside the CI environment"):
        fork.verify_installed({"wheels": {}, "dlls": {}, "sources": {}})


def _installed(tmp_path, monkeypatch, versions, dll=b"good"):
    """A fake isolated venv holding netgen/ngsolve packages and one locked DLL."""
    from types import SimpleNamespace

    prefix = tmp_path / "venv"
    for name in ("netgen", "ngsolve"):
        (prefix / name).mkdir(parents=True)
    (prefix / "netgen" / "nglib.dll").write_bytes(dll)
    monkeypatch.setattr(sys, "prefix", str(prefix))
    monkeypatch.setattr(fork.importlib.metadata, "version", versions.__getitem__)
    monkeypatch.setattr(fork.importlib.util, "find_spec",
                        lambda name: SimpleNamespace(submodule_search_locations=[str(prefix / name)]))
    return {"wheels": {"ngsolve": {"version": "6.2.2608+radia.f2223d8b1.cpu"}},
            "dlls": {"nglib.dll": hashlib.sha256(b"good").hexdigest()}, "sources": {"ok": True}}


def test_verify_rejects_a_stock_or_other_variant_version(tmp_path, monkeypatch):
    """Guards: a dependency install that swaps in stock or CUDA NGSolve fails the lane."""
    for wrong in ("6.2.2608", "6.2.2608+radia.f2223d8b1.cuda"):
        lock = _installed(tmp_path / wrong, monkeypatch, {"ngsolve": wrong})
        with pytest.raises(RuntimeError, match="Incorrect installed version: ngsolve"):
            fork.verify_installed(lock)


def test_verify_rejects_a_replaced_native_library(tmp_path, monkeypatch):
    """Guards: the installed Netgen DLL bytes are bound to the locked fork build."""
    lock = _installed(tmp_path, monkeypatch, {"ngsolve": "6.2.2608+radia.f2223d8b1.cpu"}, dll=b"stock")
    with pytest.raises(RuntimeError, match="Incorrect installed DLL: nglib.dll"):
        fork.verify_installed(lock)


def test_verify_accepts_the_locked_installation(tmp_path, monkeypatch):
    lock = _installed(tmp_path, monkeypatch, {"ngsolve": "6.2.2608+radia.f2223d8b1.cpu"})
    assert fork.verify_installed(lock) == {"ok": True}


def test_a_rejected_download_leaves_no_wheel_behind(tmp_path, monkeypatch):
    """Guards: a corrupted download is not kept where a later run could reuse it."""
    monkeypatch.setattr(fork.tempfile, "tempdir", str(tmp_path))
    monkeypatch.setattr(fork, "urlopen", lambda url, timeout: BytesIO(b"modified"))
    monkeypatch.setattr(sys, "prefix", str(tmp_path / "venv"))
    monkeypatch.setattr(fork.subprocess, "run", lambda *a, **k: pytest.fail("pip ran"))
    lock = {"wheels": {"netgen-mesher": {"filename": "w.whl", "url": "u",
                                         "sha256": hashlib.sha256(b"verified").hexdigest()}}}
    with pytest.raises(RuntimeError, match="Incorrect wheel SHA-256"):
        fork.install(lock)
    assert list(tmp_path.glob("radia-fork-wheels-*")) == []


def test_install_refuses_the_base_interpreter(monkeypatch):
    """Guards: fork wheels are never installed into a shared base Python."""
    monkeypatch.setattr(sys, "prefix", sys.base_prefix)
    with pytest.raises(RuntimeError, match="isolated CI venv"):
        fork.install({"wheels": {}})


def test_lock_urls_name_the_locked_release_and_file():
    import urllib.parse

    lock = json.loads(fork.LOCK.read_text(encoding="utf-8"))
    for package, wheel in lock["wheels"].items():
        source = lock["sources"]["netgen" if package == "netgen-mesher" else "ngsolve"]
        assert wheel["url"].endswith("/" + urllib.parse.quote(wheel["filename"], safe="")), wheel["url"]
        assert f"/releases/download/radia-win-{source[:8]}" in wheel["url"], wheel["url"]

def test_linux_fork_requires_the_locked_wheel_and_installed_payload(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("gpu_builder", ROOT / "tools/build_python_gpu.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    site = tmp_path / "venv" / "site"
    lock = {"schema": 1, "wheels": {}}
    origins = {}
    for package, module, source in (
        ("netgen-mesher", "netgen", "a" * 40),
        ("ngsolve", "ngsolve", "b" * 40),
    ):
        entries = {f"{module}/__init__.py": b"", f"{module}/native.so": module.encode()}
        wheel = tmp_path / f"{module}.whl"
        with zipfile.ZipFile(wheel, "w") as archive:
            for name, content in entries.items():
                archive.writestr(name, content)
                path = site / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(content)
        lock["wheels"][package] = {
            "file": wheel.name,
            "sha256": builder.sha(wheel),
            "source_commit": source,
            "version": f"6.2.2608+radia.{source[:8]}.cuda",
        }
        origins[module] = site / module / "__init__.py"
    path = tmp_path / "lock.json"
    path.write_text(json.dumps(lock), encoding="utf-8")
    monkeypatch.setattr(builder.sys, "prefix", str(site.parent))
    monkeypatch.setattr(builder.importlib.util, "find_spec", lambda name: SimpleNamespace(origin=origins[name]))
    monkeypatch.setattr(builder.metadata, "distribution", lambda name: SimpleNamespace(
        version=lock["wheels"][name]["version"], locate_file=lambda relative: site / relative
    ))
    assert builder.verify_linux_fork(path) == lock
    native = site / "ngsolve/native.so"
    native.write_bytes(b"wrong ABI")
    with pytest.raises(RuntimeError, match="installed|Installed"):
        builder.verify_linux_fork(path)
    native.write_bytes(b"ngsolve")
    (tmp_path / "ngsolve.whl").write_bytes(b"wrong wheel")
    with pytest.raises(RuntimeError, match="wheel hash mismatch"):
        builder.verify_linux_fork(path)


def test_matlab_consumer_uses_the_producer_fork_and_checks_abi_before_parity():
    workflow = ROOT / ".github/workflows/radia-simulink-a100.yml"
    steps = yaml.safe_load(workflow.read_text(encoding="utf-8"))["jobs"]["validate"]["steps"]
    commands = "\n".join(step.get("run", "") for step in steps)
    assert (commands.index("ci_fork_wheels.py install")
            < commands.index("ci_fork_wheels.py verify")
            < commands.index("Get-FileHash -LiteralPath $native.FullName")
            < commands.index("import radia._radia_pybind")
            < commands.index("run_device_boundary.py"))
    assert "pip install @pins" not in commands


def test_ltspice_acceptance_checks_both_checkouts_and_uploads_the_bound_log():
    workflow = ROOT / ".github/workflows/ltspice-data-safety.yml"
    document = yaml.safe_load(workflow.read_text(encoding="utf-8"))
    assert document[True]["workflow_dispatch"]["inputs"]["sha"]["required"]
    for job in ("ci-impact", "matlab-data-contracts"):
        steps = document["jobs"][job]["steps"]
        checkout = next(step for step in steps if step.get("uses", "").startswith("actions/checkout@"))
        assert "ACCEPTANCE_SHA" in checkout["with"]["ref"]
        assert any("git rev-parse HEAD" in step.get("run", "") for step in steps)
    steps = document["jobs"]["matlab-data-contracts"]["steps"]
    commands = "\n".join(step.get("run", "") for step in steps)
    assert '"source_commit": os.environ["ACCEPTANCE_SHA"]' in commands
    assert '"log_sha256": hashlib.sha256(output).hexdigest()' in commands
    assert any(step.get("with", {}).get("name") == "ltspice-acceptance" for step in steps)
