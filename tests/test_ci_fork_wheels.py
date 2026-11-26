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
    assert (commands.index("verify_public_fe_sdk.py install")
            < commands.index("verify_public_fe_sdk.py verify")
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


def test_hdiv_host_configure_uses_the_verified_cuda_toolkit(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("gpu_builder", ROOT / "tools/build_python_gpu.py")
    builder = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(builder)
    import shutil

    monkeypatch.setattr(shutil, "which", lambda executable: executable)
    monkeypatch.setattr(builder, "installed_config", lambda package, name: tmp_path / package / name)
    monkeypatch.setattr(builder, "cuda_libraries", lambda cuda, windows: {})

    class ConfigureObserved(Exception):
        pass

    cuda = tmp_path / "CUDA Toolkit" / "12.8"

    def verify_configure(command, out, label, commands):
        assert label == "_hdiv_cuda-host-configure"
        assert f"-DCUDAToolkit_ROOT={cuda.as_posix()}" in command
        raise ConfigureObserved

    monkeypatch.setattr(builder, "execute", verify_configure)
    args = SimpleNamespace(host_compiler=tmp_path / "compiler")
    with pytest.raises(ConfigureObserved):
        builder.build_hdiv_host(args, tmp_path / "ext.so", tmp_path / "device.o",
                                tmp_path / "include", cuda, tmp_path, [])


def test_matlab_child_uses_locked_mkl_directory_and_restores_parent_path(tmp_path, monkeypatch):
    import os

    spec = importlib.util.spec_from_file_location(
        "matlab_runtime", ROOT / "validation_test/ams_cuda/matlab_runtime.py"
    )
    runtime = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runtime)
    prefix, base = tmp_path / "venv", tmp_path / "python"
    net = prefix / "Lib/site-packages/netgen"
    mkl = prefix / "Library/bin/mkl_rt.3.dll"
    occt = prefix / "bin/TKOffset.dll"
    python = base / f"python{sys.version_info.major}{sys.version_info.minor}.dll"
    for path in (net / "__init__.py", mkl, occt, python):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.touch()
    lock = tmp_path / "fork-lock.json"
    lock.write_text("{}", encoding="utf-8")
    verified = []
    monkeypatch.setattr(runtime.runpy, "run_path", lambda path: {
        "LOCK": lock, "verify_installed": verified.append,
    })
    monkeypatch.setattr(runtime.sys, "prefix", str(prefix))
    monkeypatch.setattr(runtime.sys, "base_prefix", str(base))
    monkeypatch.setattr(runtime.importlib.util, "find_spec",
                        lambda name: SimpleNamespace(origin=str(net / "__init__.py")))
    files = {"mkl": mkl, "netgen-occt": occt}
    monkeypatch.setattr(runtime.metadata, "distribution", lambda name: SimpleNamespace(
        files=[files[name]], locate_file=lambda path: path
    ))
    monkeypatch.setenv("PATH", "original")

    child = {"PATH": "original"}

    def setenv(name, value, **kwargs):
        child[name] = value
        return _matlab_future(None)

    def launch(options, *, background):
        assert os.environ["PATH"] == "original" and background
        return _matlab_future(SimpleNamespace(
            getenv=lambda name, **kwargs: _matlab_future(child[name]), setenv=setenv
        ))

    engine = SimpleNamespace(start_matlab=launch)
    runtime.start_matlab(engine, "options")
    assert child["PATH"].split(os.pathsep) == [
        str(net), str(mkl.parent), str(occt.parent), str(base), "original"
    ]
    assert verified == [{}]
    assert os.environ["PATH"] == "original"
    mkl.unlink()
    with pytest.raises(FileNotFoundError, match=r"mkl_rt\.3\.dll"):
        runtime.start_matlab(engine, "options")
    assert os.environ["PATH"] == "original"


def test_public_wheel_gate_rejects_local_fork_versions(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("public_sdk", ROOT / "tools/verify_public_fe_sdk.py")
    sdk = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(sdk)
    monkeypatch.setattr(sdk.sys, "prefix", str(tmp_path))
    monkeypatch.setattr(sdk.metadata, "version", lambda name: "6.2.2608+radia.fork")
    with pytest.raises(RuntimeError, match="required public version"):
        sdk.verify()


def test_release_wheels_exercise_pip_resolution_in_fresh_venvs():
    workflow = ROOT / ".github/workflows/build-test.yml"
    jobs = yaml.safe_load(workflow.read_text(encoding="utf-8"))["jobs"]
    for job, gpu in (("build-test", False), ("windows-gpu-wheel", True)):
        commands = "\n".join(step.get("run", "") for step in jobs[job]["steps"])
        assert "ci_fork_wheels.py" not in commands
        gate = next(step["run"] for step in jobs[job]["steps"]
                    if step.get("name") == "Verify exact wheel through clean public pip resolution")
        assert "Expected a fresh acceptance path" in gate
        assert "verify_public_fe_sdk.py smoke --wheel" in gate
        assert ("--gpu" in gate) == gpu


def _matlab_future(value=None, *, error=None):
    def result(*, timeout):
        assert timeout > 0
        if error is not None:
            raise error
        return value

    return SimpleNamespace(result=result, cancel=lambda: None)


@pytest.mark.parametrize("phase", ["startup", "mex"])
def test_matlab_timeout_only_kills_the_owned_process(tmp_path, monkeypatch, capsys, phase):
    spec = importlib.util.spec_from_file_location(
        "matlab_runtime", ROOT / "validation_test/ams_cuda/matlab_runtime.py"
    )
    runtime = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runtime)
    killed = []

    class Process:
        def __init__(self, pid):
            self.pid = pid

        def kill(self):
            killed.append(self.pid)

    ours, other = Process(123), Process(999)
    class EngineTimeout(Exception):
        pass

    monkeypatch.setitem(sys.modules, "matlab.engine", SimpleNamespace(TimeoutError=EngineTimeout))
    stalled = _matlab_future(error=EngineTimeout())
    if phase == "startup":
        children = iter(({other}, {ours, other}))
        monkeypatch.setattr(runtime, "_matlab_children", lambda: next(children))
        monkeypatch.setattr(runtime, "matlab_dll_directories", lambda: ())
        engine = SimpleNamespace(start_matlab=lambda *args, **kwargs: stalled)
        with pytest.raises(TimeoutError, match="MATLAB engine start"):
            runtime.start_matlab(engine, "options")
    else:
        monkeypatch.setattr(runtime, "_matlab_children", lambda: {ours, other})
        mex = tmp_path / "radia_cuda_ams.mexw64"
        mex.touch()
        engine = SimpleNamespace(
            feature=lambda *args, **kwargs: _matlab_future(ours.pid),
            addpath=lambda *args, **kwargs: _matlab_future(),
            which=lambda *args, **kwargs: _matlab_future(str(mex)),
            feval=lambda *args, **kwargs: stalled,
        )
        with pytest.raises(TimeoutError, match="radia_cuda_ams.mexw64 feval"):
            runtime.probe_ams_mex(engine, mex)
    assert killed == [123]
    assert "TIMEOUT" in capsys.readouterr().out
