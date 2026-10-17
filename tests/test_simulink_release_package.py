from pathlib import Path
import io
import zipfile

import pytest


ROOT = Path(__file__).resolve().parents[1]








import importlib.util as _scratch_import
_scratch_spec = _scratch_import.spec_from_file_location(
    "repository_contract_scratch", Path(__file__).resolve().parents[1] / "tools/repository_contract_scratch.py")
_scratch_module = _scratch_import.module_from_spec(_scratch_spec)
_scratch_spec.loader.exec_module(_scratch_module)
isolated_engine_contract_scratch = _scratch_module.isolated_controller_scratch


def test_mex_build_verifier_rejects_missing_and_changed_binaries(tmp_path):
    verifier = load_module(
        "verify_radia_mex_build_test",
        ROOT / "packages/radia-gpu/matlab" / "verify_radia_mex_provenance.py",
    )
    mex = tmp_path / "radia_mex.mexw64"
    mex.write_bytes(fake_x64_pe())
    with pytest.raises(RuntimeError, match="provenance is missing"):
        verifier.verify(mex)
    write_test_mex_provenance(tmp_path, "a" * 40)
    assert verifier.verify(mex) == "a" * 40
    mex.write_bytes(b"changed")
    with pytest.raises(RuntimeError, match="does not match"):
        verifier.verify(mex)










def test_full_library_verifier_gates_optional_optimization_toolbox():
    verifier = (
        ROOT / "packages/radia-gpu/matlab" / "verify_radia_simulink_release.m"
    ).read_text(encoding="utf-8")
    assert "electromagnetTopologyExecuted = hasOptimizationToolbox();" in verifier
    assert "electromagnetStatus(end) ~= -1" in verifier
    assert '"electromagnet_topology_executed"' in verifier
    assert 'electromagnetCheck = "dependency-gated";' in verifier
    assert verifier.count(
        'max(abs(maglevForce(:,3,:)), [], "all")'
    ) == 2
    assert '"optuna_mex." + mexext' in verifier
    assert '"api.commands"' in verifier
    assert "optunaInfo.command_count ~= numel(optunaCommands)" in verifier


def test_matlab_smoke_decodes_utf8_without_cp932(monkeypatch, tmp_path):
    verify_module = load_module(
        "verify_simulink_release_encoding",
        ROOT / "tools" / "verify_simulink_release.py",
    )
    archive = tmp_path / "historical.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("README.txt", "Historical archive worker fixture")
    monkeypatch.setattr(verify_module, "verify_archive", lambda _: {"version": "5.2.1"})
    matlab = tmp_path / "matlab.exe"
    matlab.write_bytes(b"MZ")

    class Result:
        returncode = 0
        stdout = "\u691c\u8a3c\u5b8c\u4e86 RADIA_IH_RELEASE_OK\n".encode("utf-8")
        stderr = b""

    def run_worker(command, **kwargs):
        assert command[0] == verify_module.sys.executable
        assert "--engine-worker" in command
        assert "-batch" not in command
        return Result()

    monkeypatch.setattr(verify_module, "_run_matlab_process", run_worker)
    monkeypatch.setattr(verify_module, "_require_candidate_interpreter", lambda manifest: None)
    output = verify_module.run_matlab_smoke(archive, matlab)
    assert "RADIA_IH_RELEASE_OK" in output
    assert verify_module._console_safe("bad:\ufffd", "cp932") == "bad:\\ufffd"


@pytest.mark.parametrize("installed,passes", [("5.2.1", True), ("5.1.0", False), (None, False)])
def test_matlab_smoke_requires_the_package_release_in_its_interpreter(
        monkeypatch, installed, passes):
    # MATLAB runs Radia through RADIA_PYTHON_EXECUTABLE = this interpreter; a
    # host default still on an older release must fail, not pass silently.
    import importlib.metadata
    module = load_module("verify_candidate_interpreter", ROOT / "tools" / "verify_simulink_release.py")

    def version(name):
        assert name == "radia"
        if installed is None:
            raise importlib.metadata.PackageNotFoundError(name)
        return installed

    monkeypatch.setattr(importlib.metadata, "version", version)
    manifest = {"version": "5.2.1", "radia_version": "5.2.1"}
    if passes:
        module._require_candidate_interpreter(manifest)
    else:
        with pytest.raises(RuntimeError, match="candidate interpreter"):
            module._require_candidate_interpreter(manifest)


def test_matlab_timeout_stops_only_owned_windows_tree(monkeypatch):
    module = load_module("verify_simulink_timeout", ROOT / "tools" / "verify_simulink_release.py")
    calls = []

    class Process:
        pid = 43210
        waits = 0

        def wait(self, timeout):
            self.waits += 1
            if self.waits == 1:
                raise module.subprocess.TimeoutExpired("matlab", timeout)
            return 1

    def popen(command, stdout, stderr):
        stdout.write(b"startup diagnostics")
        assert stderr == module.subprocess.STDOUT
        return Process()

    monkeypatch.setattr(module.sys, "platform", "win32")
    monkeypatch.setattr(module.subprocess, "Popen", popen)
    monkeypatch.setattr(module.subprocess, "run", lambda command, **kw: calls.append(command))
    with pytest.raises(RuntimeError, match="startup diagnostics"):
        module._run_matlab_process(["matlab.exe", "-batch", "test"], 1)
    assert calls == [["taskkill", "/PID", "43210", "/T", "/F"]]


def test_matlab_process_keeps_success_output(monkeypatch):
    module = load_module("verify_simulink_process", ROOT / "tools" / "verify_simulink_release.py")

    class Process:
        def wait(self, timeout):
            return 0

    def popen(command, stdout, stderr):
        stdout.write(b"RADIA_SIMULINK_RELEASE_OK")
        return Process()

    monkeypatch.setattr(module.subprocess, "Popen", popen)
    result = module._run_matlab_process(["matlab.exe"], 1)
    assert result.returncode == 0
    assert result.stdout == b"RADIA_SIMULINK_RELEASE_OK"


@pytest.mark.parametrize("failure", [False, True])
def test_engine_worker_quits_its_session_on_success_and_failure(monkeypatch, tmp_path, failure):
    import types
    import sys
    module = load_module("verify_engine_worker", ROOT / "tools" / "verify_simulink_release.py")
    calls = []

    class Engine:
        def matlabroot(self):
            return str(tmp_path)

        def eval(self, expression, **kwargs):
            assert kwargs["nargout"] == 0
            kwargs["stdout"].write("engine output")
            calls.append(expression)
            if failure:
                raise RuntimeError("simulation failed")

        def setenv(self, name, value, **kwargs):
            calls.append((name, value))

        def quit(self):
            calls.append("quit")

    api = types.ModuleType("matlab.engine")
    api.start_matlab = lambda options: Engine()
    api.find_matlab = lambda: ()
    monkeypatch.setattr(module, "_matlab_process_ids", lambda: set())
    parent = types.ModuleType("matlab")
    parent.engine = api
    monkeypatch.setitem(sys.modules, "matlab", parent)
    monkeypatch.setitem(sys.modules, "matlab.engine", api)
    if failure:
        with pytest.raises(RuntimeError, match="simulation failed"):
            module._engine_worker(str(tmp_path), "verify()")
    else:
        assert module._engine_worker(str(tmp_path), "verify()") == 0
    assert calls == [("RADIA_PYTHON_EXECUTABLE", sys.executable), "verify()", "quit"]


@pytest.mark.parametrize("receipt", [None, {"pid": 124, "restored": True}, {"pid": 123, "restored": False}, {"pid": 123, "restored": True}])
def test_com_worker_requires_completed_restored_receipt(monkeypatch, tmp_path, receipt):
    import json
    import re
    import sys
    import types
    from pathlib import Path
    module = load_module("verify_com_worker", ROOT / "tools/verify_simulink_release.py")
    monkeypatch.setattr(module, "_matlab_process_ids", lambda: {123})
    factory = module.tempfile.TemporaryDirectory
    monkeypatch.setattr(module.tempfile, "TemporaryDirectory", lambda **kw: factory(dir=tmp_path))

    class Application:
        def Execute(self, command):
            script = Path(command.removeprefix("run('").removesuffix("');"))
            source = script.read_text(encoding="utf-8")
            destination = Path(re.search(r"fid=fopen\('([^']+)'", source)[1])
            if receipt is not None:
                destination.write_text(json.dumps(receipt), encoding="utf-8")
            return "MATLAB reported an error" if receipt is None else "MATLAB completed"

    client = types.ModuleType("win32com.client")
    def attach(progid):
        assert progid == "Matlab.Desktop.Application"
        return Application()
    client.GetActiveObject = attach
    parent = types.ModuleType("win32com"); parent.client = client
    monkeypatch.setitem(sys.modules, "win32com", parent)
    monkeypatch.setitem(sys.modules, "win32com.client", client)
    if receipt == {"pid": 123, "restored": True}:
        assert module._engine_worker(str(tmp_path), "verify()", "com:123") == 0
    else:
        with pytest.raises(RuntimeError, match="COM"):
            module._engine_worker(str(tmp_path), "verify()", "com:123")
    assert not list(tmp_path.iterdir())


def test_com_worker_refuses_absent_pid_without_connecting(monkeypatch):
    module = load_module("verify_com_missing_pid", ROOT / "tools/verify_simulink_release.py")
    monkeypatch.setattr(module, "_matlab_process_ids", lambda: set())
    with pytest.raises(RuntimeError, match="not running"):
        module._engine_worker("unused", "verify()", "com:123")
    with pytest.raises(ValueError, match="positive PID"):
        module._engine_worker("unused", "verify()", "com:0")


def test_unshared_existing_matlab_prevents_substitute_start(monkeypatch, tmp_path):
    import types
    import sys
    module = load_module("verify_existing_matlab", ROOT / "tools" / "verify_simulink_release.py")
    api = types.ModuleType("matlab.engine")
    api.find_matlab = lambda: ()
    api.start_matlab = lambda *_: pytest.fail("must not start a substitute MATLAB")
    parent = types.ModuleType("matlab")
    parent.engine = api
    monkeypatch.setitem(sys.modules, "matlab", parent)
    monkeypatch.setitem(sys.modules, "matlab.engine", api)
    monkeypatch.setattr(module, "_matlab_process_ids", lambda: {123})
    with pytest.raises(RuntimeError, match="No new MATLAB"):
        module._engine_worker(str(tmp_path), "verify()")


@pytest.mark.parametrize("failure,blocked", [
    (False, None), (True, None), (False, "pid"),
    (False, "model"), (False, "mex"),
])
def test_named_borrowed_engine_is_preserved_and_environment_restored(monkeypatch, tmp_path, failure, blocked):
    import types
    import sys
    module = load_module("verify_borrowed_matlab", ROOT / "tools" / "verify_simulink_release.py")
    calls = []

    class Engine:
        def feature(self, name):
            assert name == "getpid"
            return 999 if blocked == "pid" else 123
        def find_system(self, *args): return ["user_model"] if blocked == "model" else []
        def inmem(self, **kwargs): return [], ["radia_mex"] if blocked == "mex" else []
        def path(self, *args, **kwargs):
            if args: calls.append(("path", args[0]))
            return "original path"
        def pwd(self): return "original folder"
        def cd(self, value, **kwargs): calls.append(("pwd", value))
        def getenv(self, name): return "original " + name
        def setenv(self, name, value, **kwargs): calls.append((name, value))
        def matlabroot(self): return str(tmp_path)
        def eval(self, expression, **kwargs):
            calls.append(expression)
            if expression == "verify()" and failure:
                raise RuntimeError("simulation failed")
        def quit(self): pytest.fail("borrowed MATLAB must remain alive")

    api = types.ModuleType("matlab.engine")
    api.find_matlab = lambda: ("owned_by_user",)
    api.start_matlab = lambda *_: pytest.fail("must reuse the named Engine")
    def connect(name):
        assert name == "owned_by_user"
        return Engine()
    api.connect_matlab = connect
    parent = types.ModuleType("matlab")
    parent.engine = api
    monkeypatch.setitem(sys.modules, "matlab", parent)
    monkeypatch.setitem(sys.modules, "matlab.engine", api)
    monkeypatch.setattr(module, "_matlab_process_ids", lambda: {123, 456})
    if blocked:
        with pytest.raises(RuntimeError):
            module._engine_worker(str(tmp_path), "verify()", "owned_by_user")
        assert calls == []
        return
    if failure:
        with pytest.raises(RuntimeError, match="simulation failed"):
            module._engine_worker(str(tmp_path), "verify()", "owned_by_user")
    else:
        assert module._engine_worker(str(tmp_path), "verify()", "owned_by_user") == 0
    assert calls == [("RADIA_PYTHON_EXECUTABLE", sys.executable),
                     "verify()", "clear radia_mex optuna_mex",
                     ("path", "original path"), ("pwd", "original folder"),
                     ("PATH", "original PATH"),
                     ("MKL_THREADING_LAYER", "original MKL_THREADING_LAYER"),
                     ("PYTHONPATH", "original PYTHONPATH"),
                     ("RADIA_PYTHON_EXECUTABLE", "original RADIA_PYTHON_EXECUTABLE")]


def test_engine_scratch_retries_transient_dll_lock(monkeypatch, tmp_path):
    module = load_module("verify_engine_cleanup", ROOT / "tools" / "verify_simulink_release.py")
    original = module.tempfile.TemporaryDirectory.cleanup
    attempts = []

    def cleanup(instance):
        attempts.append(1)
        if len(attempts) == 1:
            raise PermissionError("DLL still unloading")
        return original(instance)

    monkeypatch.setattr(module.tempfile.TemporaryDirectory, "cleanup", cleanup)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)
    with module._EngineScratch(dir=tmp_path) as folder:
        assert Path(folder).is_dir()
    assert attempts == [1, 1]
    assert not Path(folder).exists()


@pytest.mark.parametrize("damaged", ["bad \ufffd text", "bad ??? text"])
def test_slx_text_integrity_rejects_mojibake(damaged):
    verify_module = load_module(
        "verify_simulink_release_mojibake",
        ROOT / "tools" / "verify_simulink_release.py",
    )
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as model:
        model.writestr("simulink/blockdiagram.xml", damaged)
    with pytest.raises(RuntimeError, match="Replacement glyph|question-mark run"):
        verify_module._verify_slx_text_integrity(
            payload.getvalue(), "packages/radia-gpu/matlab/damaged.slx"
        )


def test_slx_text_integrity_accepts_valid_utf8_and_single_question_mark():
    verify_module = load_module(
        "verify_simulink_release_valid_text",
        ROOT / "tools" / "verify_simulink_release.py",
    )
    payload = io.BytesIO()
    with zipfile.ZipFile(payload, "w") as model:
        model.writestr(
            "simulink/blockdiagram.xml",
            '<?xml version="1.0" encoding="utf-8"?><P>Ready?</P>',
        )
    verify_module._verify_slx_text_integrity(
        payload.getvalue(), "packages/radia-gpu/matlab/valid.slx"
    )


def load_module(name, path):
    import importlib.util

    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_test_mex_provenance(mex_dir, source_commit):
    import json
    import hashlib

    mex = mex_dir / "radia_mex.mexw64"
    data = {
        "schema": "radia.native-build-provenance.v1",
        "source_commit": source_commit,
        "source_dirty": False,
        "binary_name": mex.name,
        "binary_bytes": mex.stat().st_size,
        "binary_sha256": hashlib.sha256(mex.read_bytes()).hexdigest(),
    }
    (mex_dir / "radia_mex.mexw64.build.json").write_text(
        json.dumps(data), encoding="utf-8")


def fake_x64_pe():
    payload = bytearray(2048)
    payload[:2] = b"MZ"
    payload[60:64] = (128).to_bytes(4, "little")
    payload[128:132] = b"PE\0\0"
    payload[132:134] = (0x8664).to_bytes(2, "little")
    return bytes(payload)
