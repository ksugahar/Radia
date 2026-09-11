"""Measure MATLAB setup overhead in a private Engine, writing JSON evidence.

This probes gateway overhead with api.commands, not numerical solver speed.
--mex-dir may name an existing artifact directory; its hash is recorded and
no installed package or user MATLAB session is modified.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time


def main():
    root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mex-dir", type=Path, default=root / "matlab")
    parser.add_argument("--tests", action="store_true")
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if not args.worker:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:], "--worker"])
        try:
            return child.wait(timeout=240)
        except subprocess.TimeoutExpired:
            # Only the Engine worker owned by this benchmark is terminated.
            subprocess.run(["taskkill", "/PID", str(child.pid), "/T", "/F"], check=True)
            child.wait(timeout=20)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps({"passed": False, "error": "owned worker timeout"}))
            return 124
    import matlab.engine
    mex = args.mex_dir.resolve() / "radia_mex.mexw64"
    sources = ["matlab/+radia/setup.m", "matlab/+radia/+internal/callMex.m",
               "matlab/+radia/+internal/pythonProcessPath.m",
               "matlab/+radia/+simulink/configureFileGeneration.m",
               "tests/matlab/test_mex_runtime_setup.m", "validation_test/matlab_runtime/benchmark_setup.py"]
    record = dict(schema="radia.matlab-setup-overhead.v1", passed=False,
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  host=platform.node(), python=platform.python_version(),
                  scope="api.commands gateway overhead; not numerical solver performance",
                  mex=str(mex), mex_sha256=hashlib.sha256(mex.read_bytes()).hexdigest(),
                  source_hashes={p: hashlib.sha256((root / p).read_bytes()).hexdigest() for p in sources},
                  source_base_commit=subprocess.check_output(
                      ["git", "-c", f"safe.directory={root.as_posix()}", "-C", str(root), "rev-parse", "HEAD"],
                      text=True).strip(),
                  source_state="Source hashes identify the working snapshot, not base commit alone")
    eng = None
    try:
        started = time.perf_counter()
        eng = matlab.engine.start_matlab("-nodesktop -nosplash")
        record["engine_start_seconds"] = time.perf_counter() - started
        eng.cd("C:/temp", nargout=0)
        eng.addpath(str(mex.parent), "-end", nargout=0)
        eng.addpath(str(root / "matlab"), "-begin", nargout=0)
        assert Path(eng.which("radia.setup")).resolve() == (root / "matlab/+radia/setup.m").resolve()
        assert Path(eng.which("radia_mex")).resolve() == mex
        record["matlab"] = eng.version()
        eng.eval("tic; info=radia.setup(Force=true); cold=toc;", nargout=0)
        record["cold_setup_seconds"] = eng.workspace["cold"]
        record["samples_seconds_per_call"] = {}
        for name, expression in [
                ("setup", "radia.setup();"),
                ("runtime_only", "radia.setup(ConfigureSimulinkFileGeneration=false);"),
                ("filegen", "radia.simulink.configureFileGeneration();"),
                ("call_mex", 'commands=radia.internal.callMex("api.commands");'),
                ("direct_mex", "commands=radia_mex('api.commands');")]:
            eng.eval("times=zeros(1,7); for k=1:7; tic; for j=1:100; " + expression +
                     " end; times(k)=toc/100; end;", nargout=0)
            samples = json.loads(eng.eval("jsonencode(times)"))
            record["samples_seconds_per_call"][name] = samples
            print(name, statistics.median(samples), flush=True)
        if args.tests:
            eng.workspace["testfile"] = str(root / "tests/matlab/test_mex_runtime_setup.m")
            eng.eval("r=runtests(testfile); disp(table(r));", nargout=0)
            record["tests"] = json.loads(eng.eval(
                "jsonencode(struct('names',{string({r.Name})},'passed',[r.Passed],'incomplete',[r.Incomplete]))"))
            if not eng.eval("~isempty(r) && all([r.Passed])"):
                raise RuntimeError("MATLAB runtime regression failed")
        record["passed"] = True
    except Exception as exc:
        record["error"] = str(exc)
        raise
    finally:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
        if eng is not None:
            eng.quit()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
