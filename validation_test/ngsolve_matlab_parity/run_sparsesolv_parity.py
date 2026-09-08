"""Run the focused MEX tests through a private MATLAB Engine; retain JSON evidence."""
import argparse
import datetime
import hashlib
import json
import platform
import os
import sys
from pathlib import Path
import subprocess
import tempfile

import matlab.engine


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--timeout", type=int, default=240)
    args = parser.parse_args()
    if not args.worker:
        child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()),
                                  "--worker", "--output", str(args.output.resolve())])
        try:
            raise SystemExit(child.wait(timeout=args.timeout))
        except subprocess.TimeoutExpired:
            # Kill only this runner's owned Engine process tree.
            subprocess.run(["taskkill", "/PID", str(child.pid), "/T", "/F"], check=True)
            child.wait(timeout=20)
            args.output.parent.mkdir(parents=True, exist_ok=True)
            args.output.write_text(json.dumps(dict(passed=False, error="owned Engine worker timeout")))
            raise SystemExit(124)
    root = Path(__file__).resolve().parents[2]
    if not (root/"src/radia/_radia_pybind.pyd").is_file():
        raise RuntimeError("Build -RadiaOnly in this checkout before parity; do not borrow an editable binary")
    os.environ["PYTHONPATH"] = str(root/"src") + os.pathsep + os.environ.get("PYTHONPATH", "")
    mex = root / "matlab/radia_mex.mexw64"
    record = dict(schema="radia.sparsesolv-matlab-parity.v1",
                  utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  host=platform.node(), python=platform.python_version(),
                  mex=str(mex), mex_sha256=hashlib.sha256(mex.read_bytes()).hexdigest(),
                  source_state="Working tree snapshot; source_hashes identify tested files, not the base commit alone",
                  source_base_commit=subprocess.check_output(["git", "-c", f"safe.directory={root.as_posix()}",
                      "rev-parse", "HEAD"], cwd=root, text=True).strip(),
                  source_hashes={str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                      for p in [root/"src/matlab/radia_mex.cpp", root/"tests/matlab/test_sparsesolv_mex.m",
                                root/"tests/matlab/sparsesolv_python_reference.py",
                                *sorted((root/"src/ext/sparsesolv/include").rglob("*.hpp"))]},
                  passed=False)
    eng = None
    scratch = tempfile.TemporaryDirectory(prefix="sparsesolv-matlab-", dir="C:/temp")
    metrics = Path(scratch.name)/"metrics.json"
    try:
        eng = matlab.engine.start_matlab("-nodesktop -nosplash")
        eng.cd(str(root), nargout=0)
        eng.addpath(str(root/"matlab"), nargout=0)
        eng.setenv("RADIA_SPARSESOLV_METRICS", str(metrics), nargout=0)
        eng.eval("radia.setup(Force=true);", nargout=0)
        record["matlab"] = eng.version()
        record["resolved_mex"] = eng.which("radia_mex")
        if Path(record["resolved_mex"]).resolve() != mex.resolve():
            raise RuntimeError("MATLAB resolved a different MEX")
        eng.workspace["testfile"] = str(root/"tests/matlab/test_sparsesolv_mex.m")
        eng.eval("r = runtests(testfile); disp(table(r));", nargout=0)
        record["tests"] = json.loads(eng.eval("jsonencode(struct('names',{string({r.Name})},'passed',[r.Passed],'failed',[r.Failed],'incomplete',[r.Incomplete],'duration',[r.Duration]))"))
        record["passed"] = bool(eng.eval("~isempty(r) && all([r.Passed])"))
        if metrics.is_file():
            record["numerical"] = json.loads(metrics.read_text())
    except Exception as exc:
        record["error"] = str(exc)
        raise
    finally:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2), encoding="utf-8")
        if eng is not None:
            eng.quit()
        scratch.cleanup()
    if not record["passed"]:
        raise SystemExit("MATLAB sparsesolv parity failed")


if __name__ == "__main__":
    main()
