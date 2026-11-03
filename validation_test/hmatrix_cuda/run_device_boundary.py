"""Same-kernel pybind11/MEX bitwise parity and the MATLAB device boundary test on one exact build.

Small fixed cases with no skip path; this is not an accuracy or timing measurement. Expects the
`_radia_pybind` and `_hdiv_cuda` builds with their receipts in src/radia, and the MEX build that links
the same kernel object. Writes result.json in --out, holding hashes and pass flags only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import sysconfig

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(name: str, argv: list[str], work: Path, env: dict) -> None:
    log = work / (name + ".log")
    with log.open("wb") as stream:
        code = subprocess.run(argv, cwd=work, env=env, stdout=stream, stderr=subprocess.STDOUT, check=False).returncode
    if code != 0:
        tail = log.read_bytes()[-2000:].decode("utf-8", "replace")
        raise RuntimeError(f"{name} failed with status {code}:\n{tail}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mex-build", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--matlab", default="matlab")
    parser.add_argument("--binary-source-sha", help="Build commit when it differs from HEAD; the native closure must then be checked by the caller")
    args = parser.parse_args()
    work, mex_dir = args.out.resolve(), args.mex_build.resolve()
    work.mkdir(parents=True)
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    binary_sha = args.binary_source_sha or head
    pyd = ROOT / "src/radia" / ("_hdiv_cuda" + sysconfig.get_config_var("EXT_SUFFIX"))
    mex_receipt = json.loads((mex_dir / "radia_hdiv_device.mexw64.build.json").read_bytes())
    pyd_receipt = json.loads(Path(str(pyd) + ".build.json").read_bytes())
    if not mex_receipt["source_commit"] == pyd_receipt["source_commit"] == binary_sha:
        raise ValueError(f"Build source differs: MEX {mex_receipt['source_commit']}, pybind {pyd_receipt['source_commit']}, expected {binary_sha}")
    if mex_receipt["python_binary_sha256"] != sha256(pyd):
        raise ValueError(f"MEX receipt names pybind {mex_receipt['python_binary_sha256']}, staged {sha256(pyd)}")
    compiler = next((pyd_receipt[key] for key in ("cuda_compiler", "nvcc_version", "cuda_version", "nvcc") if pyd_receipt.get(key)), None)
    if not compiler:
        raise ValueError(f"No CUDA compiler identity in the pybind receipt: {sorted(pyd_receipt)}")
    common = dict(source_sha=binary_sha, kernel_object_sha256=mex_receipt["cuda_object"]["sha256"], cuda_compiler=compiler,
                  python_binary_sha256=mex_receipt["python_binary_sha256"], mex_binary_sha256=mex_receipt["binary_sha256"])
    (work / "common.json").write_text(json.dumps(common, indent=1), encoding="utf-8")
    env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(ROOT / "src"), str(ROOT / "packages/radia-simulink/src")]),
               PYTHONNOUSERSITE="1", PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1")
    python, parity = sys.executable, str(HERE / "hdiv_device_parity.py")

    def batch(mode: str) -> list[str]:
        call = f"runDeviceBoundary('{mode}','{ROOT.as_posix()}','{mex_dir.as_posix()}','{work.as_posix()}','{Path(python).as_posix()}')"
        return [args.matlab, "-batch", f"addpath('{HERE.as_posix()}');{call}"]

    run("setup-read", batch("setup-read"), work, env)
    run("parity-input", [python, "-B", parity, "make-input", "--out", str(work / "input.h5")], work, env)
    run("parity-pybind", [python, "-B", parity, "pybind", "--input", str(work / "input.h5"), "--out", str(work / "pybind.h5"),
                          "--common-receipt", str(work / "common.json")], work, env)
    run("parity-mex", batch("parity"), work, env)
    run("parity-compare", [python, "-B", str(HERE / "check_pybind_mex_bitwise.py"), "--python-record", str(work / "pybind.h5"),
                           "--mex-record", str(work / "mex.h5"), "--report", str(work / "bitwise.json")], work, env)
    run("device-test", batch("device-test"), work, env)
    bitwise = json.loads((work / "bitwise.json").read_bytes())
    device = json.loads((work / "device-test.json").read_bytes())
    setup_read = json.loads((work / "setup-read.json").read_bytes())
    passed = bitwise["passed"] is True and device["accepted"] is True and setup_read["accepted"] is True
    result = dict(schema="radia.hdiv-device-boundary.v1", passed=passed, source_sha=head, binary_source_sha=binary_sha,
                  python_binary_sha256=common["python_binary_sha256"], mex_binary_sha256=common["mex_binary_sha256"],
                  kernel_object_sha256=common["kernel_object_sha256"], parity_case=bitwise["case"], parity_outputs=bitwise["outputs"],
                  device_test_scope=device["scope"], device_final_live_count=device["final_live_count"],
                  setup_read_cases=setup_read["cases"], scope="Frontend parity and boundary guards; not accuracy or timing")
    (work / "result.json").write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    if not passed:
        raise RuntimeError(f"Device boundary not accepted: {result}")
    print(json.dumps(dict(passed=passed, parity_outputs=bitwise["outputs"], setup_read_cases=setup_read["cases"])))


if __name__ == "__main__":
    main()
