"""Guards: cubic hex export edge orientation and the actual MSH 4.2 handoff.

Assumes: manager-admitted licensed batch execution on the build role, the
candidate command plugin, and an isolated Gmsh 5 environment with install report.
Calls the production C++ exporter, then an independent Gmsh Jacobian evaluation;
this is geometry interoperability validation, not solver accuracy certification.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from cubit_mesh_export._layout import runtime_environment
from cubit_mesh_export._gmsh_runtime import probe_gmsh5
from cubit_mesh_export.mcp._support.mesh_quality import mesh_total_volume


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cubit-exe", type=Path, required=True)
    parser.add_argument("--plugin", type=Path, required=True)
    parser.add_argument("--install-report", type=Path, required=True)
    parser.add_argument("--wheel", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    result = {"passed": False, "state": "running", "scope": "cubic hex geometry export"}
    try:
        exe, plugin, report = (value.resolve(strict=True) for value in
                               (args.cubit_exe, args.plugin, args.install_report))
        os.environ["CUBIT_GMSH_INSTALL_REPORT"] = str(report)
        os.environ["CUBIT_GMSH_WHEEL"] = str(args.wheel.resolve(strict=True))
        result["reference_runtime"] = probe_gmsh5()
        stage = out / "plugins"
        stage.mkdir()
        candidate = stage / "cubit_mesh_export.dll"
        shutil.copy2(plugin, candidate)
        digest = hashlib.sha256(candidate.read_bytes()).hexdigest()
        result["candidate_sha256"] = digest
        script = f'''import cubit,ctypes,json,hashlib,traceback
from pathlib import Path
out=Path({str(out)!r})
result={{'passed':False}}
try:
 def command(text):
  before=cubit.get_error_count();cubit.cmd(text)
  if cubit.get_error_count()!=before:raise RuntimeError(text)
 command('help export gmsh')
 kernel=ctypes.WinDLL('kernel32',use_last_error=True)
 kernel.GetModuleHandleW.argtypes=[ctypes.c_wchar_p];kernel.GetModuleHandleW.restype=ctypes.c_void_p
 kernel.GetModuleFileNameW.argtypes=[ctypes.c_void_p,ctypes.c_wchar_p,ctypes.c_uint32]
 handle=kernel.GetModuleHandleW('cubit_mesh_export.dll')
 if not handle:raise RuntimeError('candidate DLL not loaded')
 buffer=ctypes.create_unicode_buffer(32768)
 if not kernel.GetModuleFileNameW(handle,buffer,len(buffer)):raise RuntimeError('DLL identity unavailable')
 actual=Path(buffer.value).resolve()
 if actual!=(out/'plugins/cubit_mesh_export.dll').resolve():raise RuntimeError('wrong loaded DLL')
 loaded=hashlib.sha256(actual.read_bytes()).hexdigest()
 if loaded!={digest!r}:raise RuntimeError('loaded DLL bytes differ')
 for text in ['reset','create brick x 1 y 1 z 1','volume 1 scheme map',
              'curve all interval 2','mesh volume 1','block 1 add hex all',
              'block 1 name "cubic_hex"']:
  command(text)
 count=len(cubit.get_block_hexes(1))
 if not count:raise RuntimeError('empty input')
 command('export gmsh "'+(out/'hex-p3.msh').as_posix()+'" order 3 overwrite')
 result.update(passed=True,loaded_sha256=loaded,primary_elements=count)
except Exception:result['error']=traceback.format_exc()
(out/'cubit-result.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
'''
        batch = out / "batch.py"
        batch.write_text(script, encoding="utf-8")
        env = runtime_environment(exe.parent)
        env["CUBIT_PLUGIN_DIR"] = str(stage)
        for key in list(env):
            if key.startswith("PYTHON"):
                env.pop(key)
        env.update(OMP_NUM_THREADS="2", NGS_NUM_THREADS="2",
                   MKL_NUM_THREADS="2", OPENBLAS_NUM_THREADS="2")
        command = [str(exe), "-nographics", "-batch", "-nojournal", "-noinitfile", str(batch)]
        with (out / "cubit.log").open("wb") as log:
            process = subprocess.run(command, cwd=out, env=env, stdout=log,
                                     stderr=subprocess.STDOUT, timeout=300, check=False)
        receipt = json.loads((out / "cubit-result.json").read_text(encoding="utf-8"))
        result.update(cubit=receipt, process_exit_code=process.returncode)
        if process.returncode or receipt.get("passed") is not True:
            raise RuntimeError("licensed exporter execution failed")
        os.environ["CUBIT_GMSH_INSTALL_REPORT"] = str(report)
        quality = mesh_total_volume(out / "hex-p3.msh", quadrature="Gauss4")
        result["gmsh"] = quality
        if not quality.get("ok") or quality.get("min_jacobian_det", 0) <= 0:
            raise RuntimeError("exported cubic hex has invalid Gmsh Jacobians")
        result["mesh_sha256"] = hashlib.sha256((out / "hex-p3.msh").read_bytes()).hexdigest()
        result["passed"] = True
    except Exception as exc:
        result["error"] = f"{type(exc).__name__}: {exc}"
    finally:
        result["state"] = "finished"
        (out / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
