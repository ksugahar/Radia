"""Build the Cubit ``export netgen`` command plugin (cubit_mesh_export.dll).

The plugin is not part of the radia-mcp wheel.  Build it on the licensed Cubit
machine from an MSVC 14.44 x64 developer environment (``vcvars64.bat
-vcvars_ver=14.44``) with CMake and Ninja on PATH:

    python build.py --cubit-sdk "<Cubit>/cmake" --netgen-dir "<venv>/Lib/site-packages/netgen/cmake" --out <plugin_dir>

``--netgen-dir`` must point at the fork Netgen whose nglib the worker Python
loads (``import ngsolve``) at run time; the plugin binds to that same DLL.
The built ``cubit_mesh_export.dll`` is copied to ``--out``; pass that directory
as ``plugin_dir`` (or ``RADIA_CUBIT_PLUGIN_DIR``) to ``cubit_export_vol``.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DLL = "cubit_mesh_export.dll"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cubit-sdk", required=True, help="Cubit CMake package directory (holds CubitConfig.cmake)")
    parser.add_argument("--netgen-dir", required=True, help="Netgen CMake package directory of the fork build")
    parser.add_argument("--out", required=True, help="directory that receives cubit_mesh_export.dll")
    parser.add_argument("--build-dir", default=None, help="CMake build directory (default <out>/build)")
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args(argv)
    out = Path(args.out).resolve()
    build = Path(args.build_dir).resolve() if args.build_dir else out / "build"
    for name, path in (("--cubit-sdk", args.cubit_sdk), ("--netgen-dir", args.netgen_dir)):
        if not Path(path).is_dir():
            raise SystemExit(f"{name} is not a directory: {path}")
    if not os.environ.get("VCToolsVersion", "").startswith("14.44"):
        raise SystemExit("run from an MSVC 14.44 x64 developer environment (vcvars64.bat -vcvars_ver=14.44)")
    out.mkdir(parents=True, exist_ok=True)
    configure = ["cmake", "-S", str(HERE), "-B", str(build), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
                 f"-DCubit_DIR={Path(args.cubit_sdk).as_posix()}", f"-DNetgen_DIR={Path(args.netgen_dir).as_posix()}"]
    subprocess.run(configure, check=True)
    subprocess.run(["cmake", "--build", str(build), "-j", str(args.jobs)], check=True)
    built = build / DLL
    if not built.is_file():
        raise SystemExit(f"build finished without {built}")
    shutil.copyfile(built, out / DLL)
    print(out / DLL)
    return 0


if __name__ == "__main__":
    sys.exit(main())
