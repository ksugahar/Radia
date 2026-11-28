"""Build the Cubit ``export vol`` command plugin (cubit_to_vol_export.dll) and _cubit_to_vol_mesh.

The plugin is not part of the radia-mcp wheel.  Build it on the licensed Cubit
machine from an MSVC 14.44 x64 developer environment (``vcvars64.bat
-vcvars_ver=14.44``) with CMake and Ninja on PATH:

    python build.py --cubit-sdk "<Cubit>/cmake" --netgen-dir "<venv>/Lib/site-packages/netgen/cmake" --out <plugin_dir>

``--netgen-dir`` must point at the fork Netgen whose nglib the worker Python
loads (``import ngsolve``) at run time; the plugin binds to that same DLL.
The built ``cubit_to_vol_export.dll`` is copied to ``--out``; pass that directory
as ``plugin_dir`` (or ``RADIA_CUBIT_PLUGIN_DIR``) to ``cubit_export_vol``.
The in-memory mesh module ``_cubit_to_vol_mesh`` (:mod:`radia_mcp.cubit.vol_mesh`)
is built for the Python that runs this script, which must therefore be the
worker Python, and is copied next to the plugin.
"""
from __future__ import annotations

import argparse
import importlib.machinery
import os
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DLL = "cubit_to_vol_export.dll"
MODULE = "_cubit_to_vol_mesh"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cubit-sdk", required=True, help="Cubit CMake package directory (holds CubitConfig.cmake)")
    parser.add_argument("--netgen-dir", required=True, help="Netgen CMake package directory of the fork build")
    parser.add_argument("--out", required=True, help="directory that receives cubit_to_vol_export.dll")
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
                 f"-DCubit_DIR={Path(args.cubit_sdk).as_posix()}", f"-DNetgen_DIR={Path(args.netgen_dir).as_posix()}",
                 f"-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}"]
    subprocess.run(configure, check=True)
    subprocess.run(["cmake", "--build", str(build), "-j", str(args.jobs)], check=True)
    built = build / DLL
    # The module this interpreter (the one CMake built for) imports; a module of
    # another Python left in a reused build or output directory is never taken.
    module = build / (MODULE + importlib.machinery.EXTENSION_SUFFIXES[0])
    for artifact in (built, module):
        if not artifact.is_file():
            raise SystemExit(f"build finished without {artifact}")
    for stale in out.glob(MODULE + "*"):
        if stale.name != module.name and stale.suffix in (".pyd", ".so"):
            stale.unlink()
            print(f"removed {stale} (built for another Python)")
    shutil.copyfile(built, out / DLL)
    print(out / DLL)
    shutil.copyfile(module, out / module.name)
    print(out / module.name)
    return 0


if __name__ == "__main__":
    sys.exit(main())
