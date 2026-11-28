"""Build and run the license-free native contracts of the Cubit export source.

Builds ``tests/`` (named periodic matcher, curve identity through the callback
geometry, periodic curved coefficients) and the ``_netgen_curving`` test
module against the Netgen SDK of the running Python, runs them with CTest, then
runs ``validation_test/radia_mcp/test_netgen_curving.py`` in strict mode.  No Cubit is
needed.  Run from an MSVC 14.44 x64 developer environment with CMake and Ninja
on PATH, using the Python whose netgen/ngsolve the tests must load:

    python run_contracts.py --build-dir <empty or new directory>

Every missing build input (compiler, Netgen SDK, test module) is an error.
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PACKAGE = HERE.parents[1]
REPO = PACKAGE.parents[1]
CURVING_TEST = REPO / "validation_test" / "radia_mcp" / "test_netgen_curving.py"


def netgen_cmake_dir() -> Path:
    spec = importlib.util.find_spec("netgen")
    if spec is None or not spec.origin:
        raise SystemExit(f"netgen is not importable from {sys.executable}")
    cmake_dir = Path(spec.origin).resolve().parent / "cmake"
    if not (cmake_dir / "NetgenConfig.cmake").is_file():
        raise SystemExit(f"Netgen SDK is missing: {cmake_dir / 'NetgenConfig.cmake'}")
    return cmake_dir


def run(command: list[str], **kwargs) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True, **kwargs)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--build-dir", required=True, help="CMake build directory")
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args(argv)
    if os.name == "nt" and not os.environ.get("VCToolsVersion", "").startswith("14.44"):
        raise SystemExit("run from an MSVC 14.44 x64 developer environment (vcvars64.bat -vcvars_ver=14.44)")
    build = Path(args.build_dir).resolve()
    run(["cmake", "-S", str(HERE), "-B", str(build), "-G", "Ninja", "-DCMAKE_BUILD_TYPE=Release",
         "-DCUBIT_TO_VOL_EXPORT_PLUGIN=OFF", "-DNETGEN_CURVING_PYMODULE=ON",
         "-DCUBIT_TO_VOL_EXPORT_CONTRACT_TESTS=ON",
         f"-DNetgen_DIR={netgen_cmake_dir().as_posix()}",
         f"-DPython3_EXECUTABLE={Path(sys.executable).as_posix()}"])
    run(["cmake", "--build", str(build), "-j", str(args.jobs)])
    run(["ctest", "--test-dir", str(build), "--output-on-failure", "--no-tests=error"])
    env = dict(os.environ, RADIA_REQUIRE_NETGEN_CURVING="1", RADIA_NETGEN_CURVING_DIR=str(build))
    # validation_test/conftest.py imports the repository's tests package.
    run([sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider", str(CURVING_TEST)],
        env=env, cwd=REPO)
    return 0


if __name__ == "__main__":
    sys.exit(main())
