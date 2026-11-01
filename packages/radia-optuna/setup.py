"""Stage the canonical monorepo MATLAB sources into the radia-optuna wheel."""

from __future__ import annotations

import json
import shutil
import hashlib
import sysconfig
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.command.build_py import build_py as _build_py

PACKAGE_ROOT = Path(__file__).resolve().parent
MATLAB_ROOT = PACKAGE_ROOT / "matlab"
OPTUNA_SOURCE = MATLAB_ROOT / "+radia" / "+optuna"
EVIDENCE = json.loads(
    (PACKAGE_ROOT / "src/radia_optuna/evidence_manifest.json").read_text(encoding="utf-8")
)
ORACLE_SOURCES = tuple(
    PACKAGE_ROOT / record["path"] for record in EVIDENCE["files"].values()
    if record["path"].startswith("oracles/")
)
MEX_SOURCE = MATLAB_ROOT / "optuna_mex.mexw64"
NATIVE_SOURCE = PACKAGE_ROOT / "src/radia_optuna" / ("_optuna_kernels" + sysconfig.get_config_var("EXT_SUFFIX"))
BUILD_RECEIPT = PACKAGE_ROOT / "src/radia_optuna/native_build.json"
CONTRACT_SOURCES = (
    MATLAB_ROOT / "optuna_upstream_compatibility.json",
    MATLAB_ROOT / "optuna50_api_coverage.json",
)
SIMULINK_SOURCES = (
    MATLAB_ROOT / "+radia" / "+simulink" / "buildOptunaBlock.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "buildOptunaStudyBlock.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "makeMonitorHeaderBusObject.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "makeOptunaMonitorBusObject.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "buildOptunaTeachingModel.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "reviewOptunaStudy.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "applyOptunaTrial.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "optunaSFunction.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "optunaRuntimeStore.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "addOptunaMonitor.m",
    MATLAB_ROOT / "radia_optuna_sfun.m",
    MATLAB_ROOT / "radia_optuna_teaching.slx",
)
DOC_SOURCES = (PACKAGE_ROOT / "OPTUNA_SIMULINK_LAB.md", PACKAGE_ROOT / "NATIVE_BUILD.md")


class build_py(_build_py):
    """Copy one canonical source tree; do not maintain a package-local fork."""

    def run(self) -> None:
        missing = [
            path for path in (
                OPTUNA_SOURCE,
                MEX_SOURCE,
                NATIVE_SOURCE,
                BUILD_RECEIPT,
                *CONTRACT_SOURCES,
                *ORACLE_SOURCES,
                *SIMULINK_SOURCES,
                *DOC_SOURCES,
            )
            if not path.exists()
        ]
        if missing:
            rendered = "\n".join(f"  - {path}" for path in missing)
            raise RuntimeError(
                "radia-optuna cannot build a partial distribution. Missing:\n"
                f"{rendered}\nRun buildOptunaMex in MATLAB from the package source directory first."
            )

        super().run()
        package_dir = Path(self.build_lib) / "radia_optuna"
        oracle_dir = package_dir / "oracles"
        oracle_dir.mkdir(parents=True, exist_ok=True)
        for source in ORACLE_SOURCES:
            shutil.copy2(source, oracle_dir / source.name)
        matlab_dir = package_dir / "matlab"
        staged_optuna = matlab_dir / "+radia" / "+optuna"
        shutil.copytree(OPTUNA_SOURCE, staged_optuna, dirs_exist_ok=True)
        shutil.copy2(MEX_SOURCE, matlab_dir / MEX_SOURCE.name)
        for source in CONTRACT_SOURCES:
            shutil.copy2(source, matlab_dir / source.name)
        for source in SIMULINK_SOURCES:
            target = matlab_dir / source.relative_to(MATLAB_ROOT)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source, target)
        for source in DOC_SOURCES:
            shutil.copy2(source, matlab_dir / source.name)
        shutil.copy2(PACKAGE_ROOT / "MATLAB_README.md", matlab_dir / "README.md")
        shutil.copy2(PACKAGE_ROOT / "LICENSE", matlab_dir / "LICENSE")
        shutil.copy2(
            PACKAGE_ROOT / "THIRD_PARTY_NOTICES.md",
            matlab_dir / "THIRD_PARTY_NOTICES.md",
        )


class build_ext(_build_ext):
    """Stage only the CMake-built exact Python ABI; never rebuild/fallback."""
    def run(self):
        receipt=json.loads(BUILD_RECEIPT.read_text(encoding="utf-8"))
        expected=receipt["python_extension_sha256s"].get(NATIVE_SOURCE.name)
        if not expected or hashlib.sha256(NATIVE_SOURCE.read_bytes()).hexdigest()!=expected:
            raise RuntimeError("Missing/tampered exact-ABI native extension")
        if hashlib.sha256(MEX_SOURCE.read_bytes()).hexdigest()!=receipt["mex_sha256"]:
            raise RuntimeError("MEX does not match the shared-object build receipt")
        target=Path(self.get_ext_fullpath("radia_optuna._optuna_kernels"))
        target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(NATIVE_SOURCE,target)


setup(ext_modules=[Extension("radia_optuna._optuna_kernels", sources=[])],
      cmdclass={"build_py": build_py, "build_ext": build_ext})
