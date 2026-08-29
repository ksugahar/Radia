"""Stage the canonical monorepo MATLAB sources into the radia-optuna wheel."""

from __future__ import annotations

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel
from setuptools.command.build_py import build_py as _build_py

PACKAGE_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACKAGE_ROOT.parents[1]
MATLAB_ROOT = REPO_ROOT / "matlab"
OPTUNA_SOURCE = MATLAB_ROOT / "+radia" / "+optuna"
MEX_SOURCE = MATLAB_ROOT / "optuna_mex.mexw64"
CONTRACT_SOURCES = (
    MATLAB_ROOT / "optuna_upstream_compatibility.json",
    MATLAB_ROOT / "optuna49_api_coverage.json",
)
SIMULINK_SOURCES = (
    MATLAB_ROOT / "+radia" / "+simulink" / "buildOptunaBlock.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "buildOptunaTeachingModel.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "optunaSFunction.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "optunaRuntimeStore.m",
    MATLAB_ROOT / "+radia" / "+simulink" / "addOptunaMonitor.m",
    MATLAB_ROOT / "radia_optuna_sfun.m",
    MATLAB_ROOT / "radia_optuna_teaching.slx",
)
DOC_SOURCES = (PACKAGE_ROOT / "OPTUNA_SIMULINK_LAB.md",)


class build_py(_build_py):
    """Copy one canonical source tree; do not maintain a package-local fork."""

    def run(self) -> None:
        missing = [
            path for path in (
                OPTUNA_SOURCE,
                MEX_SOURCE,
                *CONTRACT_SOURCES,
                *SIMULINK_SOURCES,
                *DOC_SOURCES,
            )
            if not path.exists()
        ]
        if missing:
            rendered = "\n".join(f"  - {path}" for path in missing)
            raise RuntimeError(
                "radia-optuna cannot build a partial distribution. Missing:\n"
                f"{rendered}\nRun Build.ps1 -OptunaMexOnly first."
            )

        super().run()
        package_dir = Path(self.build_lib) / "radia_optuna"
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
        shutil.copy2(REPO_ROOT / "LICENSE", matlab_dir / "LICENSE")
        shutil.copy2(
            PACKAGE_ROOT / "THIRD_PARTY_NOTICES.md",
            matlab_dir / "THIRD_PARTY_NOTICES.md",
        )


class bdist_wheel(_bdist_wheel):
    """Tag the MATLAB MEX payload as Windows x64 at build time."""

    def get_tag(self) -> tuple[str, str, str]:
        return "py3", "none", "win_amd64"


setup(cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel})
