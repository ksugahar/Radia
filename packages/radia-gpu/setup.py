"""Stage canonical, explicit GPU/Simulink/LTspice payload into a Windows wheel."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 build hosts
    import tomli as tomllib
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parents[1]
MATLAB = ROOT / "matlab"
SIMULINK_FILES = (
    "+radia/+simulink/cudaPCGSFunction.m", "+radia/+simulink/buildCudaPCGBlock.m",
    "+radia/+simulink/buildLTspiceBlock.m", "+radia/+simulink/ltspiceSFunction.m",
    "+radia/+simulink/runLTspice.m", "+radia/+simulink/runLTspiceInterval.m",
    "+radia/+simulink/runLTspiceClosedLoop.m", "+radia/+simulink/+internal/runArtifacts.m",
    "radia_gpu_pcg_sfun.m", "radia_ltspice_sfun.m",
)


class build_py(_build_py):
    def run(self):
        build = ROOT / "build-cuda-mex"
        receipt = json.loads((build / "build.json").read_text())
        commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
        dirty = subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip()
        binary = build / "radia_cuda_csr_pcg.mexw64"
        source = ROOT / "src/matlab/cuda/radia_cuda_csr_pcg.cu"
        if dirty or receipt["source_dirty"] or receipt["source_commit"] != commit:
            raise RuntimeError("radia-gpu requires a clean exact-commit mexcuda build")
        for path, field in ((binary, "binary_sha256"), (source, "source_sha256")):
            if hashlib.sha256(path.read_bytes()).hexdigest() != receipt[field]:
                raise RuntimeError(f"CUDA build provenance mismatch: {path.name}")
        super().run()
        package = Path(self.build_lib) / "radia_gpu"
        staged = package / "matlab"
        staged.mkdir(parents=True, exist_ok=True)
        expected = {"__init__.py", "cli.py", "LICENSE", "README.md", "build.json"}
        for group in ("+radia/+gpu", "+radia/+ltspice"):
            for path in (MATLAB / group).rglob("*.m"):
                if path.name in ("runRLCResonanceOptimizationDemo.m", "rlcResonanceObjective.m"):
                    continue  # Independent Optuna-dependent examples are not a runtime dependency.
                target = staged / path.relative_to(MATLAB)
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(path, target)
                expected.add(target.relative_to(package).as_posix())
        for relative in SIMULINK_FILES:
            target = staged / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(MATLAB / relative, target)
            expected.add(target.relative_to(package).as_posix())
        shutil.copy2(binary, staged / binary.name)
        expected.add("matlab/" + binary.name)
        shutil.copy2(ROOT / "LICENSE", package / "LICENSE")
        shutil.copy2(PACKAGE / "README.md", package / "README.md")
        # Public provenance deliberately omits private host/interpreter/install paths.
        public_build = {key: receipt[key] for key in (
            "schema", "source_commit", "source_dirty", "source_sha256", "binary_sha256",
            "matlab", "release", "no_cpu_solver_fallback", "cuda_targets")}
        (package / "build.json").write_text(json.dumps(public_build, indent=2)+"\n")
        actual = {p.relative_to(package).as_posix() for p in package.rglob("*")
                  if p.is_file() and p.name != "manifest.json"}
        if actual != expected:
            raise RuntimeError(f"Unexpected staging inventory; use a fresh build directory: {sorted(actual ^ expected)}")
        payload = {p.relative_to(package).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in package.rglob("*") if p.is_file() and p.name != "manifest.json"}
        metadata = tomllib.loads((PACKAGE / "pyproject.toml").read_text())
        manifest = dict(schema="radia.gpu.distribution.v1", version=metadata["project"]["version"],
                        source_commit=commit, files=payload,
                        capabilities=["CUDA FP64 CSR Jacobi-PCG", "Simulink CUDA PCG block", "LTspice sampled-data integration"],
                        runtime_requirements=["Windows x64", "MATLAB R2026a", "Parallel Computing Toolbox", "NVIDIA GPU/driver", "Simulink for blocks", "User-installed LTspice for circuit execution"])
        (package / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


class bdist_wheel(_bdist_wheel):
    def get_tag(self):
        return "py3", "none", "win_amd64"


setup(cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel})
