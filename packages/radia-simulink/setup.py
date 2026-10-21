"""Stage canonical, explicit GPU/Simulink/LTspice payload into a Windows wheel."""
from pathlib import Path
import hashlib
import json
import shutil
import subprocess
import os
import importlib.util
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 build hosts
    import tomli as tomllib
from setuptools import setup
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parents[1]
MATLAB = PACKAGE / "matlab"
SIMULINK_FILES = (
    "+radia/+simulink/cudaPCGSFunction.m", "+radia/+simulink/buildCudaPCGBlock.m",
    "+radia/+simulink/buildLTspiceBlock.m", "+radia/+simulink/ltspiceSFunction.m",
    "+radia/+simulink/runLTspice.m", "+radia/+simulink/runLTspiceInterval.m",
    "+radia/+simulink/runLTspiceClosedLoop.m", "+radia/+simulink/+internal/runArtifacts.m",
    "radia_simulink_pcg_sfun.m", "radia_ltspice_sfun.m",
)


def payload_module():
    # PEP 517 backends need not put this setup directory on sys.path.
    spec = importlib.util.spec_from_file_location('radia_simulink_payload', PACKAGE / 'hdiv_payload.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class build_py(_build_py):
    def run(self):
        build = ROOT / "build-cuda-mex"
        receipt = json.loads((build / "build.json").read_text())
        commit = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
        dirty = subprocess.check_output(["git", "-C", str(ROOT), "status", "--porcelain"], text=True).strip()
        binary = build / "radia_cuda_csr_pcg.mexw64"
        source = PACKAGE / "native/radia_cuda_csr_pcg.cu"
        if dirty or receipt["source_dirty"] or receipt["source_commit"] != commit:
            raise RuntimeError("radia-simulink requires a clean exact-commit mexcuda build")
        for path, field in ((binary, "binary_sha256"), (source, "source_sha256")):
            if hashlib.sha256(path.read_bytes()).hexdigest() != receipt[field]:
                raise RuntimeError(f"CUDA build provenance mismatch: {path.name}")
        ic_binary = build / 'radia_cuda_ic.mexw64'
        ic_source = PACKAGE / 'native/radia_cuda_ic.cu'
        ic_hash = hashlib.sha256(ic_binary.read_bytes()).hexdigest()
        ic_source_hash = hashlib.sha256(ic_source.read_bytes()).hexdigest()
        if (receipt.get('binary_sha256_by_file', {}).get(ic_binary.name) != ic_hash
                or receipt.get('source_sha256_by_file', {}).get(ic_source.name) != ic_source_hash):
            raise RuntimeError('IC CUDA build provenance mismatch')
        ams_build = ROOT / "build-ams-cuda"
        ams_receipt = json.loads((ams_build / "build.json").read_text(encoding="utf-8"))
        ams_binary = ams_build / "radia_cuda_ams.mexw64"
        ams_source = PACKAGE / "native/radia_cuda_ams.cu"
        ams_hash = hashlib.sha256(ams_binary.read_bytes()).hexdigest()
        if (ams_receipt["source_dirty"] or ams_receipt["source_commit"] != commit
                or ams_hash != ams_receipt["sha256"]["build-ams-cuda/radia_cuda_ams.mexw64"]
                or hashlib.sha256(ams_source.read_bytes()).hexdigest() !=
                ams_receipt["source_sha256"]["packages/radia-simulink/native/radia_cuda_ams.cu"]):
            raise RuntimeError("AMS CUDA build provenance mismatch")
        super().run()
        package = Path(self.build_lib) / "radia_simulink"
        staged = package / "matlab"
        staged.mkdir(parents=True, exist_ok=True)
        expected = {p.relative_to(PACKAGE / "src/radia_simulink").as_posix() for p in (PACKAGE / "src/radia_simulink").rglob("*") if p.is_file() and "__pycache__" not in p.parts}
        expected.update({"LICENSE", "README.md", "build.json", "ams-build.json", "ic-build.json"})
        hmatrix = os.environ.get('RADIA_HDIV_HMATRIX_BUILD')
        projection = os.environ.get('RADIA_HDIV_PROJECTION_BUILD')
        if bool(hmatrix) != bool(projection):
            raise RuntimeError('Select both HDiv HMatrix and projection build directories')
        if hmatrix:
            expected.update(payload_module().stage(ROOT, package, commit, hmatrix, projection))
        fmm = os.environ.get('RADIA_FMM_BUILD')
        if fmm:
            expected.update(payload_module().stage_fmm(ROOT, package, commit, fmm))
        for group in ("+radia/+gpu", "+radia/+ltspice"):
            for path in (MATLAB / group).rglob("*.m"):
                if path.name == 'LaplaceFMM.m' and not fmm:
                    continue
                if path.name in {"FieldCloud.m", "FieldElements.m", "HMatrixOperator.m", "solveHDivRT0.m", "solveCurvedHDiv.m", "curvedHDivQuadrature.m", "solveHDivElements.m", "affineHDivQuadrature.m"}:
                    continue  # Experimental field/H-matrix modules require separate acceptance.
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
        shutil.copy2(ic_binary, staged / ic_binary.name)
        expected.add('matlab/' + ic_binary.name)
        public_ic = dict(schema='radia.gpu.ic-payload.v1', source_commit=commit, source_dirty=False,
                         source_sha256=ic_source_hash, binary_sha256=ic_hash,
                         gpu_numerical_factorization=True, no_cpu_solver_fallback=True,
                         factorization='FP64 fixed-shift natural-order IC(0)',
                         triangular_backend='cuSPARSE SpSV', cuda_targets=receipt['cuda_targets'],
                         matlab=receipt['matlab'])
        (package/'ic-build.json').write_text(json.dumps(public_ic,indent=2)+'\n')
        shutil.copy2(ams_binary, staged / ams_binary.name)
        expected.add("matlab/" + ams_binary.name)
        for name in ("readAMSExchange.m", "amsPCG.m"):
            target = staged / "+radia/+gpu" / name
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(PACKAGE / "matlab/+radia/+gpu" / name, target)
            expected.add(target.relative_to(package).as_posix())
        public_ams = {key: ams_receipt[key] for key in (
            "schema", "source_commit", "source_dirty", "matlab", "no_cpu_solver_fallback")}
        public_ams.update(binary_sha256=ams_hash,
                          source_sha256=ams_receipt["source_sha256"]["packages/radia-simulink/native/radia_cuda_ams.cu"])
        (package / "ams-build.json").write_text(json.dumps(public_ams, indent=2)+"\n", encoding="utf-8")
        shutil.copy2(PACKAGE / "LICENSE", package / "LICENSE")
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
                        capabilities=["CUDA FP64 CSR Jacobi-PCG", "Experimental CUDA FP64 AMS", "Experimental reusable CUDA FP64 IC(0)-PCG", "Reusable affine linear HCurl matrix and AMS hierarchy", "AMS HDF5 exchange", "Simulink CUDA PCG block", "LTspice sampled-data integration"] + (["Bounded affine/P2 HDiv orders 1/2"] if hmatrix else []) + (["Experimental Laplace FMM and flat P1 boundary maps"] if fmm else []),
                        runtime_requirements=["Windows x64", "MATLAB R2026a", "Parallel Computing Toolbox", "NVIDIA GPU/driver", "Simulink for blocks", "User-installed LTspice for circuit execution"])
        (package / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


class bdist_wheel(_bdist_wheel):
    def get_tag(self):
        return "py3", "none", "win_amd64"


setup(cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel})
