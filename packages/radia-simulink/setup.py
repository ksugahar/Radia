"""Stage canonical, explicit GPU/Simulink/LTspice payload into a Windows wheel."""
from pathlib import Path
import hashlib
import json
import math
import shutil
import subprocess
import os
import sys
import sysconfig
import importlib.util
try:
    import tomllib
except ModuleNotFoundError:  # Python 3.10 build hosts
    import tomli as tomllib
from setuptools import setup, Distribution
from setuptools.command.build_py import build_py as _build_py
from setuptools.command.bdist_wheel import bdist_wheel as _bdist_wheel

PACKAGE = Path(__file__).resolve().parent
ROOT = PACKAGE.parents[1]
MATLAB = PACKAGE / "matlab"
SIMULINK_FILES = (
    "+radia/csrPCG.m",
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


def payload_module_native():
    spec = importlib.util.spec_from_file_location(
        'radia_simulink_native_payload', PACKAGE / 'native_payload.py')
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def payload_module_resident():
    spec = importlib.util.spec_from_file_location(
        'radia_simulink_resident_payload', PACKAGE / 'resident_amg_payload.py')
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
        python_binding = receipt.get('python_binding')
        if python_binding:
            name = '_csr_cuda'+sysconfig.get_config_var('EXT_SUFFIX')
            native = build / name
            if (python_binding.get('binary') != name
                    or hashlib.sha256(native.read_bytes()).hexdigest() != python_binding.get('binary_sha256')):
                raise RuntimeError('Python CUDA binding ABI or provenance mismatch')
            runtime_path=build/'python-result.json'
            runtime=json.loads(runtime_path.read_text())
            gpu_cases=[case for case in runtime.get('cases',[]) if case.get('backend')=='gpu']
            if (runtime.get('accepted') is not True or runtime.get('source_commit')!=commit
                    or runtime.get('binary_sha256')!=python_binding['binary_sha256'] or not gpu_cases
                    or not runtime.get('runtime_identity',{}).get('sha256')):
                raise RuntimeError('Exact-source Python GPU runtime acceptance required')
            for case in gpu_cases:
                for key,limit in (('relative_l2',1e-8),('true_relative_residual',1e-6)):
                    value=case.get(key)
                    if isinstance(value,bool) or not isinstance(value,(int,float)) or not math.isfinite(value) or value<0 or value>limit:
                        raise RuntimeError('Python GPU runtime numerical gate failed')
            shutil.copy2(native, package / name)
            public_python = {key: python_binding[key] for key in (
                'schema','binary','binary_sha256','extension_suffix','pybind11_version',
                'pybind_header_sha256','ngsolve_version','matlab_engine_user_dependency',
                'runtime_import_checked','redistribution_accepted')}
            public_python['source_commit']=commit
            public_python['runtime_import_checked']=True
            public_python['runtime_validation_sha256']=hashlib.sha256(runtime_path.read_bytes()).hexdigest()
            public_python['loaded_runtime_sha256']=runtime['runtime_identity']['sha256']
            (package/'python-build.json').write_text(json.dumps(public_python,indent=2)+'\n')
            expected.update({name,'python-build.json'})
        hmatrix = os.environ.get('RADIA_HDIV_HMATRIX_BUILD')
        projection = os.environ.get('RADIA_HDIV_PROJECTION_BUILD')
        if bool(hmatrix) != bool(projection):
            raise RuntimeError('Select both HDiv HMatrix and projection build directories')
        if hmatrix:
            expected.update(payload_module().stage(ROOT, package, commit, hmatrix, projection))
        resident = os.environ.get('RADIA_RESIDENT_AMG_BUILD')
        if resident:
            expected.update(payload_module_resident().stage(ROOT, package, commit, resident))
        fmm = os.environ.get('RADIA_FMM_BUILD')
        if fmm:
            expected.update(payload_module().stage_fmm(ROOT, package, commit, fmm))
        hdiv_wrappers = set(payload_module().WRAPPERS)
        for group in ("+radia/+gpu", "+radia/+ltspice"):
            for path in (MATLAB / group).rglob("*.m"):
                if path.name == 'LaplaceFMM.m' and not fmm:
                    continue
                if path.name in hdiv_wrappers or path.name in {"FieldCloud.m", "FieldElements.m", "solveHDivRT0.m", "solveCurvedHDiv.m", "evaluateHDivBHMaterial.m", "evaluateHDivBHChunks.m"}:
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
        expected.update(payload_module().stage_documentation(PACKAGE, package))
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
        runtime_common, runtime_interfaces = payload_module().runtime_contract(
            bool(python_binding or resident))
        manifest = dict(schema="radia.gpu.distribution.v1", version=metadata["project"]["version"],
                        source_commit=commit, files=payload,
                        capabilities=["CUDA FP64 CSR Jacobi-PCG", "Experimental CUDA FP64 AMS", "Experimental reusable CUDA FP64 IC(0)-PCG", "Reusable affine linear HCurl matrix and AMS hierarchy", "AMS HDF5 exchange", "Simulink CUDA PCG block", "LTspice sampled-data integration"] + (["Explicit CSR CPU MEX and direct Python CUDA adapter (qualification required)"] if python_binding else []) + (["Bounded affine/P2 HDiv orders 1/2"] if hmatrix else []) + (["Resident GPU AMG Python adapter"] if resident else []) + (["Experimental Laplace FMM and flat P1 boundary maps"] if fmm else []),
                        runtime_requirements=runtime_common,
                        runtime_requirements_by_interface=runtime_interfaces)
        (package / "manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")


class bdist_wheel(_bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        receipt_path=ROOT/'build-cuda-mex/build.json'
        if ((receipt_path.is_file() and json.loads(receipt_path.read_text()).get('python_binding'))
                or os.environ.get('RADIA_RESIDENT_AMG_BUILD')):
            self.root_is_pure=False

    def get_tag(self):
        receipt_path=ROOT/'build-cuda-mex/build.json'
        python_selected = bool(os.environ.get('RADIA_RESIDENT_AMG_BUILD') or
            (receipt_path.is_file() and json.loads(receipt_path.read_text()).get('python_binding')))
        return payload_module().wheel_tag(python_selected, sys.platform, sys.version_info)


class BinaryDistribution(payload_module_native().NativePayloadDistribution, Distribution):
    """MEX and optional Python extensions belong in the platform library root."""


setup(distclass=BinaryDistribution, cmdclass={"build_py": build_py, "bdist_wheel": bdist_wheel})
