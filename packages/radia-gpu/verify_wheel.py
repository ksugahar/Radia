"""Verify the independent GPU wheel without importing numerical runtimes."""
from pathlib import Path, PurePosixPath
from email.parser import BytesParser
import argparse
import hashlib
import json
import re
import zipfile


def verify_ams_distribution_eligibility(receipt):
    if (receipt.get('geometry_lifecycle_test_build') is True
            or receipt.get('distribution_eligible') is not True):
        raise ValueError('AMS test or unclassified build is not distribution eligible')


def verify_installed_identity(installed: dict, wheel: dict) -> None:
    """Reject a different build even when its source and version are identical."""
    digest = wheel.get('manifest_sha256', '')
    if (not re.fullmatch(r'[0-9a-f]{64}', digest)
            or any(installed.get(key) != wheel.get(key)
                   for key in ('source_commit', 'version', 'manifest_sha256'))):
        raise ValueError('Installed distribution differs from release wheel')


def verify(path: Path) -> dict:
    if not re.fullmatch(r"radia_gpu-[^-]+-py3-none-win_amd64\.whl",path.name):
        raise ValueError("Expected a Windows x64 radia-gpu wheel")
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive members")
        for name in names:
            parsed = PurePosixPath(name)
            if parsed.is_absolute() or ".." in parsed.parts or "\\" in name:
                raise ValueError("Unsafe archive member")
            if name.endswith((".dll", ".exe", ".pyd")) or "radia_mex.mexw64" in name:
                raise ValueError("CPU solver or external executable in GPU wheel")
            if not (name.startswith("radia_gpu/") or re.match(r"radia_gpu-[^/]+\.dist-info/",name)):
                raise ValueError("Foreign package in GPU wheel")
        metadata_paths = [n for n in names if n.endswith(".dist-info/METADATA")]
        if len(metadata_paths) != 1:
            raise ValueError("Expected one distribution metadata record")
        metadata = BytesParser().parsebytes(archive.read(metadata_paths[0]))
        if metadata["Name"] != "radia-gpu" or metadata.get_all("Requires-Dist"):
            # Test extras are allowed but runtime solver dependencies are not.
            requirements = metadata.get_all("Requires-Dist",[])
            if metadata["Name"] != "radia-gpu" or any(not any(f'extra == "{extra}"' in r for extra in ('test', 'exchange', 'ltspice', 'coupling')) for r in requirements):
                raise ValueError("Unexpected GPU distribution identity/dependency")
        manifest_bytes = archive.read("radia_gpu/manifest.json")
        manifest = json.loads(manifest_bytes)
        if manifest.get("schema") != "radia.gpu.distribution.v1" or manifest["version"] != metadata["Version"]:
            raise ValueError("Manifest/metadata mismatch")
        expected = {"radia_gpu/"+n for n in manifest["files"]}
        actual = {n for n in names if n.startswith("radia_gpu/") and n != "radia_gpu/manifest.json"}
        if expected != actual:
            raise ValueError("GPU payload inventory mismatch")
        for relative, digest in manifest["files"].items():
            if hashlib.sha256(archive.read("radia_gpu/"+relative)).hexdigest() != digest:
                raise ValueError(f"GPU payload hash mismatch: {relative}")
        binary = "matlab/radia_cuda_csr_pcg.mexw64"
        receipt = json.loads(archive.read("radia_gpu/build.json"))
        if receipt["source_dirty"] or receipt["source_commit"] != manifest["source_commit"]:
            raise ValueError("GPU source provenance mismatch")
        if receipt["binary_sha256"] != manifest["files"][binary]:
            raise ValueError("CUDA binary receipt mismatch")
        ams_receipt = json.loads(archive.read("radia_gpu/ams-build.json"))
        verify_ams_distribution_eligibility(ams_receipt)
        if (ams_receipt["source_dirty"] or ams_receipt["source_commit"] != manifest["source_commit"]
                or ams_receipt["binary_sha256"] != manifest["files"]["matlab/radia_cuda_ams.mexw64"]):
            raise ValueError("AMS CUDA binary receipt mismatch")
        for required in (binary,"matlab/+radia/+gpu/csrPCG.m", "matlab/+radia/+simulink/cudaPCGSFunction.m",
                         "matlab/+radia/+ltspice/run.m", "matlab/radia_cuda_ams.mexw64",
                         "matlab/+radia/+gpu/readAMSExchange.m", "matlab/+radia/+gpu/amsPCG.m",
                         "matlab/+radia/+gpu/HCurlSystem.m", "matlab/+radia/+gpu/solveHCurl.m",
                         "ams_exchange.py", "ltspice/__init__.py", "ltspice/cli.py", "LICENSE"):
            if required not in manifest["files"]:
                raise ValueError(f"Missing GPU capability payload: {required}")
        ic_files = {'ic-build.json','matlab/radia_cuda_ic.mexw64',
                    'matlab/+radia/+gpu/ICContext.m','matlab/+radia/+gpu/HCurlICSystem.m'}
        if ic_files & set(manifest['files']):
            if not ic_files <= set(manifest['files']):
                raise ValueError('Incomplete IC CUDA payload')
            ic = json.loads(archive.read('radia_gpu/ic-build.json'))
            verify_ic_receipt(manifest, ic)
        hdiv_binaries = {'matlab/radia_cuda_hmatrix.mexw64', 'matlab/radia_cuda_projection.mexw64'}
        if hdiv_binaries.intersection(manifest['files']) or 'hdiv-build.json' in manifest['files']:
            required = hdiv_binaries | {'hdiv-build.json'} | {
                'matlab/+radia/+gpu/' + name for name in
                ('solveHDivElements.m', 'affineHDivQuadrature.m', 'curvedHDivQuadrature.m')}
            if not required.issubset(manifest['files']):
                raise ValueError('Incomplete HDiv CUDA payload')
            hdiv = json.loads(archive.read('radia_gpu/hdiv-build.json'))
            if hdiv.get('schema') != 'radia.gpu.hdiv-payload.v1':
                raise ValueError('Invalid HDiv payload schema')
            for name in ('hmatrix', 'projection'):
                native = hdiv['artifacts'][name]
                if (native['source_dirty'] or native['source_commit'] != manifest['source_commit']
                        or native['taskmanager_started']
                        or native['binary_sha256'] != manifest['files']['matlab/radia_cuda_' + name + '.mexw64']):
                    raise ValueError('HDiv CUDA binary receipt mismatch')
        fmm_files = {'fmm-build.json', 'matlab/radia_cuda_fmm.mexw64', 'matlab/+radia/+gpu/LaplaceFMM.m'}
        if fmm_files & set(manifest['files']):
            if not fmm_files <= set(manifest['files']):
                raise ValueError('Incomplete FMM CUDA payload')
            fmm = json.loads(archive.read('radia_gpu/fmm-build.json'))
            if (fmm['source_dirty'] or fmm['source_commit'] != manifest['source_commit']
                    or fmm['taskmanager_used'] or fmm['python_adapter_required']
                    or fmm['binary_sha256'] != manifest['files']['matlab/radia_cuda_fmm.mexw64']):
                raise ValueError('FMM CUDA binary receipt mismatch')
        ih_files = {'ih-build.json', 'matlab/radia_cuda_ih.mexw64',
                    'matlab/radia_ih_cuda_sfun.m',
                    'matlab/+radia/+simulink/ihResidentSFunction.m',
                    'matlab/+radia/+simulink/makeCompletePEECBEMConfig.m',
                    'matlab/+radia/+simulink/buildTemperatureBHResearchModel.m',
                    'matlab/samples/ih-temperature-bh/makeTemperatureBHConfig.m',
                    'matlab/samples/ih-temperature-bh/runTemperatureBHResearch.m',
                    'IH_TEMPERATURE_BH.md', 'IH_EXPERIMENT.md'}
        if ih_files & set(manifest['files']):
            if not ih_files <= set(manifest['files']):
                raise ValueError('Incomplete IH research payload')
            ih = json.loads(archive.read('radia_gpu/ih-build.json'))
            if (ih.get('schema') != 'radia.gpu.ih-research-payload.v1'
                    or ih.get('source_dirty') is not False
                    or ih.get('source_commit') != manifest['source_commit']
                    or ih.get('taskmanager_regions') is not False
                    or ih.get('python_runtime_required') is not False
                    or ih.get('cpu_fallback') is not False
                    or not {'element-temperature-bh','complete-flat-p1-reaction','accepted-element-impedance','transactional-temperature-range'} <= set(ih.get('features', []))
                    or ih.get('research_only') is not True
                    or ih.get('experimental_validation') is not False
                    or ih.get('independent_fem_validation') is not False
                    or ih.get('binary_sha256') != manifest['files']['matlab/radia_cuda_ih.mexw64']):
                raise ValueError('IH research binary receipt mismatch')
    solvers = []
    for relative in manifest['files']:
        if relative.endswith('.mexw64'):
            name = Path(relative).name
            if not name.startswith('radia_cuda_'):
                raise ValueError('Unclassified CUDA solver in GPU wheel: '+relative)
            solvers.append(name.removeprefix('radia_cuda_').removesuffix('.mexw64'))
    return dict(passed=True, version=metadata["Version"], source_commit=manifest["source_commit"],
                wheel_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), files=len(expected),
                manifest_sha256=hashlib.sha256(manifest_bytes).hexdigest(),
                cuda_solvers=sorted(solvers))


def verify_ic_receipt(manifest, receipt):
    if (receipt.get('schema') != 'radia.gpu.ic-payload.v1'
            or receipt.get('source_dirty') is not False
            or receipt.get('source_commit') != manifest['source_commit']
            or receipt.get('gpu_numerical_factorization') is not True
            or receipt.get('no_cpu_solver_fallback') is not True
            or receipt.get('binary_sha256') != manifest['files']['matlab/radia_cuda_ic.mexw64']):
        raise ValueError('IC CUDA binary receipt mismatch')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel",type=Path)
    options = parser.parse_args()
    print(json.dumps(verify(options.wheel),indent=2))
