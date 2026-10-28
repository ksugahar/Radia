"""Verify the independent GPU wheel without importing numerical runtimes."""
from pathlib import Path, PurePosixPath
from email.parser import BytesParser
import argparse
import hashlib
import json
import re
import zipfile


# Full LGPL 2.1 text retained from the qualified upstream dependency wheels;
# normalize only CRLF so source checkout line-ending settings do not alter it.
_LGPL21_SHA256 = 'a190dc9c8043755d90f8b0a75fa66b9e42d4af4c980bf5ddc633f0124db3cee7'
_SEPARATE_GPU_DEPENDENCIES = (
    'ngsolve', 'netgen-mesher', 'nvidia-cuda-runtime-cu12',
    'nvidia-cublas-cu12', 'nvidia-cusolver-cu12', 'nvidia-cusparse-cu12',
    'nvidia-nvjitlink-cu12',
)


def audit_dependency_distribution(archive, metadata):
    """Check the user's dependency-only distribution conditions, not legal advice.

    A success applies to this archive alone; it does not authorize other builds
    or replace native, installed-runtime, licensing-provenance or release gates.
    """
    names = archive.namelist()
    if any(name.casefold().endswith('.dll') for name in names):
        raise ValueError('Dependency DLL bundling is forbidden')
    prefix = 'radia_simulink/'
    required = ('THIRD_PARTY_NOTICES.md', 'licenses/NGSOLVE-LICENSE.txt',
                'licenses/NETGEN-LICENSE.txt')
    if not all(prefix + name in names for name in required):
        raise ValueError('Dependency notice and full licenses required')
    notice = archive.read(prefix + required[0]).decode('utf-8')
    for url in ('https://github.com/NGSolve/ngsolve', 'https://github.com/NGSolve/netgen'):
        if url not in notice:
            raise ValueError('Dependency upstream source links required')
    for name in required[1:]:
        text = archive.read(prefix + name).decode('utf-8').replace('\r\n', '\n')
        if hashlib.sha256(text.encode('utf-8')).hexdigest() != _LGPL21_SHA256:
            raise ValueError('Dependency full LGPL 2.1 license text mismatch')
    requirements = metadata.get_all('Requires-Dist', [])
    pins = {}
    for name in _SEPARATE_GPU_DEPENDENCIES:
        entries = [value for value in requirements
                   if re.match(r'^' + re.escape(name) + r'\b', value, re.I)]
        if len(entries) != 1:
            raise ValueError('One exact separate GPU dependency pin required: ' + name)
        match = re.fullmatch(re.escape(name) + r'\s*(?:\(\s*)?==\s*'
                             r'([0-9][0-9A-Za-z.+-]*)\s*\)?\s*;\s*(.+)', entries[0], re.I)
        if match is None or not re.search(r'''extra\s*==\s*["']gpu["']''', match[2]):
            raise ValueError('Exact optional GPU dependency pin required: ' + name)
        pins[name] = match[1]
    if pins['ngsolve'] != pins['netgen-mesher']:
        raise ValueError('NGSolve and Netgen dependency ABI pins differ')
    return dict(passed=True, dependency_dlls_bundled=False,
                full_lgpl21_texts_checked=True, upstream_sources_checked=True,
                separate_dependency_pins=pins)


def verify_installed_identity(installed: dict, wheel: dict) -> None:
    """Reject a different build even when its source and version are identical."""
    digest = wheel.get('manifest_sha256', '')
    if (not re.fullmatch(r'[0-9a-f]{64}', digest)
            or any(installed.get(key) != wheel.get(key)
                   for key in ('source_commit', 'version', 'manifest_sha256'))):
        raise ValueError('Installed distribution differs from release wheel')


def verify_python_binding(path, names, archive, manifest):
    """Check the optional native adapter ABI; this is not runtime qualification."""
    python_binaries=[n for n in names if n.startswith('radia_simulink/') and n.endswith('.pyd')]
    tag=re.search(r'-(cp\d+)-(cp\d+)-win_amd64\.whl$',path.name)
    resident_receipt = 'radia_simulink/resident-amg-build.json'
    if resident_receipt in names and not any(
            re.fullmatch(r'radia_simulink/_resident_amg\.cp\d+-win_amd64\.pyd', n)
            for n in python_binaries):
        raise ValueError('Resident AMG receipt requires its Python binding')
    seen=set()
    for binary in python_binaries:
        filename=PurePosixPath(binary).name
        resident_tag = re.fullmatch(r'_resident_amg\.(cp\d+)-win_amd64\.pyd', filename)
        if resident_tag:
            if (not tag or tag[1] != tag[2] or tag[1] != resident_tag[1]
                    or 'resident-amg' in seen):
                raise ValueError('Resident AMG Python wheel ABI mismatch')
            seen.add('resident-amg')
            if resident_receipt not in names:
                raise ValueError('Resident AMG Python native receipt required')
            binding = json.loads(archive.read(resident_receipt))
            if (binding.get('schema') != 'radia.resident-amg-python-build.v1'
                    or binding.get('source_commit') != manifest['source_commit']
                    or binding.get('binary') != filename
                    or binding.get('binary_sha256') != manifest['files'].get(filename)
                    or binding.get('extension_suffix') != '.' + resident_tag[1] + '-win_amd64.pyd'
                    or any(binding.get(key) is not False for key in (
                        'source_dirty', 'fault_injection', 'device_fault_injection',
                        'lifecycle_diagnostic', 'matlab_engine_user_dependency', 'matlab_dlls_loaded'))
                    or any(binding.get(key) is not True for key in (
                        'runtime_import_checked', 'redistribution_accepted', 'matlab_dll_search_excluded'))
                    or any(not isinstance(binding.get(key), str) or not re.fullmatch(
                        r'[0-9a-f]{64}', binding[key]) for key in (
                            'runtime_validation_sha256', 'loaded_runtime_sha256',
                            'redistribution_reference_sha256'))):
                raise ValueError('Resident AMG runtime or provenance qualification missing')
            continue
        native_tag=re.fullmatch(r'_(csr)_cuda\.(cp\d+)-win_amd64\.pyd',filename)
        if not tag or not native_tag or tag[1]!=tag[2] or tag[1]!=native_tag[2] or native_tag[1] in seen:
            raise ValueError('Python CUDA binding wheel ABI mismatch')
        component=native_tag[1];seen.add(component)
        receipt_name='python-build.json'
        if 'radia_simulink/'+receipt_name not in names:raise ValueError('Python CUDA native receipt required')
        binding=json.loads(archive.read('radia_simulink/'+receipt_name))
        if (binding.get('schema')!='radia.'+component+'-python-build.v1'
                or binding.get('source_commit')!=manifest['source_commit']
                or binding.get('matlab_engine_user_dependency') is not False
                or binding.get('binary')!=filename
                or binding.get('binary_sha256')!=manifest['files'].get(filename)):
            raise ValueError('Python CUDA binding wheel ABI or hash mismatch')
    if not python_binaries and '-py3-none-' not in path.name:
        raise ValueError('Native ABI wheel is missing its Python binding')


def verify(path: Path) -> dict:
    if not re.fullmatch(r"radia_simulink-[^-]+-(?:py3-none|cp\d+-cp\d+)-win_amd64\.whl",path.name):
        raise ValueError("Expected a Windows x64 radia-simulink wheel")
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive members")
        for name in names:
            parsed = PurePosixPath(name)
            if parsed.is_absolute() or ".." in parsed.parts or "\\" in name:
                raise ValueError("Unsafe archive member")
            if (name.casefold().endswith((".dll", ".exe")) or "radia_mex.mexw64" in name
                    or (name.casefold().endswith(".pyd") and not re.fullmatch(r"radia_simulink/_(?:(?:csr|hdiv)_cuda|resident_amg)\.cp\d+-win_amd64\.pyd",name))):
                raise ValueError("CPU solver or external executable in GPU wheel")
            if not (name.startswith("radia_simulink/") or re.match(r"radia_simulink-[^/]+\.dist-info/",name)):
                raise ValueError("Foreign package in GPU wheel")
        metadata_paths = [n for n in names if n.endswith(".dist-info/METADATA")]
        if len(metadata_paths) != 1:
            raise ValueError("Expected one distribution metadata record")
        metadata = BytesParser().parsebytes(archive.read(metadata_paths[0]))
        distribution_audit = audit_dependency_distribution(archive, metadata)
        if metadata["Name"] != "radia-simulink" or metadata.get_all("Requires-Dist"):
            # Test extras are allowed but runtime solver dependencies are not.
            requirements = metadata.get_all("Requires-Dist",[])
            if metadata["Name"] != "radia-simulink" or any(not any(f'extra == "{extra}"' in r for extra in ('test', 'exchange', 'ltspice', 'coupling', 'gpu')) for r in requirements):
                raise ValueError("Unexpected GPU distribution identity/dependency")
        manifest_bytes = archive.read("radia_simulink/manifest.json")
        manifest = json.loads(manifest_bytes)
        if manifest.get("schema") != "radia.gpu.distribution.v1" or manifest["version"] != metadata["Version"]:
            raise ValueError("Manifest/metadata mismatch")
        expected = {"radia_simulink/"+n for n in manifest["files"]}
        actual = {n for n in names if n.startswith("radia_simulink/") and n != "radia_simulink/manifest.json"}
        if expected != actual:
            raise ValueError("GPU payload inventory mismatch")
        for relative, digest in manifest["files"].items():
            if hashlib.sha256(archive.read("radia_simulink/"+relative)).hexdigest() != digest:
                raise ValueError(f"GPU payload hash mismatch: {relative}")
        verify_python_binding(path,names,archive,manifest)
        binary = "matlab/radia_cuda_csr_pcg.mexw64"
        receipt = json.loads(archive.read("radia_simulink/build.json"))
        if receipt["source_dirty"] or receipt["source_commit"] != manifest["source_commit"]:
            raise ValueError("GPU source provenance mismatch")
        if receipt["binary_sha256"] != manifest["files"][binary]:
            raise ValueError("CUDA binary receipt mismatch")
        ams_receipt = json.loads(archive.read("radia_simulink/ams-build.json"))
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
            ic = json.loads(archive.read('radia_simulink/ic-build.json'))
            verify_ic_receipt(manifest, ic)
        hdiv_binaries = {'matlab/radia_cuda_hmatrix.mexw64', 'matlab/radia_cuda_projection.mexw64'}
        if hdiv_binaries.intersection(manifest['files']) or 'hdiv-build.json' in manifest['files']:
            required = hdiv_binaries | {'hdiv-build.json'} | {
                'matlab/+radia/+gpu/' + name for name in
                ('solveHDivElements.m', 'affineHDivQuadrature.m', 'curvedHDivQuadrature.m',
                 'HMatrixOperator.m', 'HDivBHMaterial.m', 'PreparedHDivNonlinear.m',
                 'makeHDivBHMaterial.m', 'solveNonlinearHDivElements.m')}
            if not required.issubset(manifest['files']):
                raise ValueError('Incomplete HDiv CUDA payload')
            hdiv = json.loads(archive.read('radia_simulink/hdiv-build.json'))
            verify_hdiv_receipt(manifest, hdiv)
        fmm_files = {'fmm-build.json', 'matlab/radia_cuda_fmm.mexw64', 'matlab/+radia/+gpu/LaplaceFMM.m'}
        if fmm_files & set(manifest['files']):
            if not fmm_files <= set(manifest['files']):
                raise ValueError('Incomplete FMM CUDA payload')
            fmm = json.loads(archive.read('radia_simulink/fmm-build.json'))
            if (fmm['source_dirty'] or fmm['source_commit'] != manifest['source_commit']
                    or fmm['taskmanager_used'] or fmm['python_adapter_required']
                    or fmm['binary_sha256'] != manifest['files']['matlab/radia_cuda_fmm.mexw64']):
                raise ValueError('FMM CUDA binary receipt mismatch')
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
                cuda_solvers=sorted(solvers), distribution_audit=distribution_audit)


def verify_hdiv_receipt(manifest, receipt):
    """Unknown diagnostic mode is not evidence of a release build."""
    if receipt.get('schema') != 'radia.gpu.hdiv-payload.v1':
        raise ValueError('Invalid HDiv payload schema')
    for name in ('hmatrix', 'projection'):
        native = receipt.get('artifacts', {}).get(name, {})
        if (native.get('source_dirty') is not False
                or native.get('source_commit') != manifest['source_commit']
                or native.get('taskmanager_started') is not False
                or native.get('fault_injection') is not False
                or (name == 'projection' and native.get('lifecycle_diagnostic') is not False)
                or native.get('binary_sha256') != manifest['files']['matlab/radia_cuda_' + name + '.mexw64']):
            raise ValueError('HDiv CUDA binary receipt mismatch')


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
