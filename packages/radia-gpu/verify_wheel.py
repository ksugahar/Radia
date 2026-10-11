"""Verify the independent GPU wheel without importing numerical runtimes."""
from pathlib import Path, PurePosixPath
from email.parser import BytesParser
import argparse
import hashlib
import json
import re
import zipfile


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
            if metadata["Name"] != "radia-gpu" or any(not any(f'extra == "{extra}"' in r for extra in ('test', 'ltspice', 'coupling')) for r in requirements):
                raise ValueError("Unexpected GPU distribution identity/dependency")
        manifest = json.loads(archive.read("radia_gpu/manifest.json"))
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
        for required in (binary,"matlab/+radia/+gpu/csrPCG.m", "matlab/+radia/+simulink/cudaPCGSFunction.m",
                         "matlab/+radia/+ltspice/run.m", "ltspice/__init__.py", "ltspice/cli.py", "LICENSE"):
            if required not in manifest["files"]:
                raise ValueError(f"Missing GPU capability payload: {required}")
    solvers = []
    for relative in manifest['files']:
        if relative.endswith('.mexw64'):
            name = Path(relative).name
            if not name.startswith('radia_cuda_'):
                raise ValueError('Unclassified CUDA solver in GPU wheel: '+relative)
            solvers.append(name.removeprefix('radia_cuda_').removesuffix('.mexw64'))
    return dict(passed=True, version=metadata["Version"], source_commit=manifest["source_commit"],
                wheel_sha256=hashlib.sha256(path.read_bytes()).hexdigest(), files=len(expected),
                cuda_solvers=sorted(solvers))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel",type=Path)
    options = parser.parse_args()
    print(json.dumps(verify(options.wheel),indent=2))
