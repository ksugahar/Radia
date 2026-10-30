"""Reject GPU/product integration payloads in the independently released CPU wheel."""
from __future__ import annotations

import argparse
from email.parser import Parser
import json
from pathlib import Path, PurePosixPath
import re
import zipfile


def verify(path: Path, *, require_python_gpu=False) -> dict:
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate wheel archive member")
        metadata_paths = [name for name in names if name.endswith('.dist-info/METADATA')]
        if len(metadata_paths) != 1:
            raise ValueError('Expected exactly one CPU distribution metadata record')
        metadata = Parser().parsestr(archive.read(metadata_paths[0]).decode('utf-8'))
        if metadata['Name'] != 'radia':
            raise ValueError('CPU wheel must be the radia distribution')
        forbidden = []
        retired_equation = []
        for name in names:
            member = PurePosixPath(name.lower())
            if (member.parts[:2] == ('radia', 'equation')
                    or (member.parent == PurePosixPath('radia')
                        and member.name.startswith('_equation')
                        and member.suffix in {'.pyd', '.so', '.dylib'})):
                retired_equation.append(name)
            if (set(member.parts) & {'matlab', 'simulink', 'ltspice', 'radia_simulink'}
                    or member.stem in {'matlab', 'simulink', 'ltspice', 'radia_simulink'}
                    or member.suffix in {'.m', '.mlx', '.slx', '.mdl', '.mexw64',
                                         '.mexa64', '.mexmaci64', '.mexmaca64', '.cu', '.cuh'}
                    or any(suffix.startswith('.mex') for suffix in member.suffixes)
                    or member.name in {'matlab_bridge.py', 'vim_matlab_bridge.py',
                                       'application_interface_manifest.json'}):
                forbidden.append(name)
        gpu_records = []
        if require_python_gpu:
            import importlib.util
            spec = importlib.util.spec_from_file_location("radia_python_gpu_payload", Path(__file__).with_name("python_gpu_payload.py"))
            payload = importlib.util.module_from_spec(spec); spec.loader.exec_module(payload)
            if not {"radia/_radia_pybind.pyd", "radia/_radia_pybind.pyd.build.json"}.issubset(names):
                raise ValueError("CPU native binary/provenance is missing")
            core = json.loads(archive.read("radia/_radia_pybind.pyd.build.json"))
            source = core.get("source_commit", "")
            if core.get("source_dirty") is not False or core.get("binary_sha256") != __import__("hashlib").sha256(archive.read("radia/_radia_pybind.pyd")).hexdigest():
                raise ValueError("CPU build provenance is missing or mismatched")
            tag = re.search(r"-(cp[0-9]+)-\1-win_amd64\.whl$", path.name)
            if not tag:
                raise ValueError("Python GPU payload requires a Windows CPython ABI wheel")
            for module in payload.MODULES:
                expected = "radia/" + module + "." + tag[1] + "-win_amd64.pyd"
                candidates = [n for n in names if PurePosixPath(n).name.startswith(module + ".") and n.endswith((".pyd", ".so"))]
                if candidates != [expected]:
                    raise ValueError("Missing or ambiguous required Python GPU binding: " + module)
                record = json.loads(archive.read(expected + ".build.json"))
                payload.verify_payload(PurePosixPath(expected).name, archive.read(expected), record, source)
                gpu_records.append(expected)
            if any(n.endswith(".private-build.json") or PurePosixPath(n).name in ("admission-private.json", "commands.json", "hdiv-sdk-private.json") for n in names):
                raise ValueError("Private GPU build originals must not ship")
        for requirement in metadata.get_all('Requires-Dist', []):
            dependency = re.split(r'[\s\[<>=!~;(]', requirement, maxsplit=1)[0]
            dependency = re.sub(r'[-_.]+', '-', dependency).lower()
            if dependency in {'radia-simulink', 'matlab', 'matlabengine', 'matlab-engine-for-python'}:
                forbidden.append('Requires-Dist: ' + requirement)
        if retired_equation:
            raise ValueError('CPU wheel contains retired equation payload: ' + ', '.join(retired_equation))
        if forbidden:
            raise ValueError('CPU wheel contains radia-simulink integration payload: ' + ', '.join(forbidden))
    return dict(passed=True, distribution='radia', version=metadata['Version'],
                gpu_integration_absent=True, python_gpu_required=require_python_gpu,
                python_gpu_binaries=gpu_records, runtime_qualified=False,
                retired_equation_absent=True, files=len(names))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    parser.add_argument('--require-python-gpu', action='store_true')
    args = parser.parse_args()
    print(json.dumps(verify(args.wheel, require_python_gpu=args.require_python_gpu), indent=2))
