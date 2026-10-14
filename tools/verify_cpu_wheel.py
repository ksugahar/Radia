"""Reject GPU/product integration payloads in the independently released CPU wheel."""
from __future__ import annotations

import argparse
from email.parser import Parser
import json
from pathlib import Path, PurePosixPath
import re
import zipfile


def verify(path: Path) -> dict:
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
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
            if (set(member.parts) & {'matlab', 'simulink', 'ltspice', 'radia_gpu'}
                    or member.stem in {'matlab', 'simulink', 'ltspice', 'radia_gpu'}
                    or member.suffix in {'.m', '.mlx', '.slx', '.mdl', '.mexw64',
                                         '.mexa64', '.mexmaci64', '.mexmaca64', '.cu', '.cuh'}
                    or any(suffix.startswith('.mex') for suffix in member.suffixes)
                    or member.name in {'matlab_bridge.py', 'vim_matlab_bridge.py',
                                       'application_interface_manifest.json'}):
                forbidden.append(name)
        for requirement in metadata.get_all('Requires-Dist', []):
            dependency = re.split(r'[\s\[<>=!~;(]', requirement, maxsplit=1)[0]
            dependency = re.sub(r'[-_.]+', '-', dependency).lower()
            if dependency in {'radia-gpu', 'matlab', 'matlabengine', 'matlab-engine-for-python'}:
                forbidden.append('Requires-Dist: ' + requirement)
        if retired_equation:
            raise ValueError('CPU wheel contains retired equation payload: ' + ', '.join(retired_equation))
        if forbidden:
            raise ValueError('CPU wheel contains radia-gpu integration payload: ' + ', '.join(forbidden))
    return dict(passed=True, distribution='radia', version=metadata['Version'],
                gpu_integration_absent=True, retired_equation_absent=True, files=len(names))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('wheel', type=Path)
    args = parser.parse_args()
    print(json.dumps(verify(args.wheel), indent=2))
