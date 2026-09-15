"""File-only bridge for legacy Cubit hooks; never import either native package.

New installs point directly at the exporter-owned assets. This bridge keeps
old startup paths usable until the user runs cubit-plugin-install again.
"""
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess


def load_exporter_file(relative, namespace):
    try:
        distribution = importlib.metadata.distribution('cubit-mesh-export')
        root = Path(distribution.locate_file('cubit_mesh_export'))
        # PEP 660 editable distributions need their recorded source location.
        direct = distribution.read_text('direct_url.json')
        if not root.is_dir() and direct:
            from urllib.parse import unquote, urlparse
            source = unquote(urlparse(json.loads(direct)['url']).path)
            if os.name == 'nt' and source.startswith('/') and source[2:3] == ':':
                source = source[1:]
            root = Path(source) / 'src' / 'cubit_mesh_export'
    except importlib.metadata.PackageNotFoundError:
        # Cubit's private Python must not import external CPython extensions.
        python = os.environ.get('CUBIT_MESH_EXPORT_PYTHON') or os.environ.get('RADIA_PYTHON')
        if not python:
            raise ModuleNotFoundError('Legacy Cubit GUI hook: run cubit-plugin-install with the external Python, then restart Cubit.') from None
        query = ('import importlib.util,json; '
                 's=importlib.util.find_spec("cubit_mesh_export"); '
                 'print(json.dumps(s.origin if s else None))')
        result = subprocess.run([python, '-c', query], check=True,
                                capture_output=True, text=True, encoding='utf-8', timeout=15)
        origin = json.loads(result.stdout)
        if not origin:
            raise RuntimeError('cubit-mesh-export is missing from the configured external Python')
        root = Path(origin).parent
    path = root / relative
    if not path.is_file():
        raise ImportError('Exporter-owned GUI assets are missing; update cubit-mesh-export and run cubit-plugin-install')
    namespace['__file__'] = str(path)
    exec(compile(path.read_text(encoding='utf-8'), str(path), 'exec'), namespace)
