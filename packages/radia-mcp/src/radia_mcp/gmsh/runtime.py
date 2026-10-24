"""Qualified Gmsh 5 runtime and strict MSH 4.2 boundary.

Record a downloaded official wheel before use with
``python -m radia_mcp.gmsh.runtime --wheel WHEEL --record FILE`` and set
``RADIA_GMSH_PROVENANCE`` to that file. No installation occurs here.
"""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import tempfile
import zipfile

_INSTALL_GUIDANCE = (
    'Install Gmsh first: pip install --pre --upgrade -i '
    'https://gmsh.info/python-packages-dev gmsh; '
    'then pip install "radia-mcp[gmsh]". '
    'See docs/operations/gmsh-development-runtime.md.'
)


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def require_msh42(path, _seen=None) -> None:
    """Reject other mesh formats before passing a file to native Gmsh."""
    path = Path(path)
    if path.suffix.lower() == '.geo':
        seen = set() if _seen is None else _seen
        resolved = path.resolve()
        if resolved in seen:
            return
        seen.add(resolved)
        text = re.sub(r'//[^\n]*|/\*.*?\*/', '', path.read_text(encoding='utf-8'), flags=re.S)
        merges = re.findall(r'\b(?:Merge|Include)\s+([^;]+);', text)
        for expression in merges:
            literal = re.fullmatch(r'"([^"\n]+)"', expression.strip())
            if not literal:
                raise ValueError(f'{path.name}: dynamic Merge/Include target cannot be checked for MSH 4.2')
            require_msh42(path.parent / literal.group(1), seen)
        return
    if path.suffix.lower() != '.msh':
        return
    with path.open('rb') as stream:
        opening, header, closing = [stream.readline().strip() for _ in range(3)]
    if opening != b'$MeshFormat' or header.split() != [b'4.2', b'0', b'8'] or closing != b'$EndMeshFormat':
        raise ValueError(f'{path.name}: expected ASCII MSH 4.2 (4.2 0 8), got {header!r}')

    # Field tags are dense native indices in the adopted development build.
    # Reuse the pure-Python field decoder: mesh-only size_t tags remain valid.
    from .msh_inspect import _split_sections, _parse_data_section
    sections, errors = _split_sections(path.read_text(encoding='utf-8').splitlines())
    if errors:
        raise ValueError(f'{path.name}: {"; ".join(errors)}')
    for kind, body, line in sections:
        if kind not in ('NodeData', 'ElementData', 'ElementNodeData'):
            continue
        view, errors = _parse_data_section(kind, body, line)
        if errors:
            raise ValueError(f'{path.name}: {"; ".join(errors)}')
        for tag in view['tags']:
            if tag > 2**31 - 1:
                raise ValueError(f'{path.name}: Gmsh native {kind} field tag {tag} '
                                 'exceeds 2147483647; no remapping or fallback is used')


def _gmsh5_distribution_version() -> str:
    """Keep distribution-major validation independent of release channel."""
    version = importlib.metadata.version('gmsh')
    if version.split('.', 1)[0] != '5':
        raise RuntimeError(f'Gmsh 5 required; installed {version}. See the canonical Gmsh runtime install guide.')
    return version


def record_wheel(wheel: Path, destination: Path) -> dict:
    """Bind the downloaded wheel hash to the actually installed Python/DLL."""
    import gmsh
    version = _gmsh5_distribution_version()
    payload = {}
    with zipfile.ZipFile(wheel) as archive:
        metadata = [n for n in archive.namelist() if n.endswith('.dist-info/METADATA')]
        if len(metadata) != 1 or f'Version: {version}' not in archive.read(metadata[0]).decode().splitlines():
            raise RuntimeError('wheel version does not match installed Gmsh')
        for path in (Path(gmsh.__file__).resolve(), Path(gmsh.libpath).resolve()):
            members = [n for n in archive.namelist() if Path(n).name == path.name]
            if len(members) != 1:
                raise RuntimeError(f'wheel has ambiguous/missing runtime payload {path.name}')
            digest = hashlib.sha256(archive.read(members[0])).hexdigest()
            if _sha(path) != digest:
                raise RuntimeError(f'installed runtime differs from wheel: {path.name}')
            payload[str(path)] = digest
    result = {'distribution_version': version, 'wheel_sha256': _sha(wheel),
              'runtime_payload_sha256': payload}
    destination.write_text(json.dumps(result, indent=2), encoding='utf-8')
    return result


def verify_runtime_identity() -> dict:
    """Verify installed payload provenance without initializing or writing files."""
    try:
        version = _gmsh5_distribution_version()
    except importlib.metadata.PackageNotFoundError as error:
        raise RuntimeError(_INSTALL_GUIDANCE) from error
    record = os.environ.get('RADIA_GMSH_PROVENANCE')
    if not record:
        raise RuntimeError('RADIA_GMSH_PROVENANCE is required; record the official wheel with python -m radia_mcp.gmsh.runtime --wheel WHEEL --record FILE')
    identity = json.loads(Path(record).read_text(encoding='utf-8'))
    if identity.get('distribution_version') != version or not re.fullmatch(r'[0-9a-f]{64}', identity.get('wheel_sha256', '')):
        raise RuntimeError('Gmsh wheel provenance/version mismatch')
    try:
        import gmsh
    except ModuleNotFoundError as error:
        if error.name != 'gmsh':
            raise
        raise RuntimeError(_INSTALL_GUIDANCE) from error
    for path in (Path(gmsh.__file__).resolve(), Path(gmsh.libpath).resolve()):
        if identity.get('runtime_payload_sha256', {}).get(str(path)) != _sha(path):
            raise RuntimeError(f'Gmsh runtime payload/provenance mismatch: {path.name}')
    return identity


def qualify_runtime() -> dict:
    """Verify native identity, MSH default and provenance; initialize no session."""
    identity = verify_runtime_identity()
    import gmsh
    gmsh.initialize(['radia-gmsh-qualification'], readConfigFiles=False)
    try:
        if gmsh.option.getNumber('Mesh.MshFileVersion') != 4.2 or gmsh.__version__.split('.')[0] != '5':
            raise RuntimeError('Gmsh major 5 and default MSH 4.2 required')
        with tempfile.TemporaryDirectory(prefix='gmsh42_probe_') as directory:
            gmsh.option.setNumber('General.Terminal', 0)
            gmsh.model.add('runtime_probe')
            gmsh.model.addDiscreteEntity(1, 1)
            gmsh.model.mesh.addNodes(1, 1, [1, 2], [0, 0, 0, 1, 0, 0])
            gmsh.model.mesh.addElementsByType(1, 1, [1], [1, 2])
            view = gmsh.view.add('probe')
            gmsh.view.addModelData(view, 0, 'runtime_probe', 'NodeData', [1, 2], [[7.], [-2.]], numComponents=1)
            probe = str(Path(directory) / 'probe.msh')
            gmsh.option.setNumber('Mesh.Binary', 0)
            gmsh.write(probe)
            gmsh.option.setNumber('PostProcessing.SaveMesh', 0)
            gmsh.view.write(view, probe, append=True)
            require_msh42(probe)
            gmsh.clear()
            gmsh.open(probe)
            kind, tags, data, _, ncomp = gmsh.view.getModelData(gmsh.view.getTags()[0], 0)
            if kind != 'NodeData' or ncomp != 1 or dict(zip(map(int, tags), map(list, data))) != {1: [7.], 2: [-2.]}:
                raise RuntimeError('Gmsh MSH 4.2 field read/write probe failed')
            gmsh.option.setNumber('PostProcessing.SaveMesh', 1)
        identity = {**identity, 'gmsh_version': gmsh.option.getString('General.Version'), 'msh_version': 4.2}
    finally:
        gmsh.finalize()
    install_io_guards(gmsh)
    return {**identity, 'msh42_roundtrip_probe': True}


def install_io_guards(gmsh) -> None:
    """Guard native mesh reads/writes, including explicit session open/merge."""
    if getattr(gmsh, '_radia_msh42_guarded', False):
        return
    def wrap_read(original):
        def checked(path, *args, **kwargs):
            require_msh42(path)
            return original(path, *args, **kwargs)
        return checked
    gmsh.open = wrap_read(gmsh.open)
    gmsh.merge = wrap_read(gmsh.merge)
    original_write, original_view_write = gmsh.write, gmsh.view.write
    def write(path, *args, **kwargs):
        if Path(path).suffix.lower() == '.msh':
            gmsh.option.setNumber('Mesh.MshFileVersion', 4.2)
            gmsh.option.setNumber('Mesh.Binary', 0)
        result = original_write(path, *args, **kwargs)
        require_msh42(path)
        return result
    def view_write(tag, path, append=False):
        mesh = Path(path).suffix.lower() == '.msh'
        previous = gmsh.option.getNumber('PostProcessing.SaveMesh')
        if mesh:
            gmsh.option.setNumber('Mesh.MshFileVersion', 4.2)
            gmsh.option.setNumber('Mesh.Binary', 0)
            gmsh.option.setNumber('PostProcessing.SaveMesh', 0 if append else 1)
        try:
            result = original_view_write(tag, path, append)
            require_msh42(path)
            return result
        finally:
            if mesh:
                gmsh.option.setNumber('PostProcessing.SaveMesh', previous)
    gmsh.write, gmsh.view.write = write, view_write
    gmsh._radia_msh42_guarded = True


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--wheel', type=Path, required=True)
    parser.add_argument('--record', type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(record_wheel(args.wheel, args.record)))


if __name__ == '__main__':
    main()
