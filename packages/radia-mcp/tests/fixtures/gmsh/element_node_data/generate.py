"""Regenerate tiny post-only fixtures with the qualified native Gmsh API.

No solver input or generated solver mesh is created. Run from a checkout with
the declared Gmsh runtime and provenance configured. Native writes, then native
reads the file back; JSON expectations come from getModelData, not our decoder.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json

from radia_mcp.gmsh.runtime import qualify_runtime


def main():
    runtime = qualify_runtime()
    import gmsh
    root = Path(__file__).resolve().parent
    fixtures = []
    cases = [(0, 15, 1, 1, 'point_scalar'), (0, 15, 1, 3, 'point_vector'),
             (1, 1, 2, 1, 'line_scalar'), (1, 1, 2, 3, 'line_vector'),
             (1, 1, 1, 3, 'line_sample_vector')]
    for dim, etype, samples, components, name in cases:
        gmsh.initialize(['-noconfig'])
        try:
            gmsh.option.setNumber('General.Verbosity', 0)
            gmsh.option.setNumber('Mesh.MshFileVersion', 4.2)
            gmsh.model.add(name)
            gmsh.model.addDiscreteEntity(dim, 1)
            nodes = [17, 29] if dim == 0 else [17, 29, 41]
            xyz = [coordinate for i in range(len(nodes)) for coordinate in (i, 0, 0)]
            gmsh.model.mesh.addNodes(dim, 1, nodes, xyz)
            connectivity = nodes if dim == 0 else [17, 29, 29, 41]
            gmsh.model.mesh.addElementsByType(1, etype, [11, 23], connectivity)
            tag = gmsh.view.add(name)
            values = [[(element + 1) * 10 + sample + .25 * component
                       for sample in range(samples) for component in range(components)]
                      for element in range(2)]
            # Nonmonotone association detects tag/value ordering mistakes.
            gmsh.view.addModelData(tag, 0, name, 'ElementNodeData', [23, 11],
                                   values, time=.5, numComponents=components)
            path = root / (name + '.msh')
            gmsh.view.write(tag, str(path))
            original = path.read_bytes()
            normalized = original.replace(b'\r\n', b'\n')
            path.write_bytes(normalized)
            gmsh.clear()
            gmsh.open(str(path))
            actual = gmsh.view.getModelData(gmsh.view.getTags()[-1], 0)
            fixtures.append({
                'file': path.name, 'element_type': etype, 'samples_per_element': samples,
                'components': int(actual[4]), 'time': float(actual[3]),
                'data_type': actual[0], 'tags': [int(t) for t in actual[1]],
                'values': [list(map(float, row)) for row in actual[2]],
                'sha256': hashlib.sha256(normalized).hexdigest(),
                'native_bytes_sha256': hashlib.sha256(original).hexdigest()})
        finally:
            gmsh.finalize()
    provenance = {
        'schema': 'radia.native-gmsh-fixtures.v1',
        'purpose': 'Post-only ElementNodeData parser contract; no solver mesh or numerical baseline',
        'generator': 'generate.py',
        'generator_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'gmsh_distribution': importlib.metadata.version('gmsh'),
        'native_version': runtime['gmsh_version'],
        'wheel_sha256': runtime['wheel_sha256'],
        'text_storage': 'Native output normalized CRLF to LF only; fixtures use local eol=lf attributes',
        'fixtures': fixtures}
    (root / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(provenance), flush=True)


if __name__ == '__main__':
    main()
