"""Seal a native linear-mesh capture outside the CAD process.

Guards: ordered connectivity, independent outward skin, explicit units and bytes.
Assumes: qualified native acquisition; producer claims do not authenticate it.
No CAD binding, interpreter embedding, application connection or solver call.
"""
from __future__ import annotations

import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path


FACES = {
    'TET4': ((1, 2, 3), (0, 3, 2), (0, 1, 3), (0, 2, 1)),
    'HEX8': ((0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
             (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)),
}


def encoded(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'),
                      ensure_ascii=False, allow_nan=False).encode('utf-8')


def digest(data):
    return hashlib.sha256(data).hexdigest()


def strict_json(data):
    def pairs(items):
        result = {}
        for key, value in items:
            if key in result:
                raise ValueError('Duplicate JSON key: ' + key)
            result[key] = value
        return result
    def invalid(value):
        raise ValueError('Nonfinite JSON constant: ' + value)
    return json.loads(data.decode('utf-8'), object_pairs_hook=pairs,
                      parse_constant=invalid)


def _positive(value):
    if type(value) is not int or value <= 0:
        raise ValueError('Positive integer ID required')
    return value


def decode_exodus(path):
    """Decode classic Exodus independently of the native collector.

    A missing node/element numbering map has the Exodus-defined sequential
    identity. Unsupported NetCDF encodings fail; no alternative file is guessed.
    """
    from scipy.io import netcdf_file
    with netcdf_file(path, mode='r', mmap=False) as db:
        var = db.variables
        n = int(db.dimensions['num_nodes'])
        node_ids = var['node_num_map'][:].tolist() if 'node_num_map' in var else list(range(1, n + 1))
        if 'coord' in var:
            xyz = var['coord'][:].T.tolist()
        else:
            xyz = list(zip(*(var['coord' + axis][:].tolist() for axis in 'xyz')))
        if len(node_ids) != n or len(set(node_ids)) != n or len(xyz) != n:
            raise ValueError('Native node numbering map is not bijective')
        nodes = {int(k): list(v) for k, v in zip(node_ids, xyz)}
        count = int(db.dimensions['num_elem'])
        element_ids = var['elem_num_map'][:].tolist() if 'elem_num_map' in var else list(range(1, count + 1))
        if len(element_ids) != count or len(set(element_ids)) != count:
            raise ValueError('Native element numbering map is not bijective')
        block_ids = var['eb_prop1'][:].tolist()
        cells, offset = {}, 0
        for index, region in enumerate(block_ids, 1):
            conn = var['connect' + str(index)]
            raw_type = conn.elem_type.decode('ascii').strip().upper()
            width = conn.shape[1]
            if raw_type in ('HEX', 'HEX8') and width == 8:
                kind = 'HEX8'
            elif raw_type in ('TETRA', 'TETRA4', 'TET4') and width == 4:
                kind = 'TET4'
            else:
                raise ValueError('Unsupported native element: ' + raw_type)
            for row in conn[:].tolist():
                eid = int(element_ids[offset]); offset += 1
                key = (kind, eid)
                if key in cells:
                    raise ValueError('Duplicate native typed element')
                if any(not 1 <= int(i) <= n for i in row):
                    raise ValueError('Native connectivity index outside node table')
                cells[key] = (int(region), [int(node_ids[int(i) - 1]) for i in row])
        if offset != count:
            raise ValueError('Incomplete native block coverage')
        return nodes, cells


def _cycle(face):
    face = tuple(face)
    i = face.index(min(face))
    return face[i:] + face[:i]


def _det(a, b, c):
    return (a[0] * (b[1]*c[2] - b[2]*c[1])
            - a[1] * (b[0]*c[2] - b[2]*c[0])
            + a[2] * (b[0]*c[1] - b[1]*c[0]))


def _positive_jacobians(kind, coordinates):
    origin = coordinates[0]
    p = [[v[k] - origin[k] for k in range(3)] for v in coordinates]
    if kind == 'TET4':
        values = [_det(p[1], p[2], p[3])]
    else:
        signs = ((-1,-1,-1), (1,-1,-1), (1,1,-1), (-1,1,-1),
                 (-1,-1,1), (1,-1,1), (1,1,1), (-1,1,1))
        g = 1 / math.sqrt(3)
        samples = [(x,y,z) for x in (-g,g) for y in (-g,g) for z in (-g,g)]
        samples += [(x,y,z) for x in (-1,0,1) for y in (-1,0,1) for z in (-1,0,1)]
        values = []
        for x, y, z in samples:
            dn = [(sx*(1+sy*y)*(1+sz*z)/8, sy*(1+sx*x)*(1+sz*z)/8,
                   sz*(1+sx*x)*(1+sy*y)/8) for sx,sy,sz in signs]
            matrix = [[math.fsum(v[k]*d[j] for v,d in zip(p,dn)) for k in range(3)] for j in range(3)]
            values.append(_det(*matrix))
    if any(not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Nonpositive sampled element Jacobian')


def check_capture(capture, native_path):
    if capture['schema'] != 'cubit.native-mesh-capture.v1' or capture['units'] not in ('m', 'mm'):
        raise ValueError('Explicit supported native capture/unit required')
    nodes = {}
    for node in capture['nodes']:
        nid = _positive(node['id'])
        if nid in nodes or len(node['xyz']) != 3 or any(type(x) not in (int, float) or not math.isfinite(x) for x in node['xyz']):
            raise ValueError('Duplicate/nonfinite native node')
        nodes[nid] = node['xyz']
    regions = [_positive(x['id']) for x in capture['regions']]
    if not regions or len(set(regions)) != len(regions):
        raise ValueError('Empty/duplicate regions')
    cells, owners, used = {}, defaultdict(list), set()
    for cell in capture['elements']:
        key = (cell['type'], _positive(cell['source_id']))
        conn = cell['node_ids']; kind = cell['type']
        if kind not in FACES or len(conn) != (4 if kind == 'TET4' else 8) or len(set(conn)) != len(conn):
            raise ValueError('Unsupported/repeated connectivity')
        if key in cells or cell['region_id'] not in regions or any(_positive(x) not in nodes for x in conn):
            raise ValueError('Duplicate/unresolved cell or region')
        cells[key] = (cell['region_id'], conn); used.update(conn)
        _positive_jacobians(kind, [nodes[n] for n in conn])
        for indices in FACES[kind]:
            face = [conn[i] for i in indices]
            owners[tuple(sorted(face))].append(_cycle(face))
    if not cells or used != set(nodes):
        raise ValueError('Empty/incomplete node coverage')
    for faces in owners.values():
        if len(faces) > 2 or (len(faces) == 2 and faces[0] != _cycle(tuple(reversed(faces[1])))):
            raise ValueError('Nonmanifold/same-oriented shared face')
    exterior = {key: faces[0] for key, faces in owners.items() if len(faces) == 1}
    skin = {}
    identities = set()
    for face in capture['faces']:
        identity = (_positive(face['source_surface_id']), _positive(face['source_face_id']))
        _positive(face['source_volume_id'])
        if identity in identities or face['source_volume_sense'] not in ('Forward', 'Reversed'):
            raise ValueError('Unknown/duplicate independent skin identity or sense')
        identities.add(identity)
        key = tuple(sorted(face['node_ids']))
        if key in skin:
            raise ValueError('Duplicate independent skin face')
        skin[key] = _cycle(face['node_ids'])
    if not skin or skin != exterior:
        raise ValueError('Independent outward skin does not match volume boundary')
    native_nodes, native_cells = decode_exodus(native_path)
    if native_nodes != nodes or native_cells != cells:
        raise ValueError('Independent native coordinates/typed ordered connectivity differ')
    return nodes, cells


def _write_bdf(nodes, cells, regions):
    lines = ['$ Neutral linear mesh; materials are placeholders.', 'CEND', 'BEGIN BULK']
    for nid, xyz in sorted(nodes.items()):
        values = [format(float(x), '.17g') for x in xyz]
        values = [s if '.' in s or 'e' in s.lower() else s + '.' for s in values]
        lines.append('GRID,' + str(nid) + ',,' + ','.join(values))
    for rid in sorted(regions):
        lines.extend([f'PSOLID,{rid},{rid}', f'MAT1,{rid},1.,,0.'])
    mapping = []
    for eid, ((kind, source), (region, conn)) in enumerate(sorted(cells.items()), 1):
        if eid > 999999:
            raise ValueError('Continuation label width exceeded')
        if kind == 'HEX8':
            lines.extend([f'CHEXA,{eid},{region},' + ','.join(map(str, conn[:6])) + f',+N{eid}',
                          f'+N{eid},' + ','.join(map(str, conn[6:]))])
        else:
            lines.append(f'CTETRA,{eid},{region},' + ','.join(map(str, conn)))
        mapping.append({'source_id': source, 'source_type': kind, 'eid': eid, 'pid': region})
    data = ('\n'.join(lines + ['ENDDATA']) + '\n').encode('ascii')
    # This parser has no relation to either native extraction or Exodus decode.
    from pyNastran.bdf.bdf import read_bdf
    parsed = read_bdf(io.StringIO(data.decode('ascii')), validate=True, xref=True, debug=False)
    if {i: n.xyz.tolist() for i, n in parsed.nodes.items()} != nodes:
        raise ValueError('Independent BDF node readback differs')
    if len(parsed.elements) != len(mapping) or set(parsed.properties) != set(regions):
        raise ValueError('Independent BDF element/property inventory differs')
    for row in mapping:
        element = parsed.elements[row['eid']]
        key = (row['source_type'], row['source_id'])
        if element.node_ids != cells[key][1] or element.pid != row['pid'] or not element.Volume() > 0:
            raise ValueError('Independent BDF connectivity/property/volume check failed')
    return data, mapping


def seal_capture(capture_path, receipt_path, directory, *, schema_path):
    """Validate and seal six immutable members; never connect an application."""
    capture_path, receipt_path = Path(capture_path), Path(receipt_path)
    capture = strict_json(capture_path.read_bytes()); receipt = strict_json(receipt_path.read_bytes())
    sets = capture.get('sets', {'sidesets': [], 'nodesets': []})
    if sets != {'sidesets': [], 'nodesets': []}:
        raise ValueError('Sets are preserved by the label-catalog profile; this bundle profile cannot silently drop them')
    if not receipt.get('complete') or receipt.get('batch_nographics') is not True:
        raise ValueError('Completed batch acquisition receipt required')
    native_path = capture_path.parent / 'alternate-mesh.e'
    if (receipt['capture_sha256'] != digest(capture_path.read_bytes())
            or receipt['alternate_sha256'] != digest(native_path.read_bytes())):
        raise ValueError('Acquisition output changed after receipt')
    nodes, cells = check_capture(capture, native_path)
    if receipt['loaded_module_path'] != capture['exporter_module_path']:
        raise ValueError('Acquisition module identity mismatch')
    if receipt['candidate_sha256'] != receipt['candidate_after_sha256'] or digest(Path(capture['exporter_module_path']).read_bytes()) != receipt['candidate_sha256']:
        raise ValueError('Acquisition binary changed')
    regions = [x['id'] for x in capture['regions']]
    bdf, element_maps = _write_bdf(nodes, cells, regions)
    unit = capture['units']; source = receipt['collector_source_sha256']
    mesh = {'schema': 'cubit.mesh-snapshot.v1', 'source_length_unit': unit, 'coordinate_scale': 1.0,
            'units': unit, 'dimension': 3, 'order': 1,
            'nodes': [{'id': n, 'xyz': xyz, 'xyz_input': xyz} for n, xyz in sorted(nodes.items())],
            'elements': capture['elements'], 'regions': capture['regions'],
            'sets': {'sidesets': [], 'nodesets': []},
            'transform': {'mode': 'identity', 'origin_input': [0, 0, 0], 'origin_output': [0, 0, 0],
                          'basis_u': '+x', 'basis_v': '+y', 'plane_tolerance': 0,
                          'plane_tolerance_unit': unit, 'require_nonnegative_u': False},
            'options': {'sideset_mode': 'surface_elements', 'marker_cards': None,
                        'spc_components': None, 'selected_sideset_ids': []},
            'independent_boundary': {'acquisition': 'independent_geometry_query', 'orientation': 'outward',
                                     'faces': capture['faces'], 'edges': [], 'orientation_status': 'verified',
                                     'collector_source_sha256': source},
            'correspondence': {'nodes': [{'source_id': n, 'grid_id': n} for n in sorted(nodes)],
                               'elements': element_maps,
                               'regions': [{'source_id': n, 'pid': n, 'mid': n} for n in sorted(regions)],
                               'markers': []},
            'validation': {'structural_passed': True, 'independent_boundary_matched': 'oriented_match',
                           'ordering_certified': False, 'ordering_evidence_sha256': None,
                           'producer_self_report': True}}
    snapshot_fields = ('source_length_unit', 'coordinate_scale', 'units', 'dimension', 'order',
                       'nodes', 'elements', 'regions', 'sets', 'options', 'transform')
    snapshot = encoded({'schema': 'cubit.mesh-input-snapshot.v1', **{k: mesh[k] for k in snapshot_fields}})
    alternate = native_path.read_bytes()
    oracle = {'schema': 'cubit.mesh-ordering-oracle.v1', 'native_alternate_format': 'Exodus-classic-NetCDF',
              'alternate_export_sha256': digest(alternate), 'source_journal_sha256': receipt['source_journal_sha256'],
              'ccm_sha256': receipt['candidate_sha256'], 'native_source_digest': receipt['native_source_digest'],
              'node_identity_map_sha256': digest(encoded(mesh['correspondence']['nodes'])),
              'typed_ordered_connectivity_match': True, 'node_coordinate_match': True}
    evidence = {'schema': 'cubit.mesh-extraction-evidence.v1', 'execution_id': receipt['execution_id'],
                'started_utc': receipt['started_utc'], 'finished_utc': receipt['finished_utc'],
                'batch_nographics': True, 'cubit_version': capture['cubit_version'],
                'loaded_ccm_sha256': receipt['candidate_sha256'], 'collector_source_sha256': source,
                'input_snapshot_sha256': digest(snapshot), 'bdf_sha256': digest(bdf),
                'ordering_oracle_sha256': digest(encoded(oracle)), 'source_length_unit': unit,
                'coordinate_scale': 1.0, 'units': unit, 'native_source_digest': receipt['native_source_digest'],
                'pyd_used': False, 'pyd_sha256': None}
    members = {'mesh.json': encoded(mesh), 'mesh.bdf': bdf, 'input-snapshot.json': snapshot,
               'extraction-evidence.json': encoded(evidence), 'ordering-oracle.json': encoded(oracle),
               'alternate-mesh.e': alternate}
    manifest = {'schema': 'cubit.mesh-bundle.v1', 'units': unit, 'dimension': 3, 'order': 1,
                'files': [{'path': name, 'bytes': len(data), 'sha256': digest(data)} for name, data in sorted(members.items())],
                'provenance': {'exporter_version': 'neutral-native-candidate', 'source_commit': None,
                               'ccm_sha256': receipt['candidate_sha256'], 'pyd_sha256': None},
                'placeholder_physics': True, 'analysis_ready': False, 'downstream_identity_verified': False,
                'source_length_unit': unit, 'coordinate_scale': 1.0, 'input_snapshot_sha256': digest(snapshot),
                'extraction_execution_evidence_sha256': digest(encoded(evidence)), 'sideset_sid_map': []}
    from jsonschema import Draft202012Validator, FormatChecker
    schema = strict_json(Path(schema_path).read_bytes())
    for kind, obj in [('manifest', manifest), ('mesh', mesh), ('input_snapshot', strict_json(snapshot)),
                      ('extraction_evidence', evidence), ('ordering_oracle', oracle)]:
        Draft202012Validator({'$ref': '#/$defs/' + kind, '$defs': schema['$defs']}, format_checker=FormatChecker()).validate(obj)
    out = Path(directory)
    out.mkdir(exist_ok=False)
    for name, data in {**members, 'manifest.json': encoded(manifest)}.items():
        with (out / name).open('xb') as stream:
            stream.write(data)
    return out / 'manifest.json'


def label_catalog(capture_path, receipt_path):
    """Verify label membership independently; transfer no physical conditions.

    The catalog is distinct from the six-member no-set bundle. It does not
    weaken that consumer profile or assert downstream boundary assignment.
    """
    from scipy.io import netcdf_file
    capture_path, receipt_path = Path(capture_path), Path(receipt_path)
    capture = strict_json(capture_path.read_bytes()); receipt = strict_json(receipt_path.read_bytes())
    alternate = capture_path.parent / 'alternate-mesh.e'
    if (not receipt.get('complete') or receipt.get('batch_nographics') is not True
            or receipt['capture_sha256'] != digest(capture_path.read_bytes())
            or receipt['alternate_sha256'] != digest(alternate.read_bytes())):
        raise ValueError('Fresh hash-bound native capture required for labels')
    check_capture(capture, alternate)
    with netcdf_file(alternate, mode='r', mmap=False) as db:
        v = db.variables
        node_map = v['node_num_map'][:].tolist() if 'node_num_map' in v else list(range(1, int(db.dimensions['num_nodes']) + 1))
        def names(variable):
            if variable not in v:
                return []
            return [row.tobytes().split(b'\0', 1)[0].decode('utf-8') for row in v[variable][:]]
        block_ids = [int(i) for i in v['eb_prop1'][:].tolist()]
        blocks = dict(zip(block_ids, names('eb_names')))
        if blocks != {r['id']: r['name'] for r in capture['regions']}:
            raise ValueError('Native block labels differ')
        side_patterns = {
            'HEX8': ((0,1,5,4), (1,2,6,5), (2,3,7,6), (0,4,7,3), (0,3,2,1), (4,5,6,7)),
            'TET4': ((0,1,3), (1,2,3), (0,3,2), (0,2,1)),
        }
        native_cells = []
        for index in range(1, len(block_ids) + 1):
            c = v['connect' + str(index)]
            kind = 'HEX8' if c.shape[1] == 8 else 'TET4'
            native_cells.extend((kind, [int(node_map[int(n)-1]) for n in row]) for row in c[:].tolist())
        actual_sides = {}
        side_ids = [int(i) for i in v['ss_prop1'][:].tolist()] if 'ss_prop1' in v else []
        side_names = names('ss_names')
        if len(side_ids) != len(side_names):
            raise ValueError('Missing native sideset names')
        for index, (sid, name) in enumerate(zip(side_ids, side_names), 1):
            faces = []
            for eid, side in zip(v['elem_ss'+str(index)][:].tolist(), v['side_ss'+str(index)][:].tolist()):
                if not 1 <= int(eid) <= len(native_cells):
                    raise ValueError('Native sideset element outside mesh')
                kind, conn = native_cells[int(eid)-1]
                if not 1 <= int(side) <= len(side_patterns[kind]):
                    raise ValueError('Native sideset face outside cell')
                faces.append(tuple(sorted(conn[i] for i in side_patterns[kind][int(side)-1])))
            if len(set(faces)) != len(faces):
                raise ValueError('Duplicate native sideset face')
            actual_sides[sid] = (name, set(faces))
        desired_sides = {}
        for s in capture['sets']['sidesets']:
            faces = [tuple(sorted(f['node_ids'])) for f in s['faces']]
            nodes = sorted({n for f in s['faces'] for n in f['node_ids']})
            if s['node_ids'] != nodes or len(set(faces)) != len(faces) or s['id'] in desired_sides:
                raise ValueError('Invalid captured sideset membership')
            desired_sides[s['id']] = (s['name'], set(faces))
        if actual_sides != desired_sides:
            raise ValueError('Native sideset labels or face membership differ')
        node_ids = [int(i) for i in v['ns_prop1'][:].tolist()] if 'ns_prop1' in v else []
        node_names = names('ns_names')
        if len(node_ids) != len(node_names):
            raise ValueError('Missing native nodeset names')
        actual_nodes = {}
        for index, (sid, name) in enumerate(zip(node_ids, node_names), 1):
            indices = v['node_ns'+str(index)][:].tolist()
            if any(not 1 <= int(n) <= len(node_map) for n in indices):
                raise ValueError('Native nodeset index outside mesh')
            actual_nodes[sid] = (name, sorted(int(node_map[int(n)-1]) for n in indices))
        expected_nodes = {s['id']: (s['name'], s['node_ids']) for s in capture['sets']['nodesets']}
        if actual_nodes != expected_nodes:
            raise ValueError('Native nodeset labels or membership differ')
    return {'schema': 'cubit.mesh-label-catalog.v1', 'units': capture['units'],
            'native_capture_sha256': digest(capture_path.read_bytes()),
            'alternate_sha256': digest(alternate.read_bytes()),
            'regions': capture['regions'], 'sets': capture['sets'],
            'native_label_membership_verified': True, 'physical_conditions_transferred': False,
            'downstream_labels_verified': False, 'analysis_ready': False}


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--capture', type=Path, required=True)
    parser.add_argument('--receipt', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--schema', type=Path)
    parser.add_argument('--labels-only', action='store_true')
    args = parser.parse_args()
    if args.labels_only:
        data = encoded(label_catalog(args.capture, args.receipt))
        with args.output.open('xb') as stream:
            stream.write(data)
    else:
        if args.schema is None:
            parser.error('--schema is required for a sealed mesh bundle')
        seal_capture(args.capture, args.receipt, args.output, schema_path=args.schema)


if __name__ == '__main__':
    main()
