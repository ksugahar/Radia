"""Standalone Gmsh sparse-tag field reproducer; no Radia inputs or imports.

Run each case in a fresh process with an explicit external memory limit.
The Gmsh MSH manual documents that post-processing storage scales with maxTag:
https://gmsh.info/dev/doc/texinfo/gmsh.html#MSH-file-format
An exception alone does not distinguish allocation failure from a native bug.
"""
import argparse
import importlib.metadata
import json
from pathlib import Path
import traceback

import gmsh


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--first-tag', type=int, default=2**31 + 17)
    parser.add_argument('--phase', choices=['mesh', 'field', 'read'], default='field')
    parser.add_argument('--mesh', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    tags = [args.first_tag, args.first_tag + 1]
    report = {'distribution_version': importlib.metadata.version('gmsh'),
              'first_tag': args.first_tag, 'phase': args.phase,
              'route': 'Gmsh API only; two nodes, one line',
              'passed': False, 'stage': 'initialize'}
    gmsh.initialize(['-noconfig'])
    try:
        report['native_version'] = gmsh.option.getString('General.Version')
        report['msh_version'] = gmsh.option.getNumber('Mesh.MshFileVersion')
        if args.phase == 'read':
            report['stage'] = 'read_mesh'
            gmsh.open(str(args.mesh))
        else:
            gmsh.model.add('sparse_field_reproducer')
            gmsh.model.addDiscreteEntity(1, 1)
            report['stage'] = 'add_nodes'
            gmsh.model.mesh.addNodes(1, 1, tags, [0, 0, 0, 1, 0, 0])
            report['stage'] = 'add_elements'
            gmsh.model.mesh.addElementsByType(1, 1, [1], tags)
            report['stage'] = 'write_mesh'
            gmsh.write(str(args.mesh))
            if args.phase == 'field':
                view = gmsh.view.add('prescribed_vector')
                report['stage'] = 'add_model_data'
                gmsh.view.addModelData(view, 0, 'sparse_field_reproducer',
                                       'NodeData', tags, [[1, -2, 3], [4, 5, -6]],
                                       numComponents=3)
                report['stage'] = 'read_model_data'
                _, actual_tags, values, _, _ = gmsh.view.getModelData(view, 0)
                report['field_values'] = dict(zip(map(int, actual_tags), map(list, values)))
        report['node_tags'] = list(map(int, gmsh.model.mesh.getNodes()[0]))
        report['passed'] = report['node_tags'] == tags and (
            args.phase != 'field' or report['field_values'] == {
                tags[0]: [1., -2., 3.], tags[1]: [4., 5., -6.]})
    except Exception as exc:
        report['error'] = repr(exc)
        report['traceback'] = traceback.format_exc()
    finally:
        gmsh.finalize()
        args.output.write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report))
    return 0 if report['passed'] else 1


if __name__ == '__main__':
    raise SystemExit(main())
