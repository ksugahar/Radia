"""Staged replay of the existing shape lane; never changes a live install.

Run prepare/evaluate on an idle compute host and mesh on licensed LAB/100.
The prepare phase records its staircase reference from that exact design run;
only checked Cubit meshes cross back to the compute host for evaluation.
An explicit topopt-cad source overlay permits testing a Python-only fix against
an identified wheel; both identities are recorded, never called a wheel release.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import platform
import sys
import uuid
from contextlib import contextmanager
from pathlib import Path
from types import SimpleNamespace


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


MESH_NAMES = {"tet_reference", "hex_coarse", "hex_fine"}


@contextmanager
def phase_receipt(out, phase):
    """A failed/repeated phase cannot inherit a prior success receipt."""
    state_path = out / f'{phase}.state.json'
    result_name = 'field.json' if phase == 'evaluate' else f'{phase}.json'
    if (out / result_name).exists():
        raise FileExistsError(f'Use a fresh phase directory: {out / result_name}')
    state = {'run_id': str(uuid.uuid4()), 'phase': phase, 'status': 'in_progress'}
    with state_path.open('x', encoding='utf-8') as stream:
        json.dump(state, stream)
    try:
        yield state
    except BaseException as exc:
        state.update(status='failed', error=f'{type(exc).__name__}: {exc}')
        raise
    else:
        state['status'] = 'completed'
    finally:
        state_path.write_text(json.dumps(state, indent=2), encoding='utf-8')


def verify_mesh_parent(out, preparation, meshes):
    if meshes['prepare_sha256'] != digest(out / 'prepare.json'):
        raise ValueError('Mesh receipt belongs to another preparation')
    if meshes['input_sha256'] != preparation['inputs']:
        raise ValueError('Mesh inputs differ from preparation')
    for field in ('mesh_results', 'checks', 'vol_sha256'):
        if set(meshes[field]) != MESH_NAMES:
            raise ValueError(f'Incomplete/unexpected mesh set: {field}')


def run_phase(args):
    out = args.directory.resolve()
    import radia
    if args.topopt_cad:
        load(args.topopt_cad.resolve(), "radia.topopt_cad")
    lane = load(Path(__file__).with_name("test_shape_regen_lane.py"), "shape_lane")
    import ngsolve
    import radia.topopt_cad as cad

    provenance = {
        "host": platform.node(), "python": sys.version,
        "interpreter": sys.executable, "radia_file": radia.__file__,
        "radia_version": radia.__version__, "ngsolve_version": ngsolve.__version__,
        "topopt_cad_file": cad.__file__, "topopt_cad_sha256": digest(cad.__file__),
        "lane_sha256": digest(lane.__file__), "driver_sha256": digest(__file__),
        "native_files": {p.name: digest(p) for p in Path(radia.__file__).parent.glob("*.pyd")},
    }
    print(json.dumps(provenance), flush=True)
    if args.phase == "prepare":
        regen = lane.prepare_shape(out)
        lane.test_grid_iso_stl_is_watertight_with_small_drift(regen)
        lane.write_levelset_exodus(regen.mesh, regen.nodal, out / "design_lsd.exo", level=0.5)
        data = {
            "run_id": args.run_id,
            "provenance": provenance, "stl_info": regen.stl_info,
            "staircase": lane.evaluate_staircase(regen),
            "coarse_stl_info": regen.coarse_stl_info, "timings": regen.timings,
            "inputs": {name: digest(out / name) for name in
                       ("design.stl", "design_300.stl", "design_lsd.exo")},
        }
        (out / "prepare.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    elif args.phase == "mesh":
        data = json.loads((out / "prepare.json").read_text(encoding="utf-8"))
        for name, expected in data["inputs"].items():
            assert digest(out / name) == expected, name
        regen = SimpleNamespace(out=out, stl=out / "design.stl",
                                coarse_stl=out / "design_300.stl", timings={})
        if args.phase == "mesh":
            from unittest.mock import patch
            from radia_mcp.cubit import session
            journal = session.run_headless_journal
            journal_records = []

            def checked_journal(commands, **kwargs):
                commands = list(commands)
                if args.tet_curve_interval and 'volume all scheme tetmesh' in commands:
                    commands.insert(commands.index('mesh volume all'),
                                    f'curve all interval {args.tet_curve_interval}')
                result = journal(commands, **kwargs,
                                 command_plugin_directory=args.command_plugin_directory)
                journal_records.append({'commands': commands, 'process': result})
                (out / 'journals.json').write_text(json.dumps(journal_records, indent=2), encoding='utf-8')
                if result.get('exit_code') != 0:
                    raise RuntimeError(f'Cubit did not exit cleanly: {result}')
                return result

            # Process-local instrumentation only; no daemon, install or user init is edited.
            with patch.object(session, 'run_headless_journal', checked_journal):
                lane.mesh_shape(regen)
                # Use the exact same console/plugin contract for the ISO operation.
                def iso_batch(_step, commands, timeout_s=900):
                    checked_journal(commands, timeout_s=timeout_s, working_directory=out)
                    return {'status': 'ok'}
                with patch('radia_mcp.cubit.server._run_batch', iso_batch):
                    iso_union = lane.validate_cubit_iso_union(out / 'design_lsd.exo', out)
            (out / "mesh-export.json").write_text(json.dumps({
                "mesh_results": regen.mesh_results, "timings": regen.timings}, indent=2), encoding="utf-8")
        lane.test_stl_to_vol_gates(regen)
        from radia_mcp.cubit.server import cubit_check_vol
        checks = {}
        for name in regen.mesh_results:
            checks[name] = json.loads(cubit_check_vol(
                str(out / f"{name}.vol"),
                threshold_pct=regen.mesh_results[name]["closure_tolerance"] * 100,
                report_json=str(out / f"{name}.check.json")))
            assert checks[name].get("passed") is True, checks[name]
        data = {"run_id": args.run_id, "provenance": provenance, "mesh_results": regen.mesh_results,
                "timings": regen.timings, "checks": checks,
                "iso_union": iso_union,
                "prepare_sha256": digest(out / 'prepare.json'),
                "input_sha256": data['inputs'],
                "vol_sha256": {n: digest(out / f"{n}.vol") for n in regen.mesh_results}}
        (out / "mesh.json").write_text(json.dumps(data, indent=2), encoding="utf-8")
    else:
        preparation = json.loads((out / "prepare.json").read_text(encoding="utf-8"))
        meshes = json.loads((out / "mesh.json").read_text(encoding="utf-8"))
        verify_mesh_parent(out, preparation, meshes)
        for name, expected in meshes["vol_sha256"].items():
            assert digest(out / f"{name}.vol") == expected, name
            assert meshes["checks"][name]["passed"] is True, name
        for name, expected in preparation["inputs"].items():
            assert digest(out / name) == expected, name
        assert preparation['provenance']['lane_sha256'] == provenance['lane_sha256']
        assert preparation['provenance']['native_files'] == provenance['native_files']
        assert preparation['provenance']['topopt_cad_sha256'] == provenance['topopt_cad_sha256']
        ngsolve.SetNumThreads(4)
        state, objective, _ = lane.shape_load_builders()
        regen = SimpleNamespace(out=out, state_builder=state, objective_builder=objective,
                                staircase=preparation['staircase'],
                                stl_info=preparation['stl_info'],
                                coarse_stl_info=preparation['coarse_stl_info'],
                                mesh_results=meshes['mesh_results'], timings=preparation['timings'])
        regen.timings.update(meshes["timings"])
        lane.RESULTS = out / "field.json"
        lane.test_functional_reevaluation_on_regenerated_body(regen)
        result = json.loads(lane.RESULTS.read_text(encoding="utf-8"))
        result.update(completed=True, provenance=provenance,
                      preparation=preparation, mesh_identity=meshes,
                      mesh_receipt_sha256=digest(out / 'mesh.json'),
                      run_id=args.run_id)
        lane.RESULTS.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(f"completed {args.phase}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('phase', choices=('prepare', 'mesh', 'evaluate'))
    parser.add_argument('--directory', type=Path, required=True)
    parser.add_argument('--topopt-cad', type=Path)
    parser.add_argument('--command-plugin-directory', type=Path)
    parser.add_argument('--tet-curve-interval', type=int, default=0,
                        help='Resolve every imported facet boundary curve; no hardcoded IDs')
    args = parser.parse_args()
    out = args.directory.resolve()
    out.mkdir(parents=True, exist_ok=True)
    with phase_receipt(out, args.phase) as receipt:
        args.run_id = receipt['run_id']
        if args.phase != 'prepare':
            parents = ('prepare',) if args.phase == 'mesh' else ('prepare', 'mesh')
            for parent in parents:
                state = json.loads((out / f'{parent}.state.json').read_text(encoding='utf-8'))
                if state['status'] != 'completed':
                    raise ValueError(f'Parent phase incomplete: {parent}')
        run_phase(args)


if __name__ == "__main__":
    main()
