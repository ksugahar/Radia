"""Standalone vacuum split diagnostic; NEVER a three-engine lane pass."""
import json
import hashlib
import platform
import sys
import faulthandler
import importlib.metadata
import psutil
from pathlib import Path
import numpy as np
import ngsolve as ng
import radia as rad
import run_three_engine as driver
import run_finer_three as cache

root = Path(__file__).resolve().parent
output = root/'vacuum.json'

def main():
    if output.exists():
        raise FileExistsError(output)
    ng.SetNumThreads(8)
    environment = dict(host=platform.node(), python=sys.version, executable=sys.executable,
        versions={name:importlib.metadata.version(name) for name in ('radia','ngsolve','numpy','psutil')},
        radia_path=rad.__file__, ngsolve_path=ng.__file__,
        logical_cpus=psutil.cpu_count(), memory=psutil.virtual_memory()._asdict(),
        scripts={p.name:driver.sha256(p) for p in root.glob('*.py')})
    (root/'environment.json').write_text(json.dumps(environment,indent=2))
    process=psutil.Process()
    driver.solve_omega=cache.instrument(driver.solve_omega,cache.event,lambda:process.memory_info().rss)
    report = json.loads((root/'meshes/mesh_result.json').read_text())
    vol = root/'meshes/kelvin_domain.vol'
    assert report['passed']
    assert driver.sha256(vol) == report['artifacts']['kelvin_domain_vol_sha256']
    mesh = ng.Mesh(str(vol))
    assert driver.has_kelvin_identification(mesh)
    offset = np.asarray(driver.detect_kelvin_offset(mesh))
    center = tuple((np.asarray(report['kelvin_physical_center_m'])+offset).tolist())
    rad.UtiDelAll()
    coil, manifest = driver.build_coil()
    points = driver.observation_points()
    cache.kelvin.project_source_total_hodge = cache.cached_hodge
    cache.kelvin.project_source_interface_potential = cache.measured_interface
    # Exact prescribed coil B is the vacuum target, not another FE solution.
    target = driver.MU0 * np.asarray(rad.Fld(coil, 'h', points), dtype=float).reshape(-1,3)
    field, diagnostics = driver.solve_omega(
        mesh, coil, 1.0, nonlinear=False, order=2,
        nonlinear_tolerance=2e-5, nonlinear_maximum_iterations=80,
        nonlinear_verbose=False, kelvin_center=center,
        kelvin_radius=float(report['kelvin_radius_m']), points=points,
        source_trace_tolerance=0.05, source_projection_order=2,
        bonus_intorder=4)
    mask = abs(points[:,0]) <= 0.010+1e-12
    def metric(a,b):
        return float(np.linalg.norm(a-b)/np.linalg.norm(b))
    result = dict(schema='standalone-vacuum-omega-diagnostic/v1',
        three_engine_acceptance=False, host=platform.node(),
        python=sys.version, radia_path=rad.__file__, ngsolve_path=ng.__file__,
        threads=8, mu_r=1.0, source_order=2, response_order=2,
        mesh_sha256=driver.sha256(vol), coil=manifest,
        points_m=points.tolist(), target_B_T=target.tolist(),
        computed_B_T=np.asarray(field).tolist(),
        raw_all_relative_rms=metric(field,target),
        raw_core_relative_rms=metric(np.asarray(field)[mask],target[mask]),
        diagnostics=diagnostics,
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        interpretation='Raw vacuum source reproduction only. No magnetic-material performance or three-engine acceptance claim. No pass threshold inferred after seeing the result.')
    output.write_text(json.dumps(result,indent=2))
    print(json.dumps({k:result[k] for k in ['schema','raw_all_relative_rms','raw_core_relative_rms']}),flush=True)

if __name__ == '__main__':
    faulthandler.dump_traceback_later(1800,exit=True)
    try:
        main()
    except BaseException as exc:
        (root/'failure.json').write_text(json.dumps(dict(error=repr(exc),three_engine_acceptance=False)))
        raise
    finally:
        faulthandler.cancel_dump_traceback_later()
