"""Fresh-process P1 construction cost on a self-authored thin ring.

Run once per source/size/thread combination on an admitted idle host. The
expected route is verified, not substituted. Native and Python identities
are explicit; output contains no machine names or paths. Native storage
statistics are payload estimates; the process peak also includes allocator,
quadrature, sparse-factor and runtime storage. This is a construction-only
measurement, not the complete solve peak.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path
import ngsolve as ng
import numpy as np
from radia import _radia_pybind as native
from radia.bem_sibc_solver import ScalarBIESIBCSolver
from run_loop_work_ring import ring_mesh
from run_loop_work_performance import memory


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--level',type=int,required=True)
    parser.add_argument('--threads',type=int,default=1)
    parser.add_argument('--expected-route',choices=['dense-entry','p1-entry-on-demand'],required=True)
    parser.add_argument('--source-commit',required=True)
    parser.add_argument('--cache-bytes',type=int,default=8*1024*1024)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    if args.threads<1 or args.cache_bytes<0:
        parser.error('Threads must be positive and cache bytes nonnegative')
    if len(args.source_commit)!=40 or any(c not in '0123456789abcdef' for c in args.source_commit):
        parser.error('Source identity requires an exact commit')
    binary=Path(native.__file__)
    provenance=json.loads(Path(str(binary)+".build.json").read_text(encoding="utf-8"))
    binary_hash=hashlib.sha256(binary.read_bytes()).hexdigest()
    if (provenance.get("source_dirty") or provenance.get("source_commit")!=args.source_commit
            or provenance.get("binary_sha256")!=binary_hash):
        raise RuntimeError("Clean exact-source native provenance is required for cost evidence")
    import radia.bem_sibc_solver as body_module
    source_root=Path(body_module.__file__).resolve().parents[2]
    runtime_head=subprocess.check_output(["git","-C",str(source_root),"rev-parse","HEAD"],text=True).strip()
    modified=subprocess.check_output(["git","-C",str(source_root),"diff","HEAD","--name-only","--","src/radia/bem_sibc_solver.py"],text=True).strip()
    if runtime_head!=args.source_commit or modified:
        raise RuntimeError("Solver source differs from the declared clean revision")
    ng.SetNumThreads(args.threads)
    phases=[]
    def phase(tag,message):
        phases.append(dict(message=message,seconds=time.perf_counter()-started,**memory()))
    with ng.TaskManager():
        started=time.perf_counter()
        mesh,points,triangles=ring_mesh(args.level)
        controls=dict(order=1,assemble_dense=False,use_intree_hacapk=True,
            use_intree_bem=True,intree_geom_order=1,intree_singular_n_q=6,
            intree_regular_quad_degree=7,hacapk_aca_eps=1e-10,loop_work_backend='fmm',log_fn=phase)
        if args.expected_route=='p1-entry-on-demand':
            controls['hacapk_entry_cache_bytes']=args.cache_bytes
        solver=ScalarBIESIBCSolver(mesh,**controls)
        elapsed=time.perf_counter()-started
        route=solver._SL_hacapk.GetStats().get('construction_route','dense-entry')
        if route!=args.expected_route:
            raise RuntimeError('Actual construction route differs from requested measurement')
        result=dict(schema='radia.p1-construction-cost.v1',source_commit=args.source_commit,
            native_sha256=binary_hash,
            measurement_driver_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            solver_source_sha256=hashlib.sha256(Path(__import__('radia.bem_sibc_solver',fromlist=['__file__']).__file__).read_bytes()).hexdigest(),
            os=platform.system(),ngsolve_version=ng.__version__,numpy_version=np.__version__,
            faces=len(triangles),vertices=len(points),threads_requested=args.threads,
            native_pool_max_threads=native.HLUMaxThreads(),construction_route=route,
            regular_degree=7,singular_order=6,aca_eps=1e-10,leaf=64,eta=2.,
            requested_cache_bytes=args.cache_bytes if route=='p1-entry-on-demand' else None,
            shared_provider=solver._entry_provider.GetStats() if hasattr(solver,'_entry_provider') else None,
            build_diagnostics=getattr(solver,'_body_build_diagnostics',None),
            SL=solver._SL_hacapk.GetStats(),DL=solver._DL_hacapk.GetStats(),
            geometry_sha256=hashlib.sha256(points.tobytes()+triangles.tobytes()).hexdigest(),
            construction_seconds=elapsed,phases=phases,**memory())
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(result),flush=True)


if __name__=='__main__':
    main()
