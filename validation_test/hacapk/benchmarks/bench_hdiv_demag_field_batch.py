"""Scaling of the HDiv-VIM demag-field-at-quad-points (assemble_demag_field -> C++ _hdiv_demag_field_batch),
to size step 3.  Separates the C++ batch cost (the thing an H-matrix would accelerate) from the Python
packing cost (per-face sigma fit).  obs = element centroids = co-located with the charges (the real
self-eval regime, where the embed H-matrix hits the checkerboard 4% wall).

Run on mdx (idle).  PYTHONPATH -> C:\\temp\\radiabench\\src.  Writes results_demagbatch.json.
"""
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
import json, time, platform
from datetime import datetime
import numpy as np
import radia as rad
import radia._radia_pybind as rp
import ngsolve as ng
from ngsolve import TaskManager, SetNumThreads
from netgen.occ import Sphere, Pnt, OCCGeometry
import radia.vim._field as vf


def centroids(mesh):
    pts = []
    for el in mesh.Elements(ng.VOL):
        c = np.array([mesh[v].point for v in el.vertices]).mean(0)
        pts.append([c[0], c[1], c[2]])
    return np.asarray(pts, float)


def pack(mesh, gfM):
    """Reproduce assemble_demag_field's packing (vol_flat, surf_flat) so the C++ batch can be timed alone."""
    divM = ng.div(gfM)
    nrm = ng.specialcf.normal(mesh.dim)
    vol_flat = []
    for i in range(mesh.GetNE(ng.VOL)):
        ei = ng.ElementId(ng.VOL, i)
        V = np.array([mesh[v].point for v in mesh[ei].vertices], float)
        rho0, g = vf._fit_rho_linear(divM, mesh, V)
        vol_flat += V.ravel().tolist() + [float(rho0)] + g.tolist()
    surf_flat = []
    for i in range(mesh.GetNE(ng.BND)):
        ei = ng.ElementId(ng.BND, i)
        P = np.array([mesh[v].point for v in mesh[ei].vertices], float)
        ngeom = np.cross(P[1] - P[0], P[2] - P[0]); ngeom = ngeom / np.linalg.norm(ngeom)
        sig_fn = lambda p: float(np.array([float(gfM[k](mesh(p[0], p[1], p[2]))) for k in range(3)]) @ ngeom)
        sigma0, s, S = vf._fit_sigma_quadratic(sig_fn, P)
        surf_flat += P.ravel().tolist() + [float(sigma0)] + s.tolist() + S.ravel().tolist()
    return vol_flat, surf_flat


def run(maxh, nthreads):
    SetNumThreads(nthreads)
    with TaskManager():
        geo = OCCGeometry(Sphere(Pnt(0, 0, 0), 1.0))
        mesh = ng.Mesh(geo.GenerateMesh(maxh=maxh))
        fes = ng.HDiv(mesh, order=1)
        gfM = ng.GridFunction(fes)
        gfM.Set(ng.CF((0.0, 0.0, 1.0e6)))
        NE = mesh.GetNE(ng.VOL); NB = mesh.GetNE(ng.BND)
        pts = centroids(mesh)
        Nobs = pts.shape[0]
        # total assemble_demag_field
        t0 = time.perf_counter()
        _ = vf.assemble_demag_field(mesh, gfM, pts, quantity="h")
        t_total = time.perf_counter() - t0
        # isolate: packing vs C++ batch
        t0 = time.perf_counter()
        vol_flat, surf_flat = pack(mesh, gfM)
        t_pack = time.perf_counter() - t0
        obs_flat = pts.ravel().tolist()
        t0 = time.perf_counter()
        _ = rp._hdiv_demag_field_batch(vol_flat, surf_flat, obs_flat)
        t_cpp = time.perf_counter() - t0
    Nsrc = NE + NB
    return dict(maxh=maxh, NE=NE, NB=NB, N_src=Nsrc, N_obs=Nobs, nthreads=nthreads,
                pairs=Nobs * Nsrc, t_total=t_total, t_pack=t_pack, t_cpp_batch=t_cpp,
                cpp_ns_per_pair=1e9 * t_cpp / (Nobs * Nsrc))


def main():
    cores = os.cpu_count()
    results = []
    for maxh in [0.5, 0.38, 0.30, 0.24, 0.20]:
        r = run(maxh, cores)
        results.append(r)
        print("[demag] maxh=%.2f NE=%5d NB=%4d N_obs=%5d pairs=%.2e  t_total=%.2f  t_pack=%.2f  "
              "t_cpp=%.3f (%.0f ns/pair)"
              % (maxh, r["NE"], r["NB"], r["N_obs"], r["pairs"], r["t_total"], r["t_pack"],
                 r["t_cpp_batch"], r["cpp_ns_per_pair"]), flush=True)
    # C++ batch thread-scaling at the largest size
    big = max([0.20], key=lambda h: 1)
    scale = []
    for nt in sorted({t for t in [1, 2, 4, 8, cores] if t <= cores}):
        r = run(0.24, nt)
        scale.append(dict(nthreads=nt, t_cpp_batch=r["t_cpp_batch"], N_src=r["N_src"], N_obs=r["N_obs"]))
        print("[demag-scale] threads=%2d  t_cpp=%.3f" % (nt, r["t_cpp_batch"]), flush=True)
    out = dict(timestamp=datetime.now().isoformat(), hostname=platform.node(),
               benchmark="hdiv_demag_field_batch_scaling",
               problem=dict(cores=cores, geometry="unit sphere tet, HDiv order 1, uniform Mz, obs=centroids",
                            radia_version=rad.__version__, python_version=platform.python_version(),
                            platform=platform.platform()),
               results=results, cpp_thread_scaling=scale)
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results_demagbatch.json")
    with open(fn, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print("Results saved to", fn, flush=True)


if __name__ == "__main__":
    main()
