"""Export the rotating-magnet eddy-current solve as a field-animation dataset.

Runs the full-FEM A-phi reference of ``rotating_magnet_eddy.py`` for ONE
period T (default 2 s, Yano's low-Rm case) and writes the solved conductor
current density J(x, t) together with the magnet pose used for the same
solver step.  Nothing is recomputed for display: J is the backward-Euler
eddy current ``-sigma/dt (dA_s + dA_r + grad dphi)`` of the reference solve,
sampled at fixed, identified points on the plate mid-plane.

Dataset (HDF5, schema ``radia.field-animation.v1``, SI units, index base 0):
    /time_s                 (n,)        physical time of each solver step
    /magnet/center_m        (n, 3)      magnet centre for that step
    /magnet/rotation        (n, 3, 3)   columns = magnet u, v, w axes (global)
    /magnet/size_m          (3,)        cube edge lengths
    /magnet/magnetization_A_per_m (3,)  in the magnet body frame
    /conductor/box_min_m, /conductor/box_max_m (3,)
    /samples/id             (m,)        stable identifiers
    /samples/point_m        (m, 3)
    /field/J_A_per_m2       (n, m, 3)   signed current density; step 0 is the
                                        initial state (no eddy current yet)
Attributes on / record the source script SHA-256, solver settings and the
Rm of the case.  The aggregate J_rms / force of the golden test remain the
numerical references; this file is the spatial record for visualisation.

Run:
    python rotating_magnet_eddy_animation_export.py --out <dir>/rotating_magnet.h5
"""
from __future__ import annotations

import argparse
import hashlib
import os
import time

import h5py
import numpy as np
from ngsolve import (BilinearForm, GridFunction, InnerProduct, LinearForm,
                     TaskManager, dx, grad)

import rotating_magnet_eddy as rme

SCHEMA = "radia.field-animation.v1"


def _sha256(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def sample_grid(nx, nz):
    """Fixed mid-plane grid inside the plate; id = i * nz + k."""
    xs = np.linspace(-0.0055, 0.0055, nx)
    zs = np.linspace(-0.0055, 0.0055, nz)
    pts = np.array([[x, -0.00025, z] for x in xs for z in zs])
    return np.arange(len(pts), dtype=np.int64), pts


def solve_and_sample(T=2.0, maxh=0.0005, ncfg=48, nx=25, nz=25):
    mesh, fes, fesA, fesP, fesS, tnt, Mbf, msig, Kterm = rme.build(maxh)
    A, phi, W, psi = tnt
    steps = np.linspace(0, rme.NSTEP - 1, ncfg).round().astype(int)
    As, _ = rme.project_source(fesS, steps)
    ndof = fes.ndof
    gAs = GridFunction(fesS)
    G = np.zeros((ndof, ncfg))
    for j in range(ncfg):
        gAs.vec.FV().NumPy()[:] = As[:, j]
        lf = LinearForm(fes)
        lf += rme.SIGMA * (InnerProduct(gAs, W) + InnerProduct(gAs, grad(psi))) * dx
        lf.Assemble()
        G[:, j] = lf.vec.FV().NumPy().copy()
    dt = T / ncfg
    with TaskManager():
        Abf = BilinearForm(fes, symmetric=True)
        Abf += (msig + dt * Kterm) * dx
        Abf.Assemble()
        Ainv = Abf.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")
    ids, pts = sample_grid(nx, nz)
    mpts = mesh(pts[:, 0], pts[:, 1], pts[:, 2])
    tmp = GridFunction(fes)
    out = GridFunction(fes)
    gfu = GridFunction(fes)
    gfo = GridFunction(fes)
    gAn = GridFunction(fesS)
    gAo = GridFunction(fesS)
    J = np.zeros((ncfg, len(ids), 3))
    u = np.zeros(ndof)
    for n in range(1, ncfg):
        tmp.vec.FV().NumPy()[:] = u
        out.vec.data = Mbf.mat * tmp.vec
        rhs = out.vec.FV().NumPy() - (G[:, n] - G[:, n - 1])
        tmp.vec.FV().NumPy()[:] = rhs
        out.vec.data = Ainv * tmp.vec
        u_new = out.vec.FV().NumPy().copy()
        # True residual of the linear step (policy: <= 1e-6 relative).
        tmp.vec.FV().NumPy()[:] = u_new
        out.vec.data = Abf.mat * tmp.vec
        free = np.array(list(fes.FreeDofs()), dtype=bool)
        res = np.linalg.norm((out.vec.FV().NumPy() - rhs)[free])
        if res > 1e-6 * (np.linalg.norm(rhs[free]) + 1e-300):
            raise RuntimeError(f"step {n}: linear residual {res:.3e} above limit")
        gfu.vec.FV().NumPy()[:] = u_new
        gfo.vec.FV().NumPy()[:] = u
        gAn.vec.FV().NumPy()[:] = As[:, n]
        gAo.vec.FV().NumPy()[:] = As[:, n - 1]
        Ar, Ph = gfu.components
        Ar_o, Ph_o = gfo.components
        d = (gAn - gAo) + (Ar - Ar_o) + (grad(Ph) - grad(Ph_o))
        J[n] = ((-rme.SIGMA / dt) * d)(mpts)
        u = u_new
    poses = [rme._pose(int(s)) for s in steps]
    center = np.array([p["origin"] for p in poses])
    rot = np.array([np.column_stack([p["u_axis"], p["v_axis"], p["w_axis"]]) for p in poses])
    return dict(time_s=np.arange(ncfg) * dt, center=center, rot=rot, ids=ids, pts=pts, J=J,
                steps=steps, ndof=int(ndof), dt=dt,
                Rm=float(rme.MU0 * rme.SIGMA * (2 * np.pi * 2 / T) * rme.L_CHAR ** 2))


def write_h5(path, d, T, maxh, ncfg):
    with h5py.File(path, "w") as f:
        f.attrs["schema"] = SCHEMA
        f.attrs["units"] = "SI"
        f.attrs["index_base"] = 0
        f.attrs["case"] = "rotating_magnet_over_copper_plate"
        f.attrs["source_script"] = "docs/maglev/demos/rotating_magnet_eddy.py"
        f.attrs["source_script_sha256"] = _sha256(rme.__file__)
        f.attrs["export_script_sha256"] = _sha256(__file__)
        f.attrs["period_T_s"] = T
        f.attrs["maxh_m"] = maxh
        f.attrs["ncfg"] = ncfg
        f.attrs["ndof"] = d["ndof"]
        f.attrs["magnetic_reynolds_number"] = d["Rm"]
        f.attrs["field_name"] = "J_eddy"
        f.attrs["motion"] = "magnet pose changes (translation + rotation); magnetization fixed in body frame"
        f["time_s"] = d["time_s"]
        f["solver_step"] = d["steps"]
        g = f.create_group("magnet")
        g["center_m"] = d["center"]
        g["rotation"] = d["rot"]
        g["size_m"] = np.full(3, rme.MAG_MM * 1e-3)
        g["magnetization_A_per_m"] = np.array([0.0, rme.BR / rme.MU0, 0.0])
        c = f.create_group("conductor")
        c["box_min_m"] = np.array([-0.006, -0.0005, -0.006])
        c["box_max_m"] = np.array([0.006, 0.0, 0.006])
        c.attrs["conductivity_S_per_m"] = rme.SIGMA
        s = f.create_group("samples")
        s["id"] = d["ids"]
        s["point_m"] = d["pts"]
        f.create_group("field")["J_A_per_m2"] = d["J"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--period", type=float, default=2.0)
    ap.add_argument("--maxh", type=float, default=0.0005)
    ap.add_argument("--ncfg", type=int, default=48)
    ap.add_argument("--grid", type=int, default=25)
    a = ap.parse_args()
    t0 = time.time()
    d = solve_and_sample(a.period, a.maxh, a.ncfg, a.grid, a.grid)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    write_h5(a.out, d, a.period, a.maxh, a.ncfg)
    jr = np.sqrt((d["J"] ** 2).sum(-1).mean(-1))
    print(f"wrote {a.out}: {len(d['time_s'])} steps x {len(d['ids'])} samples, "
          f"Rm={d['Rm']:.3g}, sampled J_rms peak={jr.max():.1f} A/m^2, {time.time()-t0:.0f} s")


if __name__ == "__main__":
    main()
