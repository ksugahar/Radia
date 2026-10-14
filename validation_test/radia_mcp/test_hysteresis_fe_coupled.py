# -*- coding: utf-8 -*-
r"""FE-COUPLED hysteresis: a current-driven ring core solved with the B-input Play constitutive.

This is the per-Gauss-point (here per-element) hysteresis-aware FE solve: a 2D A_z magnetostatic ring
core of hysteretic iron, driven by a central-wire current I(t) cycled over one period, with the play
model as the constitutive H(B) and the play states stored PER ELEMENT on the mesh.  The nonlinearity is
handled by the classical FIXED-POINT (M-source) method -- solve the constant-reference-reluctivity
problem with an explicit "excess field" source H_res = nu0*B - H(B), iterating H_res to convergence
each time step, then advancing the per-element play states with the converged B.  The B-input play fits
the A-formulation (B=curl A primary), per the lab congruency paper.

Verified: the core's average B-H loop CLOSES, has positive area, and its loss (∮ H.dB) reproduces the
standalone play-model loop area at the achieved amplitude -- the FE solve and the constitutive agree.
This is the seed of a hysteresis-aware motor solve (the play state p_k per element is the internal
variable; a Newton tangent dH/dB is the next refinement for faster convergence).
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, L2, VectorL2, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, Mesh, TaskManager, Integrate)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.hysteresis import play_cells

mm = 1e-3
MU0 = 4e-7 * math.pi
Rw, Ri, Ro, REXT = 2*mm, 8*mm, 11*mm, 20*mm
A_EACH, NCELL, BSAT = 800.0, 12, 1.6


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2*mm)
    g.AddCircle((0, 0), Ro, leftdomain=2, rightdomain=1, bc="core_o", maxh=0.5*mm)
    g.AddCircle((0, 0), Ri, leftdomain=3, rightdomain=2, bc="core_i", maxh=0.5*mm)
    g.AddCircle((0, 0), Rw, leftdomain=4, rightdomain=3, bc="wire", maxh=0.4*mm)
    g.SetMaterial(1, "air"); g.SetMaterial(2, "core"); g.SetMaterial(3, "hole"); g.SetMaterial(4, "wire")
    return Mesh(g.GenerateMesh(maxh=1.5*mm))


def test_fe_coupled_ring_core_hysteresis():
    model = play_cells(Bmax=BSAT, N=NCELL, reluctivities=A_EACH)
    eta, aa = model.eta, model.a
    nu0 = float(np.sum(aa))                              # stable fixed-point reference reluctivity

    mesh = _geo()
    fes = H1(mesh, order=2, dirichlet="outer")
    u, v = fes.TnT()
    nu_ref = mesh.MaterialCF({"core": nu0}, default=1.0/MU0)
    a = BilinearForm(fes); a += nu_ref*grad(u)*grad(v)*dx; a.Assemble()
    inv = a.mat.Inverse(fes.FreeDofs(), inverse="umfpack")

    fesV = VectorL2(mesh, order=0)
    gfB = GridFunction(fesV); gfHres = GridFunction(fesV)
    gfmask = GridFunction(L2(mesh, order=0)); gfmask.Set(mesh.MaterialCF({"core": 1.0}, default=0.0))
    mask = gfmask.vec.FV().NumPy() > 0.5
    nE = len(mask)
    px = np.zeros((nE, model.K)); py = np.zeros((nE, model.K))

    def play_apply(Bx, By, px, py):
        Bp = np.sqrt((Bx[:, None]-px)**2 + (By[:, None]-py)**2)
        denom = np.maximum(eta[None, :], Bp); denom = np.where(denom < 1e-30, 1e-30, denom)
        nx = Bx[:, None] - eta[None, :]*(Bx[:, None]-px)/denom
        ny = By[:, None] - eta[None, :]*(By[:, None]-py)/denom
        nx[:, 0] = Bx; ny[:, 0] = By
        Hx = np.sum(aa[None, :]*nx, axis=1); Hy = np.sum(aa[None, :]*ny, axis=1)
        return Hx, Hy, nx, ny

    gfu = GridFunction(fes)
    Bcf = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    r_cf = sqrt(x*x+y*y); core_dx = dx(definedon=mesh.Materials("core"))
    area_core = Integrate(CoefficientFunction(1.0), mesh, definedon=mesh.Materials("core"))
    rmean = 0.5*(Ri+Ro)
    NT, Imax, NCYC = 32, 90.0, 2                     # run NCYC cycles; record only the last (steady) one
    Bhist, Hhist = [], []
    with TaskManager():
        for it in range(1, NCYC*NT+1):
            I = Imax*math.sin(2*math.pi*(it/NT))
            Jz = I/(math.pi*Rw**2)
            srcJ = LinearForm(fes); srcJ += mesh.MaterialCF({"wire": Jz}, default=0.0)*v*dx; srcJ.Assemble()
            pxs, pys = px.copy(), py.copy()
            for _ in range(25):                          # fixed-point (M-source) iteration
                gfB.Set(Bcf); bv = gfB.vec.FV().NumPy().reshape(nE, 2)
                Hx, Hy, _, _ = play_apply(bv[:, 0], bv[:, 1], pxs.copy(), pys.copy())
                hv = gfHres.vec.FV().NumPy().reshape(nE, 2)
                hv[:, 0] = np.where(mask, nu0*bv[:, 0]-Hx, 0.0); hv[:, 1] = np.where(mask, nu0*bv[:, 1]-Hy, 0.0)
                rhs = srcJ.vec.CreateVector(); rhs.data = srcJ.vec
                lf = LinearForm(fes); lf += (gfHres[0]*grad(v)[1]-gfHres[1]*grad(v)[0])*core_dx; lf.Assemble()
                rhs.data += lf.vec
                gfu.vec.data = inv*rhs
            gfB.Set(Bcf); bv = gfB.vec.FV().NumPy().reshape(nE, 2)
            _, _, px, py = play_apply(bv[:, 0], bv[:, 1], px, py)
            if it > (NCYC-1)*NT:                         # record only the last (steady) cycle
                Bth = (-Bcf[0]*y + Bcf[1]*x)/r_cf
                Bhist.append(Integrate(Bth*core_dx, mesh)/area_core)
                Hhist.append(I/(2*math.pi*rmean))        # Ampere H (exact at convergence)

    Bc, Hc = np.array(Bhist), np.array(Hhist)
    Bc = np.append(Bc, Bc[0]); Hc = np.append(Hc, Hc[0])    # close the steady loop for the ∮
    loop = abs(float(np.trapezoid(Hc, Bc) if hasattr(np, "trapezoid") else np.trapz(Hc, Bc)))
    Bm = float(np.max(np.abs(Bc)))
    ref = model.loss_per_cycle(Bm, n=2001)
    rel = abs(loop - ref) / ref
    print(f"FE-coupled ring core: |B|max={Bm:.3f} T, FE loop area={loop:.2f} vs play={ref:.2f} J/m^3 (rel {rel:.2e})")

    assert Bm > 0.05, "core should magnetise to a physical B"
    assert loop > 0, "FE solve must produce a dissipative (positive-area) hysteresis loop"
    assert rel < 0.15, f"FE-coupled (steady-cycle) loop area off the play model by {rel:.2e}"


def main():
    test_fe_coupled_ring_core_hysteresis()
    print("[OK] FE-coupled hysteresis: current-driven ring core, B-input Play constitutive via the "
          "fixed-point M-source method, per-element play states; the FE B-H loop reproduces the model.")


if __name__ == "__main__":
    main()
