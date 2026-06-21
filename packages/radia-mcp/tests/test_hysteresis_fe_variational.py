# -*- coding: utf-8 -*-
r"""VARIATIONAL (energy-line-search) NEWTON for the FE-coupled B-input Play hysteresis.

The companion to tests/test_hysteresis_fe_coupled.py (which uses the robust FIXED-POINT M-source
method).  Here the same current-driven ring core is solved by a Newton-Raphson that is made GLOBALLY
CONVERGENT by a variational principle: the B-input Play field H(B) is the gradient of a CONVEX
incremental co-energy density psi*(B) (hysteresis.play_coenergy_density), so the discrete energy

    Pi(A) = INT psi*(B(A)) dx - INT J.A dx

is convex in the nodal potential A (psi* convex in B, B = curl A linear in A).  Its unique minimiser is
the FE solution, and a Newton step delta = -K^{-1} R (K = INT (dH/dB) curl.curl, the consistent SPD
tangent from play_tangent; R the residual = grad(Pi)) damped by an ARMIJO LINE SEARCH ON Pi -- not on
the residual norm -- is a guaranteed descent direction (grad(Pi).delta = -R^T K^{-1} R < 0), so it
always finds a decreasing step and converges.  This is the rate-independent-plasticity structure
(p_prev <-> back-stress, eta_k <-> yield radius, the play update <-> radial return; the incremental
energy is convex and Newton on it with the consistent tangent converges).

order-1 H1 makes B = curl A exactly piecewise constant per element, so the per-element play state is
exact; the element-average B is taken with Integrate(., element_wise) (the gfB.Set L2-projection is
inexact at order 1, which -- not an intrinsic kink pathology -- was what made an earlier Newton attempt
diverge).  Verified here: (1) Newton converges every time step in a handful of iterations (vs the
fixed-point's ~25 M-source sweeps/step), (2) the energy Pi decreases MONOTONICALLY within each solve
(the convexity signature), and (3) the steady B-H loop area matches the standalone play model TIGHTER
than the fixed-point solve does, because each step is driven to a small residual.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (H1, L2, BilinearForm, LinearForm, GridFunction, CoefficientFunction,
                     grad, dx, x, y, sqrt, Mesh, TaskManager, Integrate, InnerProduct, Projector, Norm)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.hysteresis import play_cells

mm = 1e-3
MU0 = 4e-7 * math.pi
Rw, Ri, Ro, REXT = 2*mm, 8*mm, 11*mm, 20*mm
A_EACH, NCELL, BSAT = 800.0, 12, 1.6
NU_AIR = 1.0/MU0


def _geo():
    g = SplineGeometry()
    g.AddCircle((0, 0), REXT, leftdomain=1, rightdomain=0, bc="outer", maxh=2*mm)
    g.AddCircle((0, 0), Ro, leftdomain=2, rightdomain=1, bc="core_o", maxh=0.5*mm)
    g.AddCircle((0, 0), Ri, leftdomain=3, rightdomain=2, bc="core_i", maxh=0.5*mm)
    g.AddCircle((0, 0), Rw, leftdomain=4, rightdomain=3, bc="wire", maxh=0.4*mm)
    g.SetMaterial(1, "air"); g.SetMaterial(2, "core"); g.SetMaterial(3, "hole"); g.SetMaterial(4, "wire")
    return Mesh(g.GenerateMesh(maxh=1.5*mm))


def test_fe_variational_newton_ring_core():
    model = play_cells(Bmax=BSAT, N=NCELL, reluctivities=A_EACH)
    eta, aa = model.eta, model.a

    mesh = _geo()
    fes = H1(mesh, order=1, dirichlet="outer")
    u, v = fes.TnT()
    gfu = GridFunction(fes)
    Bcf = CoefficientFunction((grad(gfu)[1], -grad(gfu)[0]))
    curlu = CoefficientFunction((grad(u)[1], -grad(u)[0]))
    curlv = CoefficientFunction((grad(v)[1], -grad(v)[0]))

    fesS = L2(mesh, order=0)
    gfHx = GridFunction(fesS); gfHy = GridFunction(fesS)
    gfTxx = GridFunction(fesS); gfTxy = GridFunction(fesS); gfTyy = GridFunction(fesS); gfPsi = GridFunction(fesS)
    gfmask = GridFunction(fesS); gfmask.Set(mesh.MaterialCF({"core": 1.0}, default=0.0))
    mask = gfmask.vec.FV().NumPy() > 0.5
    nE = len(mask); ci = np.where(mask)[0]
    areas = Integrate(CoefficientFunction(1.0), mesh, element_wise=True).NumPy()
    px = np.zeros((nE, model.K)); py = np.zeros((nE, model.K))
    proj = Projector(fes.FreeDofs(), True)
    Tmat = CoefficientFunction((gfTxx, gfTxy, gfTxy, gfTyy), dims=(2, 2))

    srcJ = LinearForm(fes); srcJ += mesh.MaterialCF({"wire": 1.0}, default=0.0)*v*dx; srcJ.Assemble()
    Rform = LinearForm(fes); Rform += (gfHx*grad(v)[1]-gfHy*grad(v)[0])*dx
    Kform = BilinearForm(fes); Kform += (Tmat*curlu)*curlv*dx

    def elem_B(dofs):
        gfu.vec.data = dofs
        bx = Integrate(Bcf[0], mesh, element_wise=True).NumPy()/areas
        by = Integrate(Bcf[1], mesh, element_wise=True).NumPy()/areas
        return bx, by

    def play_apply(Bx, By, pxl, pyl):
        Bp = np.sqrt((Bx[:, None]-pxl)**2 + (By[:, None]-pyl)**2)
        denom = np.maximum(eta[None, :], Bp); denom = np.where(denom < 1e-30, 1e-30, denom)
        nx = Bx[:, None] - eta[None, :]*(Bx[:, None]-pxl)/denom
        ny = By[:, None] - eta[None, :]*(By[:, None]-pyl)/denom
        nx[:, 0] = Bx; ny[:, 0] = By
        return np.sum(aa[None, :]*nx, axis=1), np.sum(aa[None, :]*ny, axis=1), nx, ny

    def constitutive(bx, by):
        Hx = NU_AIR*bx.copy(); Hy = NU_AIR*by.copy()
        Txx = np.full(nE, NU_AIR); Tyy = np.full(nE, NU_AIR); Txy = np.zeros(nE)
        psi = 0.5*NU_AIR*(bx*bx + by*by)
        Hxc, Hyc, _, _ = play_apply(bx[ci], by[ci], px[ci].copy(), py[ci].copy())
        txx, txy, tyy = model.play_tangent(bx[ci], by[ci], px[ci], py[ci])
        psic = model.play_coenergy_density(bx[ci], by[ci], px[ci], py[ci])
        Hx[ci], Hy[ci] = Hxc, Hyc; Txx[ci], Txy[ci], Tyy[ci] = txx, txy, tyy; psi[ci] = psic
        return Hx, Hy, Txx, Txy, Tyy, psi

    def set_fields(bx, by):
        Hx, Hy, Txx, Txy, Tyy, psi = constitutive(bx, by)
        for gf, val in ((gfHx, Hx), (gfHy, Hy), (gfTxx, Txx), (gfTxy, Txy), (gfTyy, Tyy), (gfPsi, psi)):
            gf.vec.FV().NumPy()[:] = val

    def energy(dofs, Jz):
        bx, by = elem_B(dofs)
        _, _, _, _, _, psi = constitutive(bx, by)
        gfPsi.vec.FV().NumPy()[:] = psi
        return Integrate(gfPsi, mesh) - Jz*InnerProduct(srcJ.vec, dofs)

    r_cf = sqrt(x*x+y*y)
    core_dx = dx(definedon=mesh.Materials("core"))
    area_core = Integrate(CoefficientFunction(1.0), mesh, definedon=mesh.Materials("core"))
    rmean = 0.5*(Ri+Ro)
    NT, Imax, NCYC = 32, 90.0, 2
    Bhist, Hhist = [], []
    max_newton, worst_mono = 0, 0.0

    with TaskManager():
        for it in range(1, NCYC*NT+1):
            I = Imax*math.sin(2*math.pi*(it/NT)); Jz = I/(math.pi*Rw**2)
            E_seq, niter = [], 0
            for newton in range(40):
                bx, by = elem_B(gfu.vec); set_fields(bx, by)
                Rform.Assemble()
                Rvec = Rform.vec.CreateVector(); Rvec.data = Rform.vec - Jz*srcJ.vec
                Rf = Rvec.CreateVector(); Rf.data = proj*Rvec
                if Norm(Rf) < 1e-10*max(abs(Jz), 1.0):
                    niter = newton; break
                Kform.Assemble()
                Kinv = Kform.mat.Inverse(fes.FreeDofs(), inverse="umfpack")
                delta = gfu.vec.CreateVector(); delta.data = -(Kinv*Rvec)
                dderiv = InnerProduct(Rf, delta)
                assert dderiv < 1e-12, "Newton direction must be a descent direction for the convex energy"
                E0 = energy(gfu.vec, Jz); E_seq.append(E0)
                base = gfu.vec.CreateVector(); base.data = gfu.vec
                alpha = 1.0
                while alpha > 1e-6:
                    trial = gfu.vec.CreateVector(); trial.data = base + alpha*delta
                    if energy(trial, Jz) <= E0 + 1e-4*alpha*dderiv:
                        break
                    alpha *= 0.5
                gfu.vec.data = base + alpha*delta
                niter = newton+1
            max_newton = max(max_newton, niter)
            # energy strictly non-increasing within the solve (convexity signature)
            for k in range(1, len(E_seq)):
                worst_mono = max(worst_mono, E_seq[k] - E_seq[k-1])
            # commit play states with the converged B
            bx, by = elem_B(gfu.vec)
            _, _, nxc, nyc = play_apply(bx[ci], by[ci], px[ci], py[ci])
            px[ci], py[ci] = nxc, nyc
            if it > (NCYC-1)*NT:
                Bth = (-Bcf[0]*y + Bcf[1]*x)/r_cf
                Bhist.append(Integrate(Bth*core_dx, mesh)/area_core)
                Hhist.append(I/(2*math.pi*rmean))

    Bc, Hc = np.array(Bhist), np.array(Hhist)
    Bc = np.append(Bc, Bc[0]); Hc = np.append(Hc, Hc[0])
    loop = abs(float(np.trapezoid(Hc, Bc)))
    Bm = float(np.max(np.abs(Bc)))
    ref = model.loss_per_cycle(Bm, n=2001)
    rel = abs(loop-ref)/ref
    print(f"VARIATIONAL NEWTON ring core: max Newton iters/step={max_newton}, worst energy increase "
          f"within a solve={worst_mono:.2e} J; |B|max={Bm:.3f} T, FE loop area={loop:.2f} vs "
          f"play={ref:.2f} J/m^3 (rel {rel:.2e})")

    assert Bm > 0.05, "core should magnetise to a physical B"
    assert max_newton <= 12, ("the variational (energy-line-search) Newton must converge in a handful "
                              "of iterations every step -- global convergence")
    assert worst_mono <= 1e-9 * max(abs(loop), 1.0), ("the energy must decrease monotonically within "
                                                      "each solve (the convex-energy signature)")
    assert loop > 0, "the FE solve must produce a dissipative (positive-area) hysteresis loop"
    assert rel < 0.05, ("Newton drives each step to a small residual, so the steady loop should match "
                        "the play model tighter than the fixed-point solve's ~11%")


def main():
    test_fe_variational_newton_ring_core()
    print("[OK] variational (energy-line-search) Newton FE-coupled hysteresis: convex co-energy merit "
          "function + consistent SPD play tangent -> globally convergent, a few iters/step, loop area "
          "matches the play model tighter than the fixed-point method.")


if __name__ == "__main__":
    main()
