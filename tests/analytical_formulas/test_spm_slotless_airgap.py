"""Guards: slotless SPM air-gap field (Zhu-Howe 1993) in radia.analytical_formulas.

Assumes: NGSolve.  Oracle is an independent 2-D A_z FEM of the same idealised
machine (infinitely permeable iron = natural boundary, radial square-wave
magnetisation); the closed form must agree within the FEM's own discretisation
error at this level (~0.1-0.2 % at order 4).
"""
import numpy as np
import pytest

from radia.analytical_formulas.spm_slotless_airgap import MU0, slotless_spm_airgap_field

ngsolve = pytest.importorskip("ngsolve")


def _fem_gap_field(p, Br, mur, Rr, Rm, Rs, ap, r, th, order=4, maxh=7e-4):
    from ngsolve import (BilinearForm, CF, GridFunction, H1, IfPos, LinearForm, Mesh,
                         atan2, cos, dx, grad, sqrt, x, y)
    import netgen.occ as occ
    gap = occ.Circle((0, 0), Rs).Face() - occ.Circle((0, 0), Rm).Face(); gap.faces.name = "gap"
    ring = occ.Circle((0, 0), Rm).Face() - occ.Circle((0, 0), Rr).Face(); ring.faces.name = "ring"
    mesh = Mesh(occ.OCCGeometry(occ.Glue([gap, ring]), dim=2).GenerateMesh(maxh=maxh)); mesh.Curve(order)
    phase = atan2(y, x) * p; c0 = float(np.cos(ap * np.pi / 2)); rr = sqrt(x * x + y * y)
    s = IfPos(cos(phase) - c0, 1, IfPos(-cos(phase) - c0, -1, 0))
    M = CF((Br / MU0 * s * x / rr, Br / MU0 * s * y / rr))
    fes = H1(mesh, order=order); u, v = fes.TnT()
    nu = mesh.MaterialCF({"ring": 1 / (MU0 * mur)}, default=1 / MU0)
    a = BilinearForm(nu * grad(u) * grad(v) * dx + 1e-12 * u * v * dx).Assemble()
    f = LinearForm(M / mur * CF((grad(v)[1], -grad(v)[0])) * dx(definedon=mesh.Materials("ring"))).Assemble()
    gA = GridFunction(fes); gA.vec.data = a.mat.Inverse(inverse="sparsecholesky") * f.vec
    B = CF((grad(gA)[1], -grad(gA)[0]))
    vals = np.array([B(mesh(r * np.cos(t), r * np.sin(t))) for t in th])
    return (vals[:, 0] * np.cos(th) + vals[:, 1] * np.sin(th),
            -vals[:, 0] * np.sin(th) + vals[:, 1] * np.cos(th))


def test_matches_independent_fem():
    p, Br, mur, Rr, Rm, Rs, ap = 2, 1.2, 1.05, 0.025, 0.03, 0.032, 0.8
    r = 0.5 * (Rm + Rs); th = np.linspace(0, np.pi / p, 41)[1:-1]
    br, bt = slotless_spm_airgap_field(r, th, p, Br, mur, Rr, Rm, Rs, ap)
    with ngsolve.TaskManager():          # caller-owned parallel region; the helper stays region-free
        fr, ft = _fem_gap_field(p, Br, mur, Rr, Rm, Rs, ap, r, th)
    sc = np.abs(br).max()
    err_r, err_t = np.abs(fr - br).max() / sc, np.abs(ft - bt).max() / sc
    assert err_r < 5e-3 and err_t < 5e-3, f"B_r err {err_r:.2e}, B_theta err {err_t:.2e} (of max|B_r| {sc:.4f} T)"


def test_iron_boundary_and_symmetry():
    p, Br, mur, Rr, Rm, Rs, ap = 4, 1.0, 1.0, 0.04, 0.044, 0.0455, 0.6
    th = np.linspace(-np.pi, np.pi, 73)
    br, bt = slotless_spm_airgap_field(Rs, th, p, Br, mur, Rr, Rm, Rs, ap)
    assert np.abs(bt).max() <= 1e-12 * np.abs(br).max()          # H_t = 0 on the stator iron
    br2, _ = slotless_spm_airgap_field(Rs, th + np.pi / p, p, Br, mur, Rr, Rm, Rs, ap)
    assert np.allclose(br2, -br, atol=1e-12)                        # half-period anti-symmetry


def test_rejects_unsupported_input():
    with pytest.raises(ValueError, match="pole_pairs"):
        slotless_spm_airgap_field(0.031, [0.0], 1, 1.2, 1.05, 0.025, 0.03, 0.032, 0.8)
    with pytest.raises(ValueError, match="integer"):
        slotless_spm_airgap_field(0.031, [0.0], 2.9, 1.2, 1.05, 0.025, 0.03, 0.032, 0.8)
    with pytest.raises(ValueError, match="outside the air gap"):
        slotless_spm_airgap_field(0.028, [0.0], 2, 1.2, 1.05, 0.025, 0.03, 0.032, 0.8)
