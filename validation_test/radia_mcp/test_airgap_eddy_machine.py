r"""AGE rotating-machine core, COMPLETE: two FE regions + EDDY currents + rotation + torque.

The capstone of the AGE line: steps (2) two-region coupling, (3) rotation and (4) torque were
each shown separately (and the complex gap DtN in test_airgap_eddy_coupling.py was a single
region).  Here the reusable solver (radia_ngsolve.airgap_machine) does ALL of it at once -- a
rotating machine whose STATOR is a CONDUCTOR (jw*mu*sigma eddy currents), coupled to the rotor
across the UN-MESHED gap, rotated by a pure phase, with the mesh-free closed-form torque.

Validation (independent, fully-meshed complex reference):
  AGE       : mesh ONLY rotor (R0<r<RI, Laplace) + stator (RO<r<REXT, CONDUCTING); gap un-meshed.
  reference : mesh the WHOLE R0<r<REXT (gap meshed, sigma=0; stator conducts), Dirichlet rings.
  (1) field : the AGE complex field over the conducting stator equals the fully-meshed field.
  (2) torque: the AGE closed-form rotation torque T(theta_r) equals the meshed ring-phasor torque
              and genuinely swings with the rotor angle.
The stator field is strongly complex (|Im/Re|>>1 -> genuine eddy phase), so this is the
eddy x mesh-free-gap rotating-machine capability, not a magnetostatic one.
"""
import math
import os
import sys

import numpy as np
from ngsolve import (Mesh, H1, GridFunction, BilinearForm, grad, dx, x, y, atan2, cos, sin,
                     ds, BND, Integrate, TaskManager)
from netgen.geom2d import SplineGeometry

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.airgap_machine import (airgap_coupling, airgap_factorize,
                                                    airgap_solve, airgap_torque, _inv)
from radia_mcp.radia_ngsolve.airgap_element import airgap_harmonic_torque

R0, RI, RO, REXT = 0.05, 0.10, 0.11, 0.20
N = 2
A_ROT, A_STA, LSTK = 1.0, 0.3, 0.05
MU0 = 4e-7 * math.pi
SIGMA, FREQ = 2.0e6, 200.0                # ~2.5 cm skin depth in the 9 cm stator -> strong eddy
OMEGA = 2.0 * math.pi * FREQ


def _machine_mesh(meshed_gap):
    geo = SplineGeometry()
    geo.AddCircle((0, 0), REXT, leftdomain=2, rightdomain=0, bc="outer")
    if meshed_gap:
        geo.AddCircle((0, 0), RO, leftdomain=3, rightdomain=2, bc="stator_ring")
        geo.AddCircle((0, 0), RI, leftdomain=1, rightdomain=3, bc="rotor_ring")
        geo.AddCircle((0, 0), R0, leftdomain=0, rightdomain=1, bc="rotor_inner")
        geo.SetMaterial(1, "rotor"); geo.SetMaterial(2, "stator"); geo.SetMaterial(3, "gap")
    else:
        geo.AddCircle((0, 0), RO, leftdomain=0, rightdomain=2, bc="stator_ring")
        geo.AddCircle((0, 0), RI, leftdomain=1, rightdomain=0, bc="rotor_ring")
        geo.AddCircle((0, 0), R0, leftdomain=0, rightdomain=1, bc="rotor_inner")
        geo.SetMaterial(1, "rotor"); geo.SetMaterial(2, "stator")
    m = Mesh(geo.GenerateMesh(maxh=0.006)); m.Curve(6)
    return m


def _excite(mesh, thr):
    return mesh.BoundaryCF({"rotor_inner": A_ROT * cos(N * atan2(y, x) - N * thr),
                            "outer": A_STA * cos(N * atan2(y, x))})


def _stiffness(mesh):
    sigma = mesh.MaterialCF({"stator": SIGMA}, default=0.0)
    fes = H1(mesh, order=6, complex=True, dirichlet="outer|rotor_inner")
    u, v = fes.TnT()
    a = BilinearForm(fes)
    a += grad(u) * grad(v) * dx + 1j * OMEGA * MU0 * sigma * u * v * dx
    with TaskManager():
        a.Assemble()
    return fes, a


def _age_solve(thetas):
    mesh = _machine_mesh(meshed_gap=False)
    fes, a = _stiffness(mesh)
    coupling = airgap_coupling(fes, RI, RO, "rotor_ring", "stator_ring", [N])
    factor = airgap_factorize(a.mat, coupling, fes.FreeDofs())
    with TaskManager():
        gfus = [airgap_solve(factor, fes, dirichlet_cf=_excite(mesh, thr)) for thr in thetas]
    return mesh, coupling, gfus


def _ref_solve(thetas):
    mesh = _machine_mesh(meshed_gap=True)
    fes, a = _stiffness(mesh)
    gfus = []
    with TaskManager():
        Kinv = _inv(a.mat, fes.FreeDofs())
        for thr in thetas:
            gfu = GridFunction(fes)
            gfu.Set(_excite(mesh, thr), BND)
            rhs = gfu.vec.CreateVector(); rhs.data = -(a.mat * gfu.vec)
            gfu.vec.data += Kinv * rhs
            gfus.append(gfu)
    return mesh, gfus


def _ring_phasor(mesh, gfu, bc, n, R):
    cn, sn = cos(n * atan2(y, x)), sin(n * atan2(y, x))
    norm = math.pi * R
    ac = Integrate(cn * gfu * ds(definedon=mesh.Boundaries(bc)), mesh)
    as_ = Integrate(sn * gfu * ds(definedon=mesh.Boundaries(bc)), mesh)
    return ac / norm - 1j * as_ / norm


def test_eddy_machine_age_matches_fully_meshed():
    thetas = [k * (2 * math.pi / N) / 8 for k in range(8)]
    mesh_a, coupling, gfus_a = _age_solve(thetas)
    mesh_r, gfus_r = _ref_solve(thetas)

    # (1) field over the conducting stator at one angle
    k0 = 2
    diffs, refs = [], []
    for frac in (0.2, 0.5, 0.8):
        rr = RO + frac * (REXT - RO)
        for th in (0.3, 1.2, 2.5, -0.8):
            px, py = rr * math.cos(th), rr * math.sin(th)
            diffs.append(abs(gfus_a[k0](mesh_a(px, py)) - gfus_r[k0](mesh_r(px, py))))
            refs.append(abs(gfus_r[k0](mesh_r(px, py))))
    rel_field = max(diffs) / max(refs)
    pmid = gfus_r[k0](mesh_r((RO + REXT) / 2.0, 0.0))

    # (2) torque vs rotor angle: AGE closed-form vs meshed ring-phasor (same closed form)
    Tage = np.array([airgap_torque(coupling, g, axial_length=LSTK) for g in gfus_a])
    Tref = np.array([airgap_harmonic_torque(N, RI, RO,
                                            _ring_phasor(mesh_r, g, "rotor_ring", N, RI),
                                            _ring_phasor(mesh_r, g, "stator_ring", N, RO),
                                            axial_length=LSTK) for g in gfus_r])
    rel_T = np.max(np.abs(Tage - Tref)) / np.max(np.abs(Tref))
    swing = (Tage.real.max() - Tage.real.min()) / np.max(np.abs(Tage.real))
    print(f"eddy machine AGE vs fully-meshed: field rel={rel_field:.2e}; torque rel={rel_T:.2e}; "
          f"swing={swing:.2f}; stator A={pmid:.4f} (|Im/Re|={abs(pmid.imag/pmid.real):.1f} eddy)")

    assert abs(pmid.imag) > 0.5 * abs(pmid.real), "stator should be strongly complex (genuine eddy)"
    assert rel_field < 1e-4, f"AGE eddy field off fully-meshed by {rel_field:.2e}"
    assert rel_T < 1e-2, f"AGE rotation torque off meshed reference by {rel_T:.2e}"
    assert swing > 0.3, "torque should swing with the rotor angle"


def main():
    test_eddy_machine_age_matches_fully_meshed()
    print("[OK] AGE rotating-machine core COMPLETE: rotor + EDDY-current stator coupled across the "
          "UN-MESHED gap, rotated by a pure phase, with mesh-free closed-form torque -- reproduces "
          "the fully-meshed complex reference. Eddy x mesh-free gap + rotation + torque.")


if __name__ == "__main__":
    main()
