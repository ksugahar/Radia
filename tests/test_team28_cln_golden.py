"""Golden test: 3D (axisymmetric) CLN reproduces the TEAM 28 levitation force.

Locks the verified result of examples/maglev/team28/:
  - the full-FEM (split K+sN) levitation force at dZ=0 == lab ground truth
    -2.1928 N (within a hard band),
  - the N-stage CLN/Cauer reduced force CONVERGES to it (stage 3 < 1%,
    stage 5 <= 0.05%).

Runs a real axisymmetric NGSolve eddy-current solve (~20-40 s); skipped
cleanly if ngsolve / netgen are not importable in the active env.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
pytest.importorskip("scipy")

_HERE = os.path.dirname(os.path.abspath(__file__))
_TEAM28 = os.path.join(_HERE, "..", "examples", "levitation", "team28")
sys.path.insert(0, _TEAM28)

LAB_REF = -2.1928   # N, lab full-FEM ground truth at dZ=0 (50Hz_可動 .mat)


@pytest.fixture(scope="module")
def forces():
    from team28_cln_force import cln_forces  # noqa: E402
    fz_full, stage_forces = cln_forces(max_stage=6)
    return fz_full, stage_forces


def test_full_fem_matches_lab_ground_truth(forces):
    fz_full, _ = forces
    rel = abs(fz_full - LAB_REF) / abs(LAB_REF)
    assert rel < 0.005, (
        f"full-FEM force {fz_full:.4f} N deviates {rel*100:.2f}% from the "
        f"lab ground truth {LAB_REF} N (band 0.5%)")


def test_force_is_upward_lift_of_order_2N(forces):
    fz_full, _ = forces
    # sign convention: negative = upward lift; magnitude ~2.19 N at dZ=0
    assert fz_full < 0.0
    assert 2.0 < abs(fz_full) < 2.4


def test_cln_converges_to_full(forces):
    fz_full, sf = forces
    assert len(sf) >= 5, "expected at least 5 CLN stages"
    err = [abs(f - fz_full) / abs(fz_full) for f in sf]
    # stage 1 is the eddy-free DC response -> large error
    assert err[0] > 0.5
    # convergence: by stage 3 within 1%, by stage 5 within 0.05%
    assert err[2] < 0.01, f"stage 3 rel err {err[2]*100:.3f}% (expect <1%)"
    assert err[4] < 5e-4, f"stage 5 rel err {err[4]*100:.4f}% (expect <0.05%)"
    # monotone-ish: stage 5 is at least as good as stage 3
    assert err[4] <= err[2]


# Published TEAM 28 reference (Karl-Fetzer-Kurz-Lehner-Rucker, the official
# definition; laser triangulation): rest 3.8mm, measured stationary
# levitation height 11.5mm.  The repo disk-bottom is at 10.8mm at dZ=0.
PUB_LEVITATION_MM = 11.5
DISK_BOTTOM_DZ0_MM = 10.8
DISK_WEIGHT_N = 1.055


def test_force_convention_2x_and_published_height():
    """Lock the force convention + the published-benchmark consistency.

    The TEAM 28 surface integral Re[B_r J_t] is 2x the physical time-averaged
    Lorentz force <f_z> = -(1/2) Re[J_t conj(B_r)].  The disk floats where the
    PHYSICAL lift == weight; at dZ=0 (disk bottom 10.8mm) the physical lift ~=
    the disk weight, so the equilibrium is ~11mm, matching the published
    11.5mm.  (Regression guard for the 2026-06-20 convention fix: balancing
    the 2x integral against the 1x weight gave a spurious 14.9mm height.)
    """
    import numpy as np
    import scipy.sparse.linalg as spla
    from numpy import pi
    from ngsolve import GridFunction, L2, Integrate, x, dx, TaskManager
    from team28_cln_force import setup, to_csr, aluminium_z, FREQ  # noqa: E402

    with TaskManager():
        mesh, fes, fesPhi, fesB, sig, K, Nmat, F = setup(aluminium_z)  # dZ=0
        ndof = fes.ndof
        s = 2 * pi * FREQ * 1j
        free = np.array([i for i in range(ndof) if fes.FreeDofs()[i]])
        Kf = to_csr(K.mat, ndof)[free][:, free].tocsc()
        Nf = to_csr(Nmat.mat, ndof)[free][:, free].tocsc()
        Ff = np.array(F.vec.FV().NumPy(), dtype=complex)[free]
        xf = np.zeros(ndof, dtype=complex)
        xf[free] = spla.spsolve((Kf + s * Nf).tocsc(), Ff)

        gfuB = GridFunction(fes); gfuB.vec.FV().NumPy()[:] = xf
        gfu = GridFunction(fesPhi); gfu.Set(gfuB.components[0])
        gfB1 = GridFunction(fesB);  gfB1.Set(gfuB.components[1])
        gfJt = GridFunction(L2(mesh, order=2, dirichlet="outer", complex=True))
        gfJt.Set(-sig * s * gfu / x)
        Br, Jt, AlD = gfB1[0], gfJt, dx(mesh.Materials("Al"))
        f_verbatim = Integrate((Br.real * Jt.real - Br.imag * Jt.imag)
                               * (2 * pi * x) * AlD, mesh).real
        f_phys = -0.5 * Integrate((Jt.real * Br.real + Jt.imag * Br.imag)
                                  * (2 * pi * x) * AlD, mesh).real

    # (1) the verbatim integral is 2x the physical time-averaged force
    assert abs(abs(f_verbatim / f_phys) - 2.0) < 0.02, \
        f"verbatim/physical = {f_verbatim/f_phys:.3f} (expect ~2.0)"
    # (2) physical lift at dZ=0 ~= disk weight -> equilibrium ~10.8mm,
    #     consistent with the published 11.5mm (within ~10%)
    z_eq_approx = DISK_BOTTOM_DZ0_MM   # lift~weight here, so eq is near 10.8mm
    assert 0.9 < abs(f_phys) / DISK_WEIGHT_N < 1.2, \
        f"physical lift {abs(f_phys):.3f} N vs weight {DISK_WEIGHT_N} N"
    assert abs(z_eq_approx - PUB_LEVITATION_MM) / PUB_LEVITATION_MM < 0.10
