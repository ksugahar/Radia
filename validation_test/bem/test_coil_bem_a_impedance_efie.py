"""Golden regression tests for the unified BEM-A impedance-EFIE.

Since 2026-07-02 the impedance-EFIE is the SOLE BEM-A formulation: the
Leontovich surface impedance sits inside the saddle system,

    [ jw*mu0*SL + Zs*M   D^T ] [J]   [0]
    [ D                   0  ] [p] = [g]     (complex, omega > 0)

so J is the finite-impedance surface current; R = Re(Zs)*(J^H M J),
L = mu0*(J^H SL J) (external inductance).  The historical PEC post-hoc
path (real perfect-conductor saddle + R = Re(Zs)*J^T M J evaluated
afterwards) was REMOVED: the PEC J concentrates singularly at
near-contact gaps / edges where the Leontovich integral breaks down,
over-estimating R ~3x on tightly-wound coils (kubota 3-turn pancake:
PEC 15.14 mOhm vs impedance-EFIE 4.63 mOhm, the latter confirmed by
volume PEEC 3.7 / perimeter PEEC 4.5 / analytic proximity 4.8).  On
smooth geometry both agreed (isolated wire = Bessel to <1%), so the
removal loses nothing.  See docs/peec/VOLUME_PEEC_DESIGN.md
"2026-07-02 outcome".

These tests lock the unified path:

1. ``test_impedance_efie_torus_absolute``: absolute R/L bands on the
   self-contained gapped-torus fixture at 150 kHz and 7 kHz, plus the
   mesh-robust physics lock R(150k)/R(7k) ~ sqrt(150/7) (strong-skin
   Leontovich R ~ sqrt(f)), plus the omega=0 DC branch (R == 0, same
   external L).  Reference values captured 2026-07-02 on LAB:
   150 kHz R=1.0453 mOhm L=87.270 nH; 7 kHz R=0.2231 mOhm L=87.332 nH.
2. ``test_impedance_efie_wire_bessel``: BEM-A on an isolated straight
   round wire reproduces the closed-form Bessel AC resistance -- the
   physics-anchored calibration that needs no captured numbers.
3. ``test_minres_rejected_for_complex``: the AC saddle is complex and
   scipy MINRES silently solves only its real part -- the library must
   fail fast instead.
"""
from __future__ import annotations

import math
import os
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURE = os.path.join(HERE, "fixtures", "gapped_torus_2port.vol")
MAKE_FIXTURE = os.path.join(HERE, "fixtures", "make_gapped_torus_2port.py")

SIGMA_CU = 5.8e7
MU_0 = 4e-7 * math.pi

# PEEC cross-method band for the gapped-torus L (see
# test_coil_bem_a_volume_vol.py; the .vol fixture is regenerated
# per-machine so bands are deliberately generous).
PEEC_BAND_NH = (80.0, 92.0)
# Impedance-EFIE R bands captured 2026-07-02 (LAB), +-15% for
# cross-machine mesh variation of the auto-generated fixture.
R_BAND_150K_MOHM = (0.89, 1.20)    # captured 1.0453
R_BAND_7K_MOHM = (0.19, 0.26)      # captured 0.2231


def _have_ngsbem():
    try:
        import ngsolve.bem  # noqa: F401
        return True
    except Exception:
        return False


def _zs(sigma, omega):
    delta = math.sqrt(2.0 / (omega * MU_0 * sigma))
    return (1.0 + 1.0j) / (sigma * delta)


@pytest.fixture(scope="module")
def gapped_torus_vol():
    if not os.path.isfile(FIXTURE):
        if not os.path.isfile(MAKE_FIXTURE):
            pytest.skip(f"fixture maker missing: {MAKE_FIXTURE}")
        import subprocess
        rc = subprocess.run([sys.executable, MAKE_FIXTURE],
                            capture_output=True, text=True)
        if rc.returncode != 0 or not os.path.isfile(FIXTURE):
            pytest.skip(
                f"fixture make failed (rc={rc.returncode}): "
                f"{rc.stderr.splitlines()[-3:] if rc.stderr else ''}")
    return FIXTURE


def _surface_mesh(volpath):
    from ngsolve import Mesh
    sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
    from surface_mesh_extract import _extract_surface_mesh_filtered
    mesh = Mesh(volpath)
    if mesh.ne > 0:
        mesh = _extract_surface_mesh_filtered(mesh, keep_label="")
    return mesh


@pytest.mark.slow
def test_impedance_efie_torus_absolute(gapped_torus_vol):
    """Absolute R/L bands + sqrt(f) physics lock + DC branch."""
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    from ngsolve import TaskManager
    from radia.bem.coil_inductance_ngsolve import (
        compute_inductance_source_sink)

    results = {}
    with TaskManager():
        mesh = _surface_mesh(gapped_torus_vol)
        for freq in (150e3, 7e3):
            omega = 2 * math.pi * freq
            results[freq] = compute_inductance_source_sink(
                mesh, "source", "sink",
                omega=omega, Z_s_complex=_zs(SIGMA_CU, omega))
        dc = compute_inductance_source_sink(mesh, "source", "sink",
                                            omega=0.0)

    r150 = results[150e3]
    r7 = results[7e3]

    for tag, r in (("150k", r150), ("7k", r7), ("DC", dc)):
        L_nH = float(r["L"]) * 1e9
        lo, hi = PEEC_BAND_NH
        assert lo <= L_nH <= hi, (
            f"{tag}: external L={L_nH:.3f} nH outside PEEC band "
            f"[{lo},{hi}] -- external-inductance convention regressed?")
        assert float(r["residual"]) < 1e-9

    R150 = float(r150["R"]) * 1e3
    R7 = float(r7["R"]) * 1e3
    assert R_BAND_150K_MOHM[0] <= R150 <= R_BAND_150K_MOHM[1], (
        f"impedance-EFIE R(150kHz)={R150:.4f} mOhm outside "
        f"{R_BAND_150K_MOHM} (captured 1.0453)")
    assert R_BAND_7K_MOHM[0] <= R7 <= R_BAND_7K_MOHM[1], (
        f"impedance-EFIE R(7kHz)={R7:.4f} mOhm outside "
        f"{R_BAND_7K_MOHM} (captured 0.2231)")

    # Strong-skin Leontovich scaling R ~ sqrt(f): mesh-independent
    # physics lock (captured ratio 4.686 vs sqrt(150/7)=4.629, +1.2%).
    ratio = R150 / R7
    expect = math.sqrt(150.0 / 7.0)
    assert abs(ratio / expect - 1.0) < 0.10, (
        f"R(150k)/R(7k)={ratio:.3f} deviates >10% from sqrt(f) scaling "
        f"{expect:.3f} -- Zs frequency dependence regressed?")

    # DC branch: R == 0, external L consistent with AC (same geometry
    # quantity; captured spread 87.27 vs 87.33 nH = 0.07%).  The 7 kHz
    # comparison is the SENSITIVE convention lock: if the internal
    # surface reactance Im(Zs)*(J^H M J)/omega ever gets folded back
    # into L, L(7 kHz) shifts ~+5.9% (vs only ~+1.2% at 150 kHz where
    # delta is small), so DC-vs-7k at 2% catches the regression that
    # DC-vs-150k barely sees.
    assert float(dc["R"]) == 0.0
    assert abs(float(dc["L"]) / float(r150["L"]) - 1.0) < 0.02
    assert abs(float(dc["L"]) / float(r7["L"]) - 1.0) < 0.02, (
        f"DC L={float(dc['L'])*1e9:.3f} nH vs 7kHz L="
        f"{float(r7['L'])*1e9:.3f} nH differ >2% -- was the internal "
        f"surface reactance folded into L?  L must stay the EXTERNAL "
        f"mu0*(J^H SL J) convention.")
    # Im(J) GridFunction exists and is all-zero at DC.
    import numpy as np
    assert float(np.max(np.abs(
        dc["gf_J_im"].vec.FV().NumPy()))) == 0.0


@pytest.mark.slow
def test_impedance_efie_wire_bessel():
    """Isolated straight wire: BEM-A impedance-EFIE == Bessel R_ac.

    Physics-anchored calibration (no captured mesh-dependent numbers):
    a=3.15 mm Cu wire, L=60 mm, 150 kHz (a/delta=18.5).  The Leontovich
    R of the impedance-EFIE current must match the closed-form Bessel
    cylinder impedance within 3% (measured 0.5% at maxh=1.5 mm).
    """
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")
    try:
        from netgen.occ import Cylinder, OCCGeometry, Pnt, Z
    except Exception:
        pytest.skip("netgen.occ not available")
    from ngsolve import Mesh, TaskManager
    from radia.bem.coil_inductance_ngsolve import (
        compute_inductance_source_sink)
    from radia.analytical_formulas.conductor_impedance import (
        cylinder_ac_impedance)

    a = 3.15e-3
    Lw = 60e-3
    freq = 150e3
    omega = 2 * math.pi * freq

    cyl = Cylinder(Pnt(0, 0, 0), Z, r=a, h=Lw)
    for f in cyl.faces:
        z = f.center[2]
        f.name = ("source" if abs(z) < Lw * 0.25
                  else "sink" if abs(z - Lw) < Lw * 0.25 else "wire")
    cyl.name = "cond"
    sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
    from surface_mesh_extract import _extract_surface_mesh_filtered
    with TaskManager():
        mesh = Mesh(OCCGeometry(cyl).GenerateMesh(maxh=1.5e-3))
        # compute_inductance_source_sink assumes a PURE SURFACE mesh
        # (volume tets add saddle null modes the D[:-1,:] deflation
        # cannot remove -> singular LU; same pathology the panel fixed
        # in _build_bema_coil_mesh, see test_coil_bem_a_volume_vol.py).
        mesh = _extract_surface_mesh_filtered(mesh, keep_label="")
        res = compute_inductance_source_sink(
            mesh, "source", "sink",
            omega=omega, Z_s_complex=_zs(SIGMA_CU, omega))

    R_bessel = cylinder_ac_impedance(a, SIGMA_CU, omega).real * Lw
    R = float(res["R"])
    assert abs(R / R_bessel - 1.0) < 0.03, (
        f"impedance-EFIE wire R={R*1e3:.4f} mOhm vs Bessel "
        f"{R_bessel*1e3:.4f} mOhm ({(R/R_bessel-1)*100:+.2f}%) -- "
        f"the unified BEM-A no longer reproduces the closed-form "
        f"self-skin resistance.")


def test_minres_rejected_for_complex():
    """scipy MINRES assumes a Hermitian operator; the AC saddle is
    complex-SYMMETRIC, so minres would silently return a wrong
    solution -- the saddle helper must fail fast instead."""
    import numpy as np
    from radia.bem.coil_inductance_ngsolve import _solve_saddle

    K = np.eye(4, dtype=complex) + 1j * np.eye(4)
    rhs = np.ones(4, dtype=complex)
    with pytest.raises(ValueError, match="MINRES"):
        _solve_saddle(K, rhs, method="minres")
    # Real systems keep MINRES (omega=0 vacuum-L path).
    x, info = _solve_saddle(np.eye(4), np.ones(4), method="minres")
    assert np.allclose(x, 1.0)


def test_gmres_complex_matches_lu():
    """GMRES is the panel's auto choice for large AC saddles: lock that
    the gmres path handles a COMPLEX-SYMMETRIC system correctly (same
    answer as dense LU) on a small deterministic instance."""
    import numpy as np
    from radia.bem.coil_inductance_ngsolve import _solve_saddle

    rng = np.random.default_rng(20260702)
    n = 24
    A = rng.standard_normal((n, n)) + 1j * rng.standard_normal((n, n))
    K = A + A.T + 4.0 * n * np.eye(n)      # complex SYMMETRIC (K^T=K)
    assert np.allclose(K, K.T)
    rhs = rng.standard_normal(n) + 1j * rng.standard_normal(n)

    x_lu, _ = _solve_saddle(K.copy(), rhs, method="lu")
    x_gm, _ = _solve_saddle(K.copy(), rhs, method="gmres", tol=1e-12)
    assert np.allclose(x_gm, x_lu, rtol=1e-6, atol=1e-9), (
        "gmres on a complex-symmetric system diverged from the LU "
        "reference -- complex dtype handling regressed?")


def test_parameter_validation_fails_before_assembly():
    """The omega/Z_s_complex contract is validated BEFORE the expensive
    dense assembly: a bad call must raise instantly, even with
    mesh=None (the mesh is not touched until after validation)."""
    from radia.bem.coil_inductance_ngsolve import (
        compute_inductance_source_sink)

    # AC without the surface impedance -> instant ValueError.
    with pytest.raises(ValueError, match="Z_s_complex"):
        compute_inductance_source_sink(None, omega=2 * math.pi * 1e3)
    # Negative omega must not silently select the DC vacuum branch.
    with pytest.raises(ValueError, match="omega must be > 0"):
        compute_inductance_source_sink(None, omega=-1.0,
                                       Z_s_complex=1 + 1j)
    # NaN omega likewise.
    with pytest.raises(ValueError, match="omega must be > 0"):
        compute_inductance_source_sink(None, omega=float("nan"),
                                       Z_s_complex=1 + 1j)
