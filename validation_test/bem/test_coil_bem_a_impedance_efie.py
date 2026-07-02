"""Golden regression test for the BEM-A impedance-EFIE R formulation.

Background (2026-07-02): the default BEM-A R is computed post-hoc from
the perfect-conductor (PEC) surface current, ``R = Re(Zs) * J^T M J``.
That OVER-estimates R for tightly-wound / near-contact / faceted coils
because the PEC current concentrates singularly at edges and
near-touching surfaces, where the Leontovich SIBC (which assumes J
smooth over the skin depth) breaks down.  On the kubota 3-turn coil the
PEC path gives 15.14 mOhm, a loss map shows 71% of that loss sits on 2%
of the surface area, and three independent methods (volume PEEC 3.7,
perimeter PEEC 4.5, analytic proximity 4.8) put the physical R at
~4.5-5 mOhm.

The fix (``--bema-impedance-efie``) puts Zs INTO the saddle system --
the (1,1) block becomes ``j w mu0 SL + Zs M`` (complex) -- so J is the
finite-impedance current and does NOT over-concentrate.  On the kubota
coil it gives 4.63 mOhm (matching the physics); on smooth geometry
(isolated straight wire) it reproduces the same Bessel R as the PEC
path (validated in C:/temp/volpeec_proto/imp_efie.py:
R_imp = R_pec = 0.3153 mOhm vs Bessel 0.3148).

This test locks the invariants on the self-contained gapped-torus
fixture (shared with test_coil_bem_a_volume_vol.py):
  * L is essentially unchanged (both PEC and impedance-EFIE L inside the
    PEEC band) -- the reactance is dominated by the geometry, not Zs;
  * R_impedance <= R_pec (the impedance-EFIE only ever REMOVES the
    spurious over-concentration, never adds);
  * both R are finite and positive.

See docs/peec/VOLUME_PEEC_DESIGN.md "2026-07-02 outcome".
"""
from __future__ import annotations

import os
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FIXTURE = os.path.join(HERE, "fixtures", "gapped_torus_2port.vol")
MAKE_FIXTURE = os.path.join(HERE, "fixtures", "make_gapped_torus_2port.py")

PEEC_BAND_NH = (80.0, 92.0)


def _have_ngsbem():
    try:
        import ngsolve.bem  # noqa: F401
        return True
    except Exception:
        return False


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


def _args(vol, impedance_efie):
    from argparse import Namespace
    return Namespace(
        coil_vol=vol,
        coil_source_name="source",
        coil_sink_name="sink",
        coil_rwg_quad_degree=5,
        coil_rwg_singular_nq=6,
        coil_bem_solver="auto",
        coil_aca_eps=1e-10,
        coil_gmres_tol=1e-10,
        coil_saddle_solver="auto",
        coil_sigma=5.8e7,
        frequency=150_000.0,
        current=1.0,
        n_threads=0,
        bema_impedance_efie=impedance_efie,
    )


@pytest.mark.slow
def test_impedance_efie_invariants(gapped_torus_vol):
    """PEC vs impedance-EFIE on the gapped-torus volume .vol."""
    if not _have_ngsbem():
        pytest.skip("ngsolve.bem not available")

    sys.path.insert(0, os.path.join(ROOT, "src", "radia", "panels"))
    from calc_inductance import _solve_coil_bem_a

    pec = _solve_coil_bem_a(_args(gapped_torus_vol, impedance_efie=False))
    imp = _solve_coil_bem_a(_args(gapped_torus_vol, impedance_efie=True))

    L_pec = float(pec["L_coil"]) * 1e9
    L_imp = float(imp["L_coil"]) * 1e9
    R_pec = float(pec["R_coil"]) * 1e3
    R_imp = float(imp["R_coil"]) * 1e3

    lo, hi = PEEC_BAND_NH
    assert lo <= L_pec <= hi, f"PEC L={L_pec:.3f} nH outside band [{lo},{hi}]"
    assert lo <= L_imp <= hi, f"IMP L={L_imp:.3f} nH outside band [{lo},{hi}]"

    # Both resistances finite + positive.
    assert R_pec > 0 and R_pec < 1e4, f"PEC R nonphysical: {R_pec}"
    assert R_imp > 0 and R_imp < 1e4, f"IMP R nonphysical: {R_imp}"

    # The impedance-EFIE only REMOVES the spurious PEC over-concentration:
    # R_imp <= R_pec (small numerical slack).  It must never exceed PEC.
    assert R_imp <= R_pec * 1.05, (
        f"impedance-EFIE R={R_imp:.4f} mOhm exceeds PEC R={R_pec:.4f} mOhm "
        f"-- the impedance-EFIE should only reduce over-concentration, "
        f"never increase R.  Regression in the (1,1)=jw*mu0*SL+Zs*M block "
        f"or the Z = Zs(J^H M J)+jw mu0(J^H SL J) extraction?")

    # residual clean for both dense LU solves.
    assert float(pec["bem_a_residual"]) < 1e-9
    assert float(imp["bem_a_residual"]) < 1e-9
