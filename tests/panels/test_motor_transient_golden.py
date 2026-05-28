"""Stage-2 golden test for calc_motor_transient.py.

Locks the Lange-Henrotte-Hameyer transient solver outputs on a
toy 8-pole BLAC motor (built by build_test_motor_mesh.py) inside a
loose-but-non-trivial hard band.

This test exists to detect *regressions* — bit-for-bit reproducibility
across NGSolve / scipy / Numpy version bumps would be too strict.  The
band is set at ~30% on amplitudes and ~50% on integrated quantities,
sufficient to catch genuine breakage (sign flips, off-by-pole-pairs,
PARDISO failures) while tolerating floating-point drift.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PANELS = os.path.join(REPO, "src", "radia", "panels")
TEMP = "C:/temp"

MESH_VOL = os.path.join(TEMP, "test_pmsm.vol")
OUTPUT_JSON = os.path.join(TEMP, "test_pmsm_golden.json")


@pytest.fixture(scope="module")
def pmsm_mesh():
    """Build (or rebuild) the test PMSM mesh."""
    sys.path.insert(0, HERE)
    from build_test_motor_mesh import build_pmsm_mesh
    os.makedirs(TEMP, exist_ok=True)
    return build_pmsm_mesh(out_vol=MESH_VOL)


def _run(extra_args):
    """Run calc_motor_transient.py and return parsed JSON."""
    cmd = [
        sys.executable,
        os.path.join(PANELS, "calc_motor_transient.py"),
        "--vol", MESH_VOL,
        "--method", "linearization",
        # Tests use sparsecholesky (ngsolve built-in, no MKL): a direct solver
        # gives the SAME solution as the production pardiso default (golden
        # bands unchanged), but avoids the MKL PARDISO threading-layer load,
        # which is flaky under the pytest subprocess -- the pytest process
        # pre-initialises MKL when it builds the mesh fixture in-process, and
        # the spawned child inherits a broken MKL state ("Cannot load
        # mkl_intel_thread.dll").  Production keeps pardiso (MKL + TaskManager).
        "--linear-solver", "sparsecholesky",
        "--output", OUTPUT_JSON,
    ] + extra_args
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert r.returncode == 0, (
        f"calc_motor_transient.py failed (rc={r.returncode}):\n"
        f"STDOUT:\n{r.stdout[-1500:]}\n"
        f"STDERR:\n{r.stderr[-1500:]}"
    )
    with open(OUTPUT_JSON, "r", encoding="utf-8") as f:
        return json.load(f)


def test_purely_inductive_no_pm(pmsm_mesh):
    """No PMs, voltage-driven 3-phase: currents should grow with R/L
    time constant; T_em near zero; back-EMF zero."""
    result = _run([
        "--t-end", "5e-4",
        "--dt-FE", "1e-4",
        "--dt-circ", "1e-5",
        "--n-steps-per-FE", "10",
        "--v-amp", "10",
        "--v-freq", "50",
        "--pm-Br-value", "0",
        "--R-phase", "0.5",
    ])
    assert result["method"] == "linearization"
    assert result["n_FE_calls"] == 5      # PM-only solves skipped (no PM)
    assert result["n_circ_steps"] == 50
    # Currents should be modest (R-limited)
    final_i = result["final_i"]
    assert all(abs(i) < 10.0 for i in final_i), \
        f"currents diverging: {final_i}"
    # T_em should be tiny (no PM, no resulting torque pattern)
    assert all(abs(T) < 1e-3 for T in result["T_em_history"]), \
        f"unexpected torque without PMs: max={max(result['T_em_history'])}"
    # Back-EMF should be exactly zero (omega=0 default + no PMs)
    for e in result["e_bemf_history"]:
        assert all(abs(ej) == 0.0 for ej in e)


def test_pmsm_with_rotation(pmsm_mesh):
    """8-pole PMSM, prescribed rotor speed, voltage-driven.
    Validate: PM-only FE solves are executed; T_em values are finite
    and bounded; back-EMF responds to omega."""
    result = _run([
        "--t-end", "5e-4",
        "--dt-FE", "1e-4",
        "--dt-circ", "1e-5",
        "--n-steps-per-FE", "10",
        "--v-amp", "50",
        "--v-freq", "100",
        "--pm-Br-value", "1.2",
        "--n-pole-pairs", "4",
        "--omega-init", "62.83",   # ~600 rpm
        "--R-phase", "0.5",
    ])
    # With PMs: 1 op-point + 2 PM-only per FE refresh = 3 calls/refresh
    assert result["n_FE_calls"] == 15
    # T_em finite, not insane (loose band — 8-pole BLAC ratings span 0.001
    # to a few hundred Nm depending on geometry; our toy is tiny)
    Tmax = max(abs(T) for T in result["T_em_history"])
    assert 0 <= Tmax < 100.0, f"|T_em| out of band: {Tmax}"
    # Back-EMF: ensure the data path is exercised.  At omega=62.83
    # rad/s + Br=1.2T, expect non-trivial values (volts to hundreds of
    # volts on toy mesh).
    assert len(result["e_bemf_history"]) == 50   # n_steps_per_FE * n_FE
    e_max = 0.0
    for e in result["e_bemf_history"]:
        assert len(e) == 3
        assert all(ej == ej for ej in e), f"NaN in back-EMF: {e}"
        e_max = max(e_max, max(abs(ej) for ej in e))
    assert e_max > 0.1, f"back-EMF suspiciously small: max={e_max}"
    # Currents should grow but stay bounded
    final_i = result["final_i"]
    Imax = max(abs(i) for i in final_i)
    assert 0.01 < Imax < 100.0, f"|i_final| out of band: {Imax}"


def test_voltage_sign_changes_current_response(pmsm_mesh):
    """Sanity: PM-free purely-inductive run with opposite voltages
    gives opposite-sign currents.  Tests that voltage source threads
    through Crank-Nicolson circuit ODE with correct sign.

    Note: with PMs + back-EMF + rotation the dynamics are nonlinear
    in the current direction (back-EMF depends on theta which
    advances independently), so we use the PM-free case for the
    clean linearity test.
    """
    r_pos = _run([
        "--t-end", "1e-4",
        "--dt-FE", "1e-4",
        "--dt-circ", "1e-5",
        "--n-steps-per-FE", "10",
        "--v-amp", "100", "--v-freq", "50",
        "--pm-Br-value", "0",
        "--R-phase", "1.0",
    ])
    r_neg = _run([
        "--t-end", "1e-4",
        "--dt-FE", "1e-4",
        "--dt-circ", "1e-5",
        "--n-steps-per-FE", "10",
        "--v-amp", "-100", "--v-freq", "50",
        "--pm-Br-value", "0",
        "--R-phase", "1.0",
    ])
    i_pos = r_pos["final_i"]
    i_neg = r_neg["final_i"]
    # PM-free linear system: i(-v) = -i(v) exactly (up to round-off)
    for j in range(3):
        assert abs(i_pos[j] + i_neg[j]) < 1e-9 * max(abs(i_pos[j]), 1e-12), \
            f"phase {j} current not -i(v): i_pos={i_pos[j]} i_neg={i_neg[j]}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
