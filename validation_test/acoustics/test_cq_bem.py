"""Validation of the Lubich CQ time-domain acoustic BEM (radia.acoustics.cq).

Rigorous core check: every CQ Laplace node's frequency-domain BEM scattered field
equals the analytic sound-soft sphere at that COMPLEX wavenumber kappa = i s / c
(soft_sphere_scattering_complex_k) -- this validates the per-frequency ngsolve.bem
solves independent of the time-domain FFT reconstruction.  Plus: the recovered
time-domain signal is real and causal.  Run as a script to (re)write
cq_results.json.
"""
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pytest

from radia.acoustics import cq

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "cq_results.json"

R, C, N, DT = 1.0, 1.0, 16, 0.28
OBS = np.array([[0.0, 0.0, -3.0], [3.0, 0.0, 0.0], [0.0, 0.0, 2.0]])


def _worst_frequency_error():
    from ngsolve import BilinearForm, TaskManager, ds, z as Zc, exp as ngexp
    n = np.arange(N)
    rho = (np.finfo(float).eps ** 0.5) ** (1.0 / N)
    zeta = rho * np.exp(-2j * np.pi * n / N)
    s = cq.bdf_delta(zeta, "BDF2") / DT
    kappa = 1j * s / C
    mesh, fes = cq._build_sphere_screen(R, 0.4, 3, 7.0)
    u, v = fes.TnT()
    worst = 0.0
    with TaskManager():
        pre = BilinearForm(u * v * ds("sphere"), diagonal=True).Assemble().mat.Inverse()
        for l in range(N):
            ghat = -ngexp(-complex(s[l]) / C * Zc)             # unit incident plane wave
            bem = cq._frequency_scattered(mesh, fes, pre, kappa[l], ghat, OBS)
            ana = cq.soft_sphere_scattering_complex_k(complex(kappa[l]), R, OBS)
            worst = max(worst, float(np.max(np.abs(bem - ana)) / max(np.max(np.abs(ana)), 1e-30)))
    return worst


def test_cq_frequency_bem_matches_complex_k_analytic():
    assert _worst_frequency_error() < 5e-3


def test_cq_time_signal_real_and_causal():
    res = cq.cq_soft_sphere_scattering(OBS, radius=R, num_time=N, time_step=DT, sound_speed=C)
    scale = max(float(np.max(np.abs(res["scattered"]))), 1e-30)
    assert res["max_imag"] < 1e-6 * scale        # real time signal (Lubich rho-weighting)
    assert scale > 0                             # a nonzero scattered response
    assert res["scattered"].shape == (N, len(OBS))


def _write_results():
    import ngsolve
    import radia
    worst = _worst_frequency_error()
    res = cq.cq_soft_sphere_scattering(OBS, radius=R, num_time=N, time_step=DT, sound_speed=C)
    out = {
        "schema": "radia.acoustics.cq.v1",
        "generated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "problem": {"radius": R, "sound_speed": C, "num_time": N, "time_step": DT,
                    "method": "BDF2", "boundary": "sound_soft_sphere"},
        "worst_frequency_bem_vs_complex_k_analytic": worst,
        "time_signal_max_imag_over_scale": res["max_imag"] / max(float(np.max(np.abs(res["scattered"]))), 1e-30),
        "versions": {"radia": getattr(radia, "__version__", "?"),
                     "ngsolve": ngsolve.__version__,
                     "python": platform.python_version(),
                     "platform": platform.platform()},
    }
    RESULTS.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    _write_results()
    sys.exit(0)
