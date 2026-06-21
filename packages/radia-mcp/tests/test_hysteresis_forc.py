# -*- coding: utf-8 -*-
r"""FORC (first-order reversal curves) for the B-input Play model -- physics-gated.

Learns the FORC / Preisach-identification technique that HysterSoft (Dimian-Andrei, FAMU-FSU / TU
Vienna / Cuza University -- "Scalar and vector hysteresis simulations using HysterSoft") performs on
measured loops, and brings it to the radia-ngsolve energy-based Play model.  Because the Play model is
mathematically the static Preisach model (Bobbio 1997), its FORC distribution is exactly computable:
rho(Ba,Bb) = -1/2 d2H/dBa dBb is a set of RIDGES at the coercivity (Bb-Ba)/2 = eta_k with integrated
weight a_k(Bsat-eta_k) -- i.e. the play thresholds eta_k ARE the B-space Preisach / FORC density and
the slopes a_k set the ridge weights.  Gated here against that closed form.
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.hysteresis import (PlayHysteresis, forc_distribution,
                                                forc_coercivity_weight, identify_from_forc)


def _model():
    return PlayHysteresis(eta=[0.0, 0.3, 0.6, 0.9], a=[200.0, 150.0, 100.0, 80.0])


def test_forc_family_structure():
    """The FORC family H(Ba, Bb) is defined only on the upper triangle Bb >= Ba (the up-sweep)."""
    m = _model()
    grid, H = m.forc_curves(1.3, n=61)
    n = len(grid)
    assert np.isnan(H[10, 3]) and np.isfinite(H[3, 10]), "H defined only for Bb >= Ba"
    assert np.all(np.isfinite(np.diag(H))), "the reversal points (Bb=Ba) are measured"


def test_forc_ridges_at_play_thresholds():
    """rho ridges sit at the coercivity (Bb-Ba)/2 = eta_k, with integrated weight a_k (Bsat - eta_k)."""
    m = _model(); Bsat = 1.3
    grid, H = m.forc_curves(Bsat, n=121)
    rho = forc_distribution(grid, H)
    for eta_k, w in m.analytic_forc_weights(Bsat):                 # closed-form ridge weights
        num = forc_coercivity_weight(grid, rho, eta_k, 0.06)
        print(f"ridge eta={eta_k:.2f}: weight {num:.1f} vs analytic a*(Bsat-eta)={w:.1f}")
        assert abs(num - w) / w < 0.1, f"FORC ridge at eta={eta_k} off by {abs(num-w)/w:.2f}"


def test_forc_distribution_concentrates_on_ridges():
    """Between thresholds the FORC distribution is ~0 (it is concentrated on the eta_k ridges)."""
    m = _model(); Bsat = 1.3
    grid, H = m.forc_curves(Bsat, n=121)
    rho = forc_distribution(grid, H)
    on = forc_coercivity_weight(grid, rho, 0.30, 0.04)            # on the eta_1 ridge
    off = forc_coercivity_weight(grid, rho, 0.45, 0.04)          # between eta_1 and eta_2
    print(f"on-ridge(eta=0.3)={on:.1f}, off-ridge(Bc=0.45)={off:.2f}")
    assert abs(off) < 0.05 * abs(on), "the FORC distribution must vanish away from the eta_k ridges"


def test_forc_total_weight():
    """INT INT rho dBa dBb over the whole plane = sum_k a_k (Bsat - eta_k)."""
    m = _model(); Bsat = 1.3
    grid, H = m.forc_curves(Bsat, n=121)
    rho = forc_distribution(grid, H)
    total = float(np.nansum(rho) * (grid[1] - grid[0]) ** 2)
    wsum = sum(w for _, w in m.analytic_forc_weights(Bsat))
    print(f"total FORC weight {total:.1f} vs sum a_k(Bsat-eta_k) {wsum:.1f}")
    assert abs(total - wsum) / wsum < 0.05


def test_identify_from_forc_roundtrip():
    """INVERSE FORC identification (HysterSoft workflow): recover the cell slopes a_k from a FORC family
    by reading the ridge comb, and reproduce an unseen loop."""
    eta = [0.0, 0.3, 0.6, 0.9]; a_true = [200.0, 150.0, 100.0, 80.0]
    m0 = PlayHysteresis(eta=eta, a=a_true)
    Bsat = 1.3
    grid, H = m0.forc_curves(Bsat, n=161)
    m1 = identify_from_forc(grid, H, eta, Bsat, a0=200.0)
    rel = max(abs(m1.a[k] - a_true[k]) / a_true[k] for k in range(1, len(a_true)))
    Bt = 1.1
    gen = abs(m1.loss_per_cycle(Bt) - m0.loss_per_cycle(Bt)) / m0.loss_per_cycle(Bt)
    print(f"identify from FORC: recovered a={[round(x,1) for x in m1.a]} (true {a_true}), "
          f"max rel={rel:.3f}, unseen-loop rel={gen:.3f}")
    assert rel < 0.06, f"recovered slopes off by {rel:.2f}"
    assert gen < 0.05, "the FORC-identified model must reproduce an unseen loop"
    assert m1.a[0] == 200.0, "the reversible slope a0 is the separate loop-tip input (no FORC ridge)"


def main():
    test_forc_family_structure()
    test_forc_ridges_at_play_thresholds()
    test_forc_distribution_concentrates_on_ridges()
    test_forc_total_weight()
    test_identify_from_forc_roundtrip()
    print("[OK] FORC for the B-input Play model: first-order reversal curves -> FORC distribution "
          "rho=-1/2 d2H/dBa dBb ridges at the play thresholds eta_k with weight a_k(Bsat-eta_k) "
          "(the B-space Preisach density; the HysterSoft FORC technique on the Play == Preisach model).")


if __name__ == "__main__":
    main()
