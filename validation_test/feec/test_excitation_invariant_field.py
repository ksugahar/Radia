"""Golden: EXCITATION-invariant flux lines (same field-line shape as the drive rises).

Locks the finding of the docs notebook + helper
docs/clebsch_hodograph/excitation_invariant_field.{ipynb,py} on a small-gap (24 mm),
near-knee bending end pack:

  (1) LINEARITY => invariance.  With mu forced constant the flux-line DIRECTION drift
      D_dir(B_drive) is 0 at EVERY drive -- scaling the current scales B everywhere, so
      the field-line pattern is identical.  Saturation (nonlinear mu) is the SOLE thing
      that rotates the flux lines, so D_dir grows only once the iron saturates.

  (2) Even a hard FLAT cut is already nearly excitation-invariant (sub-degree direction
      drift, D_dir < 1e-2 rad) because the high-mu pole face stays equipotential
      (gap-reluctance robustness); the residual drift is pole-tip-corner-dominated.

  (3) The end chamfer that relieves the corner keeps the flux lines invariant several
      times DEEPER into saturation (D_dir several-fold smaller at the saturated drive),
      and it does so while dropping the corner kappa -- the same corner-relief lever as
      bending_endpack_saturation_opt.py.

ngsolve only (Optuna if present, else the built-in grid + local refine).
"""
import sys
from pathlib import Path

import pytest

# The compute helper lives beside the docs notebook (docs/clebsch_hodograph).
DOCDIR = Path(__file__).resolve().parents[2] / "docs" / "clebsch_hodograph"
sys.path.insert(0, str(DOCDIR))


@pytest.mark.slow
def test_excitation_invariant_field():
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import excitation_invariant_field as ei

    r = ei.optimize(trials=8, order=2, seed=0)
    f, o = r["flat_cut"], r["optimized"]

    # (2) even the FLAT cut is already nearly excitation-invariant, but with a real,
    # positive drift under saturation (the corner it saturates).
    assert 3e-4 < f["D_dir"] < 3e-3, f                    # ~1.5 mrad (sub-degree)
    assert r["flat_flux_lines_already_invariant"] is True, r

    # (3) the corner-relieving end chamfer keeps the flux lines invariant deeper into
    # saturation: several-fold smaller direction drift at the saturated drive.
    assert o["D_dir"] < f["D_dir"], (f, o)
    assert o["D_dir"] < 6e-4, o                           # ~0.2 mrad
    assert r["invariance_factor"] > 2.0, r                # measured ~6-7x
    assert r["flux_lines_more_invariant"] is True, r

    # the optimizer's corner-relief lever agrees with the kappa lever: the optimized
    # end has a LOWER pole-tip corner concentration than the flat cut.
    assert o["corner_kappa_sat"] < f["corner_kappa_sat"], (f, o)
    assert 1.4 < f["corner_kappa_sat"] < 2.6, f
    assert o["corner_kappa_sat"] < 1.35, o

    # the optimized chamfer is a real, positive, in-range end profile.
    assert 0.001 < o["depth_m"] < ei.G2, o
    assert 0.4 <= o["exponent"] <= 3.0, o
    assert r["n_evals"] >= 8, r
    # the iron actually saturated at the top drive (so the drift is a real effect).
    assert o["mur_sat"] < ei.MUR0, o


@pytest.mark.slow
def test_excitation_invariant_linear_control_is_zero():
    """LINEARITY => invariance: with mu forced constant, the flux-line direction drift
    is 0 at every drive; with saturation ON, it grows monotonically with the drive."""
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import excitation_invariant_field as ei

    c = ei.invariance_curve(0.0, 1.0, drives=(0.15, 0.75, 1.70), order=2)

    # linear control: identical flux lines at every excitation (drive-independent).
    assert c["linear_control_max"] < 1e-5, c
    # saturated: drift starts at ~0 (low drive) and grows monotonically with excitation.
    ds = c["D_dir_saturated"]
    assert ds[0] < 1e-6, ds                               # reference drive
    assert ds[-1] > ds[0], ds
    assert ds[-1] > ds[1] > ds[0], ds                    # monotone increasing
    assert ds[-1] > 3e-4, ds                             # a real, measurable drift
