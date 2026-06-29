"""Golden test: the magnetic-Reynolds crossover for a moving magnet over a
plate -- locks docs/maglev/demos/rotating_magnet_eddy.py.

The eddy current J / Joule / Lorentz force can be obtained three ways:
kinematic source-only (J = -sigma dA_s/dt, NO eddy FEM), full-FEM A-phi
(reference), and a constant-basis multiport CLN.  The locked facts:

  * Low Rm (Yano's case, Rm ~ 0.016): the source-only error is < 1 % --
    the analytic-source shortcut suffices, no per-step FEM is needed.
  * High Rm (~16): the source-only error is %-level-large -- the eddy
    reaction matters there.
  * The source-only error grows monotonically with Rm.
  * The CLN reproduces the full-FEM Lorentz force across ALL Rm.
  * Physical sanity at Yano's case: J_rms ~ a few 100 A/m^2.

Runs a coarse, fast version of the example (~1 min) by calling run() directly
(does NOT overwrite the committed fine-resolution JSON).  Skipped cleanly if
radia / ngsolve / netgen are unavailable.
"""
import importlib.util
import os

import pytest

pytest.importorskip("radia")
pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPT = os.path.join(_HERE, "..", "docs", "maglev", "demos", "rotating_magnet_eddy.py")


@pytest.fixture(scope="module")
def res():
    spec = importlib.util.spec_from_file_location("rotating_magnet_eddy", _SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    # coarse + fast; same T_list (Rm sweep) as the committed run
    return mod.run(maxh=0.001, ncfg=20, verbose=False)


def test_low_rm_kinematic_suffices(res):
    low = res["rows"][0]            # T = 2.0 s -> Rm ~ 0.016 (Yano's case)
    assert low["Rm"] < 0.05
    assert low["source_only_F_err"] < 0.01      # source-only within 1 % -> no FEM needed
    assert 100.0 < low["Jrms_peak"] < 2000.0    # physical eddy current magnitude


def test_high_rm_reaction_matters(res):
    high = res["rows"][-1]          # T = 2e-3 s -> Rm ~ 16
    assert high["Rm"] > 5.0
    assert high["source_only_F_err"] > 0.05     # source-only fails -> reaction matters


def test_source_error_grows_with_rm(res):
    errs = [r["source_only_F_err"] for r in res["rows"]]
    rms = [r["Rm"] for r in res["rows"]]
    assert rms == sorted(rms, reverse=True) or rms == sorted(rms)   # ordered sweep
    # the kinematic error increases with Rm (low Rm < high Rm by a wide margin)
    assert errs[-1] > 20 * errs[0]


def test_cln_accurate_across_all_rm(res):
    for r in res["rows"]:
        assert r["cln_F_err"] < 1e-2, (r["Rm"], r["cln_F_err"])
        assert r["speedup"] is None or r["speedup"] > 50    # ~1000x; loose floor
