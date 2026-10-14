# -*- coding: utf-8 -*-
r"""End-to-end (coil / IH domain): the build123d helmholtz_pair makes a UNIFORM centre field.

Completes the geometry->field trio (accelerator: test_build123d_halbach_field; motor:
test_build123d_pmsm_field; coil/IH: here).  The coil mean radius and the separation are READ BACK from
the build123d `helmholtz_pair` geometry (z-centres of the two children; mean radius by inverting the
tube volume), then the on-axis Biot-Savart field is evaluated: at the Helmholtz spacing (separation =
mean radius) the centre field is FLAT (the 1st & 2nd derivatives vanish) and equals the textbook
(4/5)^(3/2) mu0 NI / R, while a non-Helmholtz spacing leaves a deep central dip.
"""
import math
import os
import sys

import numpy as np

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.archetypes import helmholtz_pair

MU0 = 4e-7 * math.pi
mm = 1e-3


def _coil_params(coil):
    """Read mean radius [m], separation [m], centre-z [m] back from the helmholtz_pair geometry."""
    zs = sorted(c.center().Z for c in coil.children)
    ch = coil.children[0]; bb = ch.bounding_box()
    r_out = bb.max.X; h = bb.max.Z - bb.min.Z
    r_in = math.sqrt(max(r_out**2 - ch.volume / (math.pi * h), 0.0))     # invert tube volume
    return 0.5 * (r_in + r_out) * mm, abs(zs[1] - zs[0]) * mm, 0.5 * (zs[0] + zs[1]) * mm


def _on_axis(coil, NI=1.0, frac=0.3, n=41):
    R, sep, zc = _coil_params(coil)
    loop = lambda z, z0: MU0 * NI * R * R / (2.0 * (R * R + (z - z0) ** 2) ** 1.5)
    zz = np.linspace(zc - frac * sep, zc + frac * sep, n)
    Bz = np.array([loop(z, zc - sep / 2) + loop(z, zc + sep / 2) for z in zz])
    return R, sep, Bz, Bz[n // 2]


def test_helmholtz_pair_uniform_centre_field():
    R_in, R_out, h = 40.0, 50.0, 8.0
    Rmean = 0.5 * (R_in + R_out)                                          # 45 mm
    hh = helmholtz_pair(R_in, R_out, h, separation=Rmean, name="hh")      # Helmholtz spacing
    R, sep, Bz, Bc = _on_axis(hh)
    assert abs(R - Rmean * mm) < 1e-4 and abs(sep - Rmean * mm) < 1e-4, "params read back from geometry"
    ripple = (Bz.max() - Bz.min()) / Bc
    Bideal = (4.0 / 5.0) ** 1.5 * MU0 * 1.0 / R                           # textbook Helmholtz centre field
    print(f"Helmholtz: R={R*1e3:.1f}mm sep={sep*1e3:.1f}mm, centre Bz={Bc:.4e} (analytic {Bideal:.4e}, "
          f"rel {abs(Bc-Bideal)/Bideal:.2e}), ripple over +/-30% gap = {100*ripple:.3f}%")
    assert ripple < 0.01, f"Helmholtz centre field must be uniform (got {100*ripple:.2f}% ripple)"
    assert abs(Bc - Bideal) / Bideal < 0.02, "centre field must match (4/5)^1.5 mu0 NI / R"


def test_non_helmholtz_spacing_has_central_dip():
    R_in, R_out, h = 40.0, 50.0, 8.0
    Rmean = 0.5 * (R_in + R_out)
    hh = helmholtz_pair(R_in, R_out, h, separation=2.5 * Rmean, name="hh")
    _, _, Bz, Bc = _on_axis(hh)
    ripple = (Bz.max() - Bz.min()) / Bc
    print(f"non-Helmholtz (sep=2.5R): ripple over +/-30% gap = {100*ripple:.1f}%")
    assert ripple > 0.2, "an over-wide pair must NOT be uniform (deep central dip)"


def main():
    test_helmholtz_pair_uniform_centre_field()
    test_non_helmholtz_spacing_has_central_dip()
    print("[OK] build123d helmholtz_pair: geometry read-back -> on-axis field uniform at the Helmholtz "
          "spacing (matches (4/5)^1.5 mu0 NI/R), central dip when over-spaced.")


if __name__ == "__main__":
    main()
