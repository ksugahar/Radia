"""test_filament_axis_orientation.py -- regression for the ObjFlmCur
axis-parallel segment bug (fixed 2026-05-29).

Bug: radTFlmLinCur's constructor special-cased a segment lying exactly
parallel to the global x-axis (LinVect.y == 0 && LinVect.z == 0) by
always using the identity rotation.  That is only correct for a +x
segment; a -x segment was silently evaluated as if its current flowed
in +x (wrong far endpoint AND wrong current direction), giving a field
that could be ~10x off.

The defect only fired when LinVect.y / LinVect.z were *bit-exactly*
zero, which happens for hand-specified coordinates (e.g. square loops),
not for arc chords whose endpoints differ by ~1 ULP.  That is why full
multi-segment rings converged but explicit square loops did not, and why
a loop at azimuth 90 deg broke while 270 deg (1 ULP off) did not.

These tests pin:
  1. on-axis Bz of a rotated filament loop is phi-independent;
  2. a single -x-aligned segment matches the analytical Biot-Savart
     reference (and equals its +x mirror);
  3. a full circular loop still matches the closed-form loop formula
     (guards against regressing the +x / general rotation paths).
"""
import math

import numpy as np
import pytest

import radia as rad
from radia.biot_savart import h_filament, h_segments

MU0 = 4.0e-7 * math.pi

# Tolerance rationale: the bug produced a sign flip + ~2-10x magnitude
# error (relative error of order 1-9).  The achievable agreement floor is
# ~1e-8: (a) the on-axis loop value is only phi-independent up to the
# floating-point rounding of each chord's cos/sin endpoints, and (b)
# Radia's internal segment formula and the numpy reference use different
# (analytically-equal) algebraic forms.  RTOL=1e-6 therefore sits two
# decades above the FP floor and six decades below the defect.
RTOL = 1e-6


def _quad_loop_pts(phi, a=0.15, dphi_deg=45.0, dz=0.025, zc=0.0):
    """Closed quadrilateral patch on a cylinder of radius a, centred at
    azimuth phi, spanning dphi_deg x dz."""
    dphi = math.radians(dphi_deg)
    plo, phih = phi - dphi / 2, phi + dphi / 2
    zlo, zhi = zc - dz / 2, zc + dz / 2
    c = [
        (a * math.cos(plo), a * math.sin(plo), zlo),
        (a * math.cos(phih), a * math.sin(phih), zlo),
        (a * math.cos(phih), a * math.sin(phih), zhi),
        (a * math.cos(plo), a * math.sin(plo), zhi),
    ]
    c.append(c[0])  # close the loop
    return c


# Azimuths chosen to include 90 and 270 deg, where the patch's constant-z
# chords are exactly +/- x aligned (the orientations that exercised the bug).
PHIS_DEG = [0, 30, 45, 90, 135, 180, 225, 270, 315]


@pytest.mark.basic
def test_filament_loop_onaxis_bz_phi_independent():
    """On-axis Bz of a rotated filament loop must not depend on phi.

    The observation point sits on the rotation (z) axis, so a radial-moment
    loop gives a phi-independent on-axis Bz.  Pre-fix, phi=90 deg was ~10x
    off; the spread below is dominated by floating-point differences between
    the per-segment rotation matrices, not by physics.
    """
    obs = [0.0, 0.0, 0.05]
    vals = []
    for pd in PHIS_DEG:
        pts = _quad_loop_pts(math.radians(pd))
        rad.UtiDelAll()
        h = rad.ObjFlmCur(pts, 1.0)
        vals.append(rad.Fld(h, 'b', obs)[2])
    rad.UtiDelAll()

    vals = np.array(vals)
    mean = vals.mean()
    assert mean != 0.0
    rel_spread = np.max(np.abs(vals - mean)) / abs(mean)
    assert rel_spread < RTOL, (
        f"on-axis Bz varies with phi: rel_spread={rel_spread:.3e}\n"
        + "\n".join(f"  phi={pd:>4}deg  Bz={v:+.12e}"
                    for pd, v in zip(PHIS_DEG, vals))
    )


@pytest.mark.basic
def test_filament_negx_segment_matches_reference():
    """A single segment exactly anti-parallel to +x must match the
    analytical Biot-Savart reference, and must equal its +x mirror."""
    obs = [0.0, 0.0, 0.05]
    # Endpoints with bit-exactly equal y and z -> LinVect is exactly -x,
    # which is precisely what triggered the identity-branch bug.
    start = [0.057402, 0.138582, -0.0125]
    end = [-0.057402, 0.138582, -0.0125]

    rad.UtiDelAll()
    h = rad.ObjFlmCur([start, end], 1.0)
    Bz = rad.Fld(h, 'b', obs)[2]
    rad.UtiDelAll()

    Bz_ref = MU0 * h_filament(start, end, obs, 1.0)[2]
    assert Bz_ref != 0.0
    assert abs(Bz - Bz_ref) <= RTOL * abs(Bz_ref), (
        f"-x segment field wrong: Radia Bz={Bz:+.6e}, ref={Bz_ref:+.6e}"
    )

    # The same wire traversed +x with reversed current must give the same field.
    rad.UtiDelAll()
    h_rev = rad.ObjFlmCur([end, start], -1.0)
    Bz_rev = rad.Fld(h_rev, 'b', obs)[2]
    rad.UtiDelAll()
    assert abs(Bz - Bz_rev) <= RTOL * abs(Bz_ref)


@pytest.mark.basic
def test_filament_axis_aligned_all_six_directions():
    """Field of a unit segment must be correct for +/-x, +/-y, +/-z
    (the +/-x cases go through the special-case branch)."""
    obs = [0.02, 0.03, 0.05]
    c = [0.01, -0.02, 0.015]   # segment midpoint, off all axes
    L = 0.08
    axes = {
        '+x': [L, 0, 0], '-x': [-L, 0, 0],
        '+y': [0, L, 0], '-y': [0, -L, 0],
        '+z': [0, 0, L], '-z': [0, 0, -L],
    }
    for name, d in axes.items():
        start = [c[i] - d[i] / 2 for i in range(3)]
        end = [c[i] + d[i] / 2 for i in range(3)]
        rad.UtiDelAll()
        h = rad.ObjFlmCur([start, end], 1.0)
        B = np.array(rad.Fld(h, 'b', obs))
        rad.UtiDelAll()
        B_ref = MU0 * h_filament(start, end, obs, 1.0)
        err = np.linalg.norm(B - B_ref)
        assert err <= RTOL * (np.linalg.norm(B_ref) + 1e-30), (
            f"axis {name}: Radia B={B}, ref={B_ref}, err={err:.3e}"
        )


@pytest.mark.basic
def test_filament_full_circle_matches_loop_formula():
    """Full circular loop converges to mu0*I*a^2 / (2*(a^2+d^2)^1.5).

    Guards against regressing the +x / general rotation paths while fixing
    the -x special case.
    """
    a = 0.05
    d = 0.04          # on-axis distance from loop plane to observation
    I = 1.0
    nseg = 360
    pts = [(a * math.cos(2 * math.pi * k / nseg),
            a * math.sin(2 * math.pi * k / nseg), 0.0)
           for k in range(nseg + 1)]
    obs = [0.0, 0.0, d]

    rad.UtiDelAll()
    h = rad.ObjFlmCur(pts, I)
    Bz = rad.Fld(h, 'b', obs)[2]
    rad.UtiDelAll()

    Bz_analytic = MU0 * I * a ** 2 / (2.0 * (a ** 2 + d ** 2) ** 1.5)
    assert abs(Bz - Bz_analytic) <= 1e-4 * abs(Bz_analytic), (
        f"loop Bz={Bz:+.6e}, analytic={Bz_analytic:+.6e}"
    )


if __name__ == "__main__":
    test_filament_loop_onaxis_bz_phi_independent()
    print("[OK] phi-independence")
    test_filament_negx_segment_matches_reference()
    print("[OK] -x segment matches reference")
    test_filament_axis_aligned_all_six_directions()
    print("[OK] all six axis directions")
    test_filament_full_circle_matches_loop_formula()
    print("[OK] full circle matches loop formula")
    print("ALL PASSED")
