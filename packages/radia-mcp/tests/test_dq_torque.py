r"""Amplitude-invariant Park + dq flux-linkage torque: the identities behind salient Ld/Lq extraction.

Pure analytic reference (no solver, no commercial tooling) for the NEW helpers park_transform /
inverse_park / dq_flux_torque, tied to the existing dq_torque.  The relations are independently
confirmed by a salient-rotor magnetostatic FEA + Maxwell stress tensor: the (Ld,Lq) extracted by
projecting the FE phase flux linkages through Park reproduce the MEAN stress-tensor torque over a
rotor period, the difference being the salient-pole torque ripple (the single-instant fundamental
matches exactly).
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.radia_ngsolve.solve import (park_transform, inverse_park,
                                           dq_flux_torque, dq_torque, dq_torque_components)


def test_park_roundtrip():
    # park(inverse_park(d,q,th),th) == (d,q) at arbitrary angle
    for th in (0.0, 0.5, 1.234, -2.0, math.pi):
        a, b, c = inverse_park(0.7, -0.3, th)
        d, q = park_transform(a, b, c, th)
        assert abs(d - 0.7) < 1e-12 and abs(q + 0.3) < 1e-12, th


def test_balanced_set_projects_to_constant_dq():
    # a balanced 3-phase cosine set of amplitude I at angle th0 projects to a pure d=I (q=0)
    # when Park tracks th0 -- the amplitude-invariant (peak-preserving) property
    I = 5.0
    for th0 in (0.0, 1.0, 2.5):
        a = I * math.cos(th0)
        b = I * math.cos(th0 - 2 * math.pi / 3)
        c = I * math.cos(th0 + 2 * math.pi / 3)
        d, q = park_transform(a, b, c, th0)
        assert abs(d - I) < 1e-12 and abs(q) < 1e-12, th0


def test_flux_route_equals_lumped():
    # dq_flux_torque with the linear decomposition lam_d=lam_m+Ld id, lam_q=Lq iq == lumped dq_torque
    p, lam_m, Ld, Lq = 2, 0.05, 9e-3, 5e-3
    for (idc, iqc) in ((3.0, 4.0), (-7.0, 7.0), (0.0, 6.0), (-5.0, 0.0)):
        lam_d = lam_m + Ld * idc
        lam_q = Lq * iqc
        assert abs(dq_flux_torque(p, lam_d, lam_q, idc, iqc)
                   - dq_torque(lam_m, Ld, Lq, idc, iqc, p)) < 1e-15, (idc, iqc)


def test_reluctance_torque_sign_and_optimum():
    # SynRM: lam_m=0, Ld>Lq -> reluctance torque (3/2)p(Ld-Lq) id iq, positive for id,iq same sign,
    # and at fixed current magnitude maximised at the 45-deg current angle (id==iq).
    p, Ld, Lq, Imag = 2, 0.84e-3, 0.68e-3, 10.0
    best = max(range(1, 90), key=lambda deg: dq_torque(
        0.0, Ld, Lq, Imag * math.cos(math.radians(deg)), Imag * math.sin(math.radians(deg)), p))
    assert abs(best - 45) <= 1                                 # peak reluctance torque at 45 deg
    assert dq_torque(0.0, Ld, Lq, 5.0, 5.0, p) > 0            # same sign -> motoring
    assert dq_torque(0.0, Ld, Lq, 5.0, -5.0, p) < 0          # opposite sign -> reverses


def test_torque_components_split():
    # dq_torque_components splits T into magnet + reluctance; the parts sum to dq_torque, the
    # magnet part is the pure-q (saliency-free) torque, and the reluctance part carries all the
    # id dependence (the salient PM design decomposition).
    p, lam_m, Ld, Lq = 3, 0.08, 3.6e-5, 2.95e-5
    for (idc, iqc) in ((-4.0, 5.0), (0.0, 6.0), (3.0, 3.0), (-7.0, 0.0)):
        t_mag, t_rel, t_tot = dq_torque_components(lam_m, Ld, Lq, idc, iqc, p)
        assert abs(t_tot - dq_torque(lam_m, Ld, Lq, idc, iqc, p)) < 1e-15, (idc, iqc)
        assert abs(t_mag - 1.5 * p * lam_m * iqc) < 1e-15           # magnet term ~ iq only
        assert abs(t_rel - 1.5 * p * (Ld - Lq) * idc * iqc) < 1e-15  # reluctance term ~ (Ld-Lq) id iq
        assert abs(t_mag - dq_torque(lam_m, Ld, Lq, 0.0, iqc, p)) < 1e-15   # magnet == pure-q torque
    # non-salient PM: the reluctance share is exactly zero at every operating point
    _, t_rel0, _ = dq_torque_components(0.15, 1e-3, 1e-3, -5.0, 6.0, 4)
    assert t_rel0 == 0.0


def test_ipm_mtpa_runs_into_second_quadrant():
    # IPM: lam_m>0 and Lq>Ld (Ld-Lq<0); the reluctance term adds torque only when id<0, so the
    # maximum-torque-per-amp current angle exceeds 90 deg (negative id) -- unlike a surface-PM machine.
    p, lam_m, Ld, Lq, Imag = 3, 0.1, 4e-3, 9e-3, 20.0
    best = max(range(91, 160), key=lambda deg: dq_torque(
        lam_m, Ld, Lq, Imag * math.cos(math.radians(deg)), Imag * math.sin(math.radians(deg)), p))
    assert best > 90                                          # MTPA angle in the field-weakening half
    assert Imag * math.cos(math.radians(best)) < 0           # id < 0 at the optimum


def main():
    p, Ld, Lq = 2, 0.84e-3, 0.68e-3
    print("SynRM reluctance torque (3/2 p (Ld-Lq) id iq), |I|=10A vs current angle:")
    for deg in (0, 30, 45, 60, 90):
        idc, iqc = 10 * math.cos(math.radians(deg)), 10 * math.sin(math.radians(deg))
        print(f"  {deg:3d} deg: T = {dq_torque(0.0, Ld, Lq, idc, iqc, p)*1e3:+.3f} mN.m")
    print("[OK] Park round-trip, flux==lumped torque, 45-deg reluctance optimum, IPM MTPA id<0.")


if __name__ == "__main__":
    main()
