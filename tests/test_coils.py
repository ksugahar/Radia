"""Filament current-loop helpers (radia.coils): circular_loop / helmholtz_pair,
validated against the single-loop and Helmholtz on-axis closed forms. Run via
PowerShell (radia C-ext needs the pwsh launch)."""
import math
import os
import sys

_RADIA_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src", "radia"))
if _RADIA_SRC not in sys.path:
    sys.path.insert(0, _RADIA_SRC)

import radia as rad
import coils

R, I = 0.05, 10.0          # 50 mm loop, 10 A


def main():
    # 1) single circular loop on-axis vs mu0 I R^2 / (2 (R^2+z^2)^{3/2})
    loop = coils.circular_loop([0, 0, 0], R, I, axis="z", nseg=96)
    worst = 0.0
    for z in [0.0, 0.03, 0.08]:
        b = rad.Fld(loop, 'bz', [0, 0, z])
        bref = coils.loop_axial_field(R, I, z)
        e = abs(b - bref) / bref
        worst = max(worst, e)
        print(f"  loop z={z:.3f}  Bz={b:.6e}  ref={bref:.6e}  err={e:.2e}")
    assert worst < 5e-4, f"loop axial off by {worst:.2e}"

    # convergence in nseg (faceting ~ 1/nseg^2)
    e48 = abs(rad.Fld(coils.circular_loop([0, 0, 0], R, I, nseg=48), 'bz', [0, 0, 0])
              - coils.loop_axial_field(R, I, 0)) / coils.loop_axial_field(R, I, 0)
    e96 = abs(rad.Fld(loop, 'bz', [0, 0, 0]) - coils.loop_axial_field(R, I, 0)) \
        / coils.loop_axial_field(R, I, 0)
    print(f"  faceting convergence: nseg48 {e48:.2e} -> nseg96 {e96:.2e}")
    assert e96 < e48, "loop faceting did not converge in nseg"

    # 2) Helmholtz pair (separation = radius) centre vs (4/5)^{3/2} mu0 I / R
    hp = coils.helmholtz_pair([0, 0, 0], R, I, axis="z", nseg=96)
    bc = rad.Fld(hp, 'bz', [0, 0, 0])
    bref = coils.helmholtz_center_field(R, I)
    eh = abs(bc - bref) / bref
    print(f"  helmholtz centre Bz={bc:.6e}  ref={bref:.6e}  err={eh:.2e}")
    assert eh < 5e-4, f"helmholtz centre off by {eh:.2e}"

    # 3) axis generality: x-axis loop reads on its own axis
    lx = coils.circular_loop([0, 0, 0], R, I, axis="x", nseg=96)
    bx = rad.Fld(lx, 'bx', [0, 0, 0])
    ex = abs(bx - coils.loop_axial_field(R, I, 0)) / coils.loop_axial_field(R, I, 0)
    print(f"  x-axis loop centre Bx={bx:.6e}  err={ex:.2e}")
    assert ex < 5e-4, f"x-axis loop off by {ex:.2e}"

    # 4) finite solenoid on-axis vs the thin finite-solenoid closed form.
    Rs, Ls, Nt, Is = 0.02, 0.04, 100, 1.0
    sol = coils.solenoid([0, 0, 0], Rs, Ls, Is, Nt, axis="z", nseg=48)
    worst_s = 0.0
    for z in (0.0, 0.01, 0.03):
        b = rad.Fld(sol, 'bz', [0, 0, z])
        bref = coils.solenoid_axial_field(Rs, Ls, Nt, Is, z)
        e = abs(b - bref) / abs(bref)
        worst_s = max(worst_s, e)
        print(f"  solenoid z={z:.3f}  Bz={b:.6e}  ref={bref:.6e}  err={e:.2e}")
    assert worst_s < 1e-2, f"solenoid axial off by {worst_s:.2e}"

    print("\n[OK] radia.coils circular_loop + helmholtz_pair + solenoid validated vs the "
          "loop / Helmholtz / finite-solenoid closed forms (<0.05% loops, <1% solenoid).")


if __name__ == "__main__":
    main()
