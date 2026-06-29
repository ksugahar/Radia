"""(B) 1-turn coil via the stream function -- the single-wire limit.

RESEARCH example (track B foundation).  The stream function method designs
a CONTINUOUS current (its equal-current contours are the wires); a 1-turn
coil is the COARSEST discretization -- ONE contour = one wire -- so it
cannot reproduce a target field that needs a distributed current.  The
design task is therefore "the single BEST wire", and the honest deliverable
is the limit of what one turn can do vs the full stream-function current.

Demonstrated axisymmetrically (the canonical solenoid case): a target
uniform on-axis B_z over a central region is produced well by the full
stream function (the turn distribution on a cylinder), and the best SINGLE
loop is measurably worse -- the 1-turn limit.

On-axis field of a unit-current loop of radius R at z = z0 (mu0 = 1):
    B_z(z) = R^2 / (2 (R^2 + (z - z0)^2)^{3/2}).
The loop axial positions are the stream-function support; the loop currents
are the discretized stream function (turn density); the equal-current
contour count = the number of turns; 1 turn = the single best loop.
"""
import numpy as np


def _bz_loop(R, z0, z):
    """On-axis B_z at z of a unit-current loop radius R centred at z0."""
    return R * R / (2.0 * (R * R + (z - z0) ** 2) ** 1.5)


def solve(R=1.0, n_loops=12, loop_span=2.5, n_target=40, target_half=1.0):
    """Full stream-function current vs the best single (1-turn) loop."""
    z_loop = np.linspace(-loop_span, loop_span, n_loops)     # SF support
    z_obs = np.linspace(-target_half, target_half, n_target)
    target = np.ones_like(z_obs)                             # uniform B_z

    # --- full stream function: solve the turn-density (loop currents) ---
    G = np.array([[_bz_loop(R, z0, zz) for z0 in z_loop] for zz in z_obs])
    I, *_ = np.linalg.lstsq(G, target, rcond=None)
    bz_multi = G @ I
    err_multiturn = (np.linalg.norm(bz_multi - target)
                     / np.linalg.norm(target))

    # --- 1-turn coil: the single best loop (axial position + current) ---
    best_err, best_z0, best_I = np.inf, None, None
    for z0 in np.linspace(-loop_span, loop_span, 400):
        g = np.array([_bz_loop(R, z0, zz) for zz in z_obs])
        I1 = (g @ target) / (g @ g)                         # best current
        err = np.linalg.norm(I1 * g - target) / np.linalg.norm(target)
        if err < best_err:
            best_err, best_z0, best_I = err, z0, I1

    return {
        "err_multiturn": float(err_multiturn),
        "err_one_turn": float(best_err),
        "one_turn_radius": float(R),
        "one_turn_z0": float(best_z0),
        "one_turn_current": float(best_I),
        "n_loops_streamfunction": int(n_loops),
        "ratio_oneturn_over_multi": float(best_err / err_multiturn),
    }


def main():
    r = solve()
    print("(B) 1-turn coil via the stream function -- the single-wire limit\n")
    print(f"  full stream function ({r['n_loops_streamfunction']} turns):")
    print(f"    uniform-B_z error = {r['err_multiturn']:.3e}")
    print(f"  best 1-turn coil (single loop):")
    print(f"    R={r['one_turn_radius']:.2f}  z0={r['one_turn_z0']:.3f}  "
          f"I={r['one_turn_current']:.3f}")
    print(f"    uniform-B_z error = {r['err_one_turn']:.3e}  "
          f"({r['ratio_oneturn_over_multi']:.1f}x worse)\n")
    print("  => one turn is the coarsest stream-function realization; the")
    print("     design task is the single best wire, here quantified vs the")
    print("     full current.  Next: a non-axisymmetric target -> the single")
    print("     wire SHAPE (via the stream-function contour) + single-stroke.")


if __name__ == "__main__":
    main()
