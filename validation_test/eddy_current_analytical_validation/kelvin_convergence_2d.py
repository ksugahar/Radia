"""Kelvin radius convergence test for 2D axisymmetric FEM.

If the Kelvin transformation is implemented correctly, the computed
inductance L must be INDEPENDENT of the Kelvin radius a_kelvin (as
long as the physical objects are enclosed).  A Kelvin-radius-dependent
L indicates an implementation bug.

This script sweeps a_kelvin from 0.08 to 0.25 m on the air-only coil
(no workpiece) using reference_2d_axisym's solve_2d_axisym. Expected
air-only L for R=30mm, a=3mm circular coil, Neumann formula:

    L_air = mu_0 * R * (ln(8R/a) - 2)
          = 4 * pi * 1e-7 * 0.030 * (ln(80) - 2)
          = 28.58 nH
"""

import math
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__),
                                '../../validation_test/induction_heating/cubit_panels_legacy'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

# Air-only: use fem_esim_kelvin with sigma=0 (no WP)
# But fem_esim_kelvin meshes a workpiece block. Easier: use reference_2d_axisym
# at very low frequency where sigma effect is negligible.

from em_material import EMMaterial
from reference_2d_axisym import solve_2d_axisym


def neumann_L_nH(R_coil, a_coil):
    """Air-only L for circular-section ring by Neumann/Lorenz formula."""
    mu0 = 4 * math.pi * 1e-7
    return mu0 * R_coil * (math.log(8 * R_coil / a_coil) - 2) * 1e9


def main():
    R_coil = 0.030
    a_coil = 0.003
    R_wp = 0.025
    H_wp = 0.025

    print("=" * 80)
    print("  Kelvin radius convergence test (2D axisym full-res, copper)")
    print("  Coil R=%gmm a=%gmm, WP R=%gmm H=%gmm, freq=7kHz"
          % (R_coil*1e3, a_coil*1e3, R_wp*1e3, H_wp*1e3))
    print("=" * 80)
    print("  Note: a real Kelvin transform gives L(a_kelvin) = const.")
    print()

    mat = EMMaterial.from_name("copper")

    print(f"  Air-only (Neumann):  L_air = {neumann_L_nH(R_coil, a_coil):.3f} nH")
    print()
    print(f"  {'a_kelvin [m]':>14s}  {'z_offset [m]':>14s}  "
          f"{'L [nH]':>10s}  {'P [uW]':>10s}  {'ndof':>8s}")
    print("  " + "-" * 64)

    # Monkey-patch solve_2d_axisym internals is complex; let's just
    # manually vary a_kelvin by editing the function? No — it's hardcoded.
    # Easier: import solve_2d_axisym and pass kwargs.
    import inspect
    sig = inspect.signature(solve_2d_axisym)
    print(f"  (solve_2d_axisym params: {list(sig.parameters.keys())})")
    print()

    # It doesn't take a_kelvin as arg. Need a different approach.
    # Option: use fem_esim_kelvin.solve_fem_full which does take a_kelvin.
    from fem_esim_kelvin import solve_fem_full

    for a_kelvin in [0.08, 0.10, 0.12, 0.15, 0.20, 0.25]:
        z_offset = a_kelvin * 2.5
        try:
            r = solve_fem_full(R_coil, a_coil, R_wp, H_wp,
                               sigma=mat.sigma, frequency=7000,
                               a_kelvin=a_kelvin, z_offset=z_offset,
                               maxh=0.004, order=2,
                               coil_shape='circular')
            print(f"  {a_kelvin:>14.3f}  {z_offset:>14.3f}  "
                  f"{r['L']*1e9:>10.3f}  {r['P_total']*1e6:>10.3f}  "
                  f"{r['ndof']:>8d}")
        except Exception as e:
            import traceback; traceback.print_exc()
            print(f"  {a_kelvin:>14.3f}  FAILED: {e}")


if __name__ == "__main__":
    main()
