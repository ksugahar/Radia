"""Radiation BC (T^4) for the IH transient heat solver
(src/radia/panels/calc_heat_axisym.py, emissivity arg).

A workpiece with a UNIFORM surface heat flux q and radiative cooling
eps*sigma*(T^4 - T_ext^4) reaches an ISOTHERMAL steady state where every
surface point balances q = eps*sigma*(T_ss^4 - T_ext^4), i.e.

    T_ss = ( q/(eps*sigma) + T_ext_K^4 )^(1/4)   [Kelvin].

We run the transient axisymmetric solver to steady state and check T_max
against this closed form.  Validates the explicit-radiation RHS term and the
internal Celsius<->Kelvin handling.  (emissivity defaults to 0 = off, so the
existing convection-only golden tests are unaffected.)
"""
import math
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_PANELS = os.path.abspath(os.path.join(_HERE, "..", "..", "src", "radia", "panels"))
if _PANELS not in sys.path:
    sys.path.insert(0, _PANELS)

SIGMA_SB = 5.670374419e-8
R, H = 0.01, 0.02            # workpiece (r in [0,R], z in [-H/2,H/2])
Q = 3.0e4                    # uniform surface flux [W/m^2]
EPS = 0.8                    # emissivity
T_EXT = 20.0                 # ambient [degC]
VOL = r"C:\temp\heat_rad_axisym.vol"


def make_axisym_mesh():
    from netgen.occ import OCCGeometry, MoveTo
    from ngsolve import TaskManager
    wp = MoveTo(0, -H / 2).Rectangle(R, H).Face()
    wp.faces.name = "workpiece"
    wp.edges.name = "outer"
    os.makedirs(os.path.dirname(VOL), exist_ok=True)
    with TaskManager():
        ngm = OCCGeometry(wp, dim=2).GenerateMesh(maxh=0.0035)
        ngm.Save(VOL)


def main():
    make_axisym_mesh()
    from calc_heat_axisym import solve_heat_axisym

    res = solve_heat_axisym(
        VOL, material="steel",
        q_uniform=Q, emissivity=EPS, h_conv=0.0, t_ext=T_EXT, t_initial=T_EXT,
        dt=5.0, t_end=900.0, fes_order=1, linear_solver="sparsecholesky")

    assert "error" not in res, res
    T_ext_K = T_EXT + 273.15
    T_ss = (Q / (EPS * SIGMA_SB) + T_ext_K ** 4) ** 0.25 - 273.15
    T_max = res["T_max_C"]
    T_min = res["T_min_C"]
    err = (T_max - T_ss) / T_ss
    print(f"  steady-state radiation balance: T_max={T_max:.2f} C  "
          f"T_min={T_min:.2f} C  exact T_ss={T_ss:.2f} C  ({100*err:+.3f}%)")
    print(f"  (isothermal check: T_max-T_min = {T_max-T_min:.3f} C)")
    assert abs(err) < 0.02, f"radiation steady-state T off {100*err:.2f}% (>2%)"
    assert (T_max - T_min) < 0.02 * T_ss, "steady state should be ~isothermal"
    print(f"\n[OK] IH heat-solver radiation BC validated: steady T_max within "
          f"{100*abs(err):.2f}% of eps-sigma balance ({T_ss:.0f} C).")


if __name__ == "__main__":
    main()
