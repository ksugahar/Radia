"""FEMM-parity capability map for the radia-ngsolve FEM path.

Documents the FEMM (Finite Element Method Magnetics, D. Meeker) analyses that
are reproduced as EXECUTABLE + TESTED NGSolve capability in this package, with
the function to call, the analytical benchmark each was validated against, and
the regression test that locks it in. Every entry is validated to <2 % (most
<0.2 %) against a closed-form / exact reference -- "real code + tests", not
accumulated benchmark numbers.

The MCP server exposes this via femm_parity_documentation(topic=...). Topics:
overview, matrix, magnetics, scalar (electrostatic/heat/current), api, validation.
"""

FEMM_OVERVIEW = """\
# FEMM-parity on NGSolve -- executable, tested

FEMM solves 2D PLANAR and AXISYMMETRIC problems in four physics: magnetics,
electrostatics, heat flow, DC current flow. radia-ngsolve now reproduces the
core of all four as standard-NGSolve capability (the axisymmetric MAGNETIC
A_phi case uses the `axihenrotte` Q-element to kill the 1/r axis singularity;
everything else is plain H1, since planar magnetics and all scalar potentials
have no axis singularity).

Design rule learned: build CAPABILITY (a tested solver function) not a number.
Each analysis below ships a ``solve_*`` / extractor function in solve.py /
force.py / scalar_fem2d.py plus a regression test asserting agreement with an
analytical reference.
"""

FEMM_MATRIX = """\
# Capability matrix (planar / axisymmetric)

| FEMM analysis          | planar | axisym | reference & error                   |
|------------------------|:------:|:------:|-------------------------------------|
| magnetostatic (A_z/A_phi)| yes  | yes*   | sphere/cyl -0.05%                   |
| permanent magnet        | yes   | yes*   | sphere/cyl -0.05%                   |
| nonlinear B-H (Picard)  | yes   | yes*   | Ampere x B-H  <0.05% (planar)       |
| eddy current / AC       | yes   | yes*   | Kelvin Rac +0.07%; Cu-disk 0.27%    |
| current-driven circuit  | yes   | yes*   | net-current constraint exact        |
| voltage-driven circuit  | yes   | yes*   | Z=V/I  +0.07% (planar)              |
| force & torque (Maxwell)| yes   | (3D)   | two-wire mu0 I^2/2pi d  1.6%        |
| inductance (energy)     | yes   | yes*   | L_int = mu0/8pi  -0.06% (planar)    |
| electrostatics (csolv)  | yes   | yes    | coaxial/sphere -0.19/-0.15%         |
| heat flow (hsolv)       | yes   | yes    | coaxial/sphere -0.19/-0.15%         |
| DC current flow         | yes   | (same) | coaxial conductance -0.19%          |

(*) axisymmetric magnetics / eddy / nonlinear use H1Henrotte (axihenrotte FESpace)
-- see `axifemm_documentation`. "(3D)" = available via the 3D eggshell force.
API: solve_axi_magnetostatic / solve_axi_eddy / solve_axi_magnetostatic_nonlinear
in solve.py; inductance_axi / ohmic_loss_axi in force.py.
"""

FEMM_MAGNETICS = """\
# Magnetics API (FEMM prob1-4big analogs -- planar AND axisymmetric)

## PLANAR (Cartesian A_z, standard H1)

```python
from radia_mcp.radia_ngsolve.solve import (
    reluctivity, solve_planar_magnetostatic, planar_magnet_source,
    solve_planar_magnetostatic_nonlinear, solve_planar_eddy)
from radia_mcp.radia_ngsolve.force import (
    eggshell_force_2d, eggshell_torque_2d, inductance_2d, ohmic_loss_2d)

nu = reluctivity(mesh, {"iron": 1000})                 # nu = nu0/mu_r per material

# static + permanent magnet (FEMM prob1big planar)
A = solve_planar_magnetostatic(mesh, nu, Jz=Jz_cf,
        magnets={"mag": (Hc, 90.0)})                   # B = CF((grad(A)[1], -grad(A)[0]))

# nonlinear saturating iron (Picard on nu(|B|))
A = solve_planar_magnetostatic_nonlinear(mesh, nu_of_B, Jz=Jz_cf, relax=0.5)

# eddy / AC, current-driven (NumberSpace net-current constraint) -> compound gfu
gfu = solve_planar_eddy(mesh, nu, sigma, omega,
        driven_region="wire", total_current=I)
Az, Vc = gfu.components
# voltage-driven (prescribe axial E-field):
Az = solve_planar_eddy(mesh, nu, sigma, omega, applied_Ez=Vc)

# post: force/torque, planar inductance [H/m], AC loss [W/m]
Fx, Fy = eggshell_force_2d(B, mesh, center, r_in, r_out)
tau    = eggshell_torque_2d(B, mesh, center, r_in, r_out, pivot=(0,0))
L_pm   = inductance_2d(B, mesh, nu, current)           # H/m (per unit length)
P_pm   = ohmic_loss_2d(Ez, mesh, sigma, region="wire") # W/m ; Rac = 2P/|I|^2
```

## AXISYMMETRIC (A_phi, H1Henrotte FESpace)

```python
from radia_mcp.radia_ngsolve.solve import (
    solve_axi_magnetostatic, solve_axi_eddy, solve_axi_magnetostatic_nonlinear)
from radia_mcp.radia_ngsolve.force import inductance_axi, ohmic_loss_axi

# Axis (r=0) MUST be in dirichlet; B_z = grad(u)[0]+u/r, B_r = -grad(u)[1]
A = solve_axi_magnetostatic(mesh, nu,
        Jr=Jr_cf,                                      # phi-direction source [A/m^2]
        magnets={"mag": (Hc, 90.0)})                   # theta=90 = axial magnetization

# nonlinear B-H Picard (same nu_of_B interface as planar)
A = solve_axi_magnetostatic_nonlinear(mesh, nu_of_B, Jr=Jr_cf, relax=0.5)

# eddy / AC, current-driven (Vc = r*E_phi, NumberSpace):
gfu = solve_axi_eddy(mesh, nu, sigma, omega,
        driven_region="wire", total_current=I)         # I = 2pi int J_phi r dr dz
Az, Vc = gfu.components
# voltage-driven (Vc = r*E_phi prescribed):
Az = solve_axi_eddy(mesh, nu, sigma, omega, applied_Vc=Vc)

# post: 3D inductance [H], AC loss [W]
B    = CoefficientFunction((grad(Az)[0] + Az/x, -grad(Az)[1]))  # (Bz, Br)
L_3D = inductance_axi(B, mesh, nu, current)            # H (full torus)
P_3D = ohmic_loss_axi(-1j*omega*Az, mesh, sigma, region="wire") # W ; Rac = 2P/|I|^2
```

Conventions: axis(r=0) MUST be Dirichlet for H1Henrotte; magnet theta is from
the r-axis (theta=90 => +z, theta=0 => +r); Vc = r*E_phi [V/turn/rad] for axi
(vs E_z [V/m] for planar). The Cu-disk eddy reference τ₁ = 224.31 µs (BEM-Foster,
see axifemm_documentation "validation" topic).
"""

FEMM_SCALAR = """\
# Scalar potentials: electrostatics / heat / current flow (FEMM csolv/hsolv)

All three are the same operator  -div(c grad u) = f  -- one core, three wrappers.

```python
from radia_mcp.radia_ngsolve.scalar_fem2d import (
    EPS0, solve_electrostatic, capacitance, solve_thermal,
    solve_current_flow, conductance,           # planar
    solve_poisson_axi, capacitance_axi)        # axisymmetric (r-weighted)

# planar (per unit length)
V = solve_electrostatic(mesh, eps, {"hi": V0, "gnd": 0.0})   # E = -grad(V)
C = capacitance(V, mesh, eps, V0)                            # F/m
T = solve_thermal(mesh, k, {"hot": Th, "cold": Tc}, heat_source=q)
Vc= solve_current_flow(mesh, sigma, {"a": V0, "b": 0.0}); G = conductance(...)

# axisymmetric (full 3D quantity): boundary arcs split -> name via Nearest at
# TWO off-axis points (upper+lower); leave the r=0 axis unnamed (Neumann).
V = solve_poisson_axi(mesh, eps, {"inner": V0, "outer": 0.0})
C = capacitance_axi(V, mesh, eps, V0)                        # Farads (4 pi eps ab/(b-a))
```
"""

FEMM_VALIDATION = """\
# Validation references (regression tests, all green)

| test file (packages/.../tests)      | reference                          | error  |
|-------------------------------------|------------------------------------|--------|
| radia-axifemm/test_magnetized_sphere| B_in=2 mu0 mu_r Hc/(mu_r+2)        | -0.05% |
| radia-mcp/test_axi_magnetostatic    | same sphere via solve_axi_magnet.  | -0.05% |
| test_planar_magnet                  | B_in=mu0 mu_r Hc/(mu_r+1) (cyl)   | -0.05% |
| test_planar_force                   | two-wire mu0 I^2/(2 pi d)          |  1.6%  |
| test_planar_eddy                    | round-wire Kelvin Rac/Rdc (q=4)   | +0.07% |
| test_planar_eddy_voltage            | same, voltage-driven Z             | +0.07% |
| test_planar_nonlinear               | Ampere H=I/2 pi r, B=BH(H)        | <0.05% |
| test_planar_inductance              | L_int = mu0/(8 pi)                | -0.06% |
| test_scalar_fem2d                   | coaxial 2 pi c/ln(b/a) (eps/sig/k)| -0.19% |
| test_axi_scalar                     | sphere 4 pi c ab/(b-a) (eps,k)    | -0.15% |
| radia-axifemm/test_disk_eigenvalue  | Cu-disk tau_1 = 224.31 us (BEM)   |  0.27% |

Axi eddy (solve_axi_eddy) validated via Cu-disk tau_1; planar eddy validated via
Kelvin Rac. Gotchas: eggshell force needs a mesh-resolved band; nonlinear
point-B overshoots at interface (sample away); axisym arc boundaries split ->
Nearest at 2 points (upper + lower).
"""


def get_femm_parity_documentation(topic: str = "all") -> str:
    """Return FEMM-parity capability documentation for the requested topic."""
    sections = {
        "overview":   FEMM_OVERVIEW,
        "matrix":     FEMM_MATRIX,
        "magnetics":  FEMM_MAGNETICS,
        "scalar":     FEMM_SCALAR,
        "validation": FEMM_VALIDATION,
    }
    if topic == "all":
        return "\n\n".join(sections[k] for k in
                           ["overview", "matrix", "magnetics", "scalar", "validation"])
    if topic in sections:
        return sections[topic]
    raise ValueError(f"Unknown femm_parity topic '{topic}'. "
                     f"Available: {', '.join(sections)}, all")
