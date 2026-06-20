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
| nonlinear magnet recoil | yes   | (*)    | self-consistent mu_r(B)  0.05%      |
| nonlinear B-H (Picard)  | yes   | yes*   | Ampere x B-H  <0.05% (planar)       |
| eddy current / AC       | yes   | yes*   | Kelvin Rac +0.07%; Cu-disk 0.27%    |
| current-driven circuit  | yes   | yes*   | net-current constraint exact        |
| voltage-driven circuit  | yes   | yes*   | Z=V/I  +0.07% (planar)              |
| force & torque (Maxwell)| yes   | yes    | 2-wire 1.6%; coaxial dM/dz 0.93%   |
| inductance (energy)     | yes   | yes*   | L_int = mu0/8pi  -0.06% (planar)    |
| electrostatics (csolv)  | yes   | yes    | coaxial/sphere -0.19/-0.15%         |
| heat flow (hsolv)       | yes   | yes    | coaxial/sphere -0.19/-0.15%         |
| DC current flow         | yes   | (same) | coaxial conductance -0.19%          |
| laminated steel (AC)    | yes   | (mat)  | complex mu_eff tanh; 1D-FE 0.00%    |
| laminated steel (DC)    | yes   | (mat)  | anisotropic; resolved stack 0.00%   |
| multi-conductor circuit | yes   | (—)    | series=2x single 0.34%; prox +20%   |
| nonlinear AC (eff. mu)  | yes   | (*)    | ->linear 3e-11; saturation sub-lin  |
| stranded / litz wire    | yes   | yes    | uniform-J DC field  0.14%           |
| AC current flow (cmplx) | yes   | (same) | coaxial Y=G+j w C  0.34%            |
| convection BC (Robin)   | yes   | yes    | slab L/k + 1/h  0.003%              |
| radiation BC (T^4)      | yes   | yes    | slab cond/rad balance  0.000%       |
| open-bdry x-check       | (—)   | yes    | vs stored open-bdry ref 0.78% of peak|
| axi forced eddy (BFI)   | (—)   | yes    | P_eddy(w): BFI+scipy; w^2 onset OK  |
| SIBC eddy force (3D)    | (3D)  | (3D)   | time-avg Maxwell on hole; 0.5*static|
| periodic/anti-periodic BC| yes  | (—)    | strip closed form  2e-7 / 2e-8     |

(*) axisymmetric magnetics / eddy / nonlinear use H1Henrotte (axihenrotte FESpace)
"(mat)" = a material model (anisotropic / complex-mu CF) usable in either path.
-- see `axifemm_documentation`. "(3D)" = available via the 3D eggshell / SIBC force.
API: solve_axi_magnetostatic / solve_axi_eddy / solve_axi_eddy_harmonic /
solve_axi_magnetostatic_nonlinear in solve.py; inductance_axi / ohmic_loss_axi /
maxwell_surface_force / maxwell_surface_force_harmonic in force.py.
Periodic / anti-periodic SECTOR BC (FEMM motor/machine cyclic symmetry):
periodic_h1(mesh, order, dirichlet, antiperiodic=False) wraps H1 with
ngsolve.Periodic -- the mesh must carry the netgen periodic identification (2D
SplineGeometry slave edge ``copy=master`` with a CCW loop, or OCC
``edge.Identify(.., IdentificationType.PERIODIC, trafo)``). antiperiodic=True =
phase -1 (A_slave=-A_master, half-period sector) -> complex space, take .real.

NOTE on the axi FORCED-harmonic eddy: the generic grad-weak-form solve_axi_eddy
canNOT assemble the complex (nu-stiffness + j w sigma-mass) system -- H1Henrotte's
AxiHenrotte DiffOps are real-coefficient only. The working path is
solve_axi_eddy_harmonic: assemble the REAL closed-form C++ integrators
AxiHenrotteStiffnessBFI(mu) + AxiHenrotteSigmaMassBFI(sigma) (the same ones behind
the validated Cu-disk tau_1 eigenvalue), symmetrise, form S=K+j w M in scipy and
solve. order=1 only (P2 AxiHenrotte has an axis-singularity NaN). Use the
matrix-based loss P_eddy=0.5 w^2 Re(x^H M x). The COMPLEX H1Henrotte field
value/gradient eval is now RELIABLE (fixed in src/ext/axifemm: the FESpace is
complex-capable and the DiffOps got complex CalcMatrix overloads; validated to
machine precision in test_axi_henrotte_complex_eval.py), so |A|/B post-processing
of the complex solution works -- this unblocks the nonlinear mu(|B|) Picard layer
(update mu per element from |B| at centroids, reassemble K each sweep, M fixed).

NOTE on the SIBC eddy FORCE: in calc_fem_kelvin the workpiece is a HOLE (no volume
mesh), so the only force handle is the time-averaged Maxwell stress over its
"sibc" surface: maxwell_surface_force_harmonic(curl(A_total), mesh, "sibc"). The
panel reports it as result["F_sibc"] (negated, since specialcf.normal points out
of the air mesh = into the hole). Reduces to 0.5*static maxwell_surface_force.
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

# CONVECTION / RADIATION BC for heat (FEMM hsolv mixed boundary):
#   -k dT/dn = h(T-Tinf)  and/or  eps*sigma_SB*(T^4 - Tinf^4)  [T in KELVIN]
T = solve_thermal(mesh, k, {"hot": T0},                      # fixed-T walls
                  convection={"surf": (h, Tinf)},            # Robin walls
                  radiation={"surf2": (eps, Tinf)})          # T^4 (Picard, Kelvin)
# AC current flow (complex sigma+j w eps); terminal admittance Y = G + j w C
from radia_mcp.radia_ngsolve.scalar_fem2d import solve_current_flow_ac, admittance_ac
V = solve_current_flow_ac(mesh, sigma, eps, omega, {"a": V0, "b": 0.0})
Y = admittance_ac(V, mesh, sigma, eps, omega, V0)            # S/m (per length)
```
"""

FEMM_LAMINATION = """\
# Laminated steel, multi-conductor circuits, FEMM cross-check

## Laminated steel (FEMM "Lamination & Wire Type" material)

A lamination stack is a MATERIAL MODEL, not a special solver -- two helpers in
solve.py build the effective property; feed it as the usual ``nu``.

```python
from radia_mcp.radia_ngsolve.solve import (
    laminated_mu_eff, laminated_reluctivity_tensor, solve_planar_eddy,
    solve_planar_magnetostatic)

# AC: complex effective permeability of in-plane laminations (eddy in each sheet)
#   mu_eff = mu0[ fill*mu_r*tanh(b)/b + (1-fill) ],  b = (d/2) sqrt(j w mu0 mu_r sigma)
mu_eff = laminated_mu_eff(mu_r=1000, sigma=2e6, omega=w, d_lam=0.5e-3, fill=0.95)
nu_lam = mesh.MaterialCF({"core": 1.0/mu_eff}, default=1.0/mu0)   # complex CF
Az = solve_planar_eddy(mesh, nu_lam, sigma=0, omega=w)  # sigma=0 in core: eddy is
#   already in Im(mu_eff) -- do NOT also mesh-resolve the sheets or set sigma there.

# DC: anisotropic reluctivity tensor (fill-factor): flux || sheets sees
#   mu_par = mu0(fill*mu_r+1-fill); across the stack mu_perp = mu0/(fill/mu_r+1-fill)
nu_t = laminated_reluctivity_tensor(mesh, {"core": (1000, 0.95)}, lam_normal="y")
A = solve_planar_magnetostatic(mesh, nu_t, Jz=Jz)   # matrix nu -> anisotropic form
```

## Multi-conductor circuit -- skin + PROXIMITY effect (FEMM "circuit")

```python
from radia_mcp.radia_ngsolve.solve import solve_planar_eddy_multi
# series: every conductor carries the SAME current (one Vc per conductor)
gfu = solve_planar_eddy_multi(mesh, nu, sigma, omega,
        conductors=["w0", "w1", "w2"], connection="series", total_current=I)
Az = gfu.components[0]
P = sum(ohmic_loss_2d(-1j*omega*Az + gfu.components[1+k], mesh, sigma, region=c)
        for k, c in enumerate(["w0","w1","w2"]))      # Rac_circuit = 2 P / |I|^2
# connection="parallel" -> one shared Vc, currents redistribute & sum to I.
```
Proximity (each conductor sees the others' AC field) is captured by the shared
A_z + per-conductor net-current constraints; it raises Rac above the isolated
sum (measured +20% for two wires at 2.4 a spacing).

## Stranded / litz conductor, nonlinear AC, axisymmetric force

```python
from radia_mcp.radia_ngsolve.solve import stranded_source, solve_planar_eddy_nonlinear
from radia_mcp.radia_ngsolve.force import eggshell_force_axi

# STRANDED / litz (FEMM "stranded" wire): imposed UNIFORM current, sigma=0 there
# -> Rac=Rdc, current stays uniform at any frequency (no skin/proximity in bundle)
Jz = stranded_source(mesh, {"bundle": I})                # Jz = I/area, int Jz = I
Az = solve_planar_eddy(mesh, nu, sigma=0, omega=w, Jz=Jz)   # sigma=0 in the bundle

# NONLINEAR AC (FEMM nonlinear harmonic, amplitude-based effective mu): Picard on
# nu(|B|) where |B| is the phasor amplitude. nu_of_B takes the SCALAR |B|.
Az = solve_planar_eddy_nonlinear(mesh, nu_of_B, sigma, omega, Jz=Jz, relax=0.4)

# AXISYMMETRIC net axial force [N] (full torus): meridional eggshell, 2*pi*r weight
B = CoefficientFunction((-grad(u)[1], grad(u)[0] + u/x))    # (B_r, B_z)
Fz = eggshell_force_axi(B, mesh, center=(rc, zc), r_inner=.., r_outer=..)
```

## Open-boundary cross-check (independent open-region reference)

``tests/test_femm_xcheck_kelvin_magnet.py`` loads a stored Kelvin-transform
(open-boundary) axisymmetric-magnet regression reference (B_z at r=7 mm) and
reproduces it with ``solve_axi_magnetostatic`` on a large truncated domain --
agreement 0.78 % of peak. A permanent magnet's field (mu_r=1) is scale-invariant,
so the SI solve matches the millimeter reference problem in tesla directly. The
reference is treated as an opaque stored data file.
"""

FEMM_VALIDATION = """\
# Validation references (regression tests, all green)

| test file (packages/.../tests)      | reference                          | error  |
|-------------------------------------|------------------------------------|--------|
| radia-core/axifem/test_magnetized_sphere| B_in=2 mu0 mu_r Hc/(mu_r+2)    | -0.05% |
| radia-mcp/test_axi_magnetostatic    | same sphere via solve_axi_magnet.  | -0.05% |
| test_planar_magnet                  | B_in=mu0 mu_r Hc/(mu_r+1) (cyl)   | -0.05% |
| test_planar_force                   | two-wire mu0 I^2/(2 pi d)          |  1.6%  |
| test_planar_eddy                    | round-wire Kelvin Rac/Rdc (q=4)   | +0.07% |
| test_planar_eddy_voltage            | same, voltage-driven Z             | +0.07% |
| test_planar_nonlinear               | Ampere H=I/2 pi r, B=BH(H)        | <0.05% |
| test_planar_inductance              | L_int = mu0/(8 pi)                | -0.06% |
| test_planar_periodic                | strip closed form (periodic+anti) | ~1e-7  |
| test_scalar_fem2d                   | coaxial 2 pi c/ln(b/a) (eps/sig/k)| -0.19% |
| test_axi_scalar                     | sphere 4 pi c ab/(b-a) (eps,k)    | -0.15% |
| radia-core/axifem/test_disk_eigenvalue  | Cu-disk tau_1 = 224.31 us (BEM)|  0.27% |
| test_laminated_steel (AC)           | mu_eff = mu_r tanh(b)/b vs 1D FE  |  0.00% |
| test_laminated_steel (DC)           | resolved stack vs anisotropic ten.|  0.00% |
| test_planar_proximity (single)      | round-wire Kelvin Rac/Rdc         | +0.07% |
| test_planar_proximity (2 far)       | series = 2x isolated Rac          | +0.34% |
| test_femm_xcheck_kelvin_magnet      | stored open-bdry ref B_z(r=7)     |  0.78% |
| test_scalar_fem2d_ext (Robin)       | slab T_L = T_inf+dT/(1+hL/k)      | +0.003%|
| test_scalar_fem2d_ext (AC I-flow)   | coaxial Y = 2pi(sig+jwe)/ln(b/a)  | -0.34% |
| test_stranded                       | uniform-J B(a/2)=mu0 I/(4 pi a)   | +0.14% |
| test_planar_eddy_nonlinear (lin)    | reduces to linear solve_planar_eddy| 3e-11 |
| test_axi_force                      | coaxial loops I1 I2 dM/dz (ellip.)| +0.93% |
| test_scalar_fem2d_ext (radiation)   | slab k(T0-TL)/L = eps sig(TL^4-..) | 0.000% |
| test_nonlinear_magnet               | B_op = mu0 mu_r(B_op)Hc/(mu_r+1)  | -0.05% |
| test_axi_eddy_harmonic              | P_eddy(w)>0, rising, w^2-onset    | physics|
| test_maxwell_surface_harmonic       | harmonic = 0.5*static (box face)  | machine|

Axi eddy (solve_axi_eddy) validated via Cu-disk tau_1; planar eddy validated via
Kelvin Rac. The FORCED-harmonic axi eddy (solve_axi_eddy_harmonic, BFI+scipy) is
validated by the eddy-loss physics P_eddy(w) (positive, rising with w, w^2 at low
f). Its matrix loss P_eddy=0.5 w^2 Re(x^H M x) is reliable, AND -- since the
complex H1Henrotte field eval was fixed (complex-capable FESpace + complex
CalcMatrix DiffOp overloads, validated to machine precision in
test_axi_henrotte_complex_eval.py) -- the complex |A|/B field post-processing of
the solution is now reliable too (this unblocks the nonlinear mu(|B|) Picard layer).
The harmonic Maxwell surface force (maxwell_surface_force_harmonic, the SIBC-hole
force handle) is validated by reduction to the independently cross-validated
static force (harmonic = 0.5*static, machine precision). The laminated DC/AC checks are exact
(uniform-field & 1D-eddy have closed forms representable to machine precision).
Gotchas: eggshell force needs a mesh-resolved band; nonlinear point-B overshoots
at interface (sample away); axisym arc boundaries split -> Nearest at 2 points
(upper + lower); for laminated AC set sigma=0 in the homogenized core (eddy is in
Im(mu_eff), not re-meshed); axi forced-harmonic eddy is order=1 only (P2 axis NaN).
"""


def get_femm_parity_documentation(topic: str = "all") -> str:
    """Return FEMM-parity capability documentation for the requested topic."""
    sections = {
        "overview":   FEMM_OVERVIEW,
        "matrix":     FEMM_MATRIX,
        "magnetics":  FEMM_MAGNETICS,
        "scalar":     FEMM_SCALAR,
        "lamination": FEMM_LAMINATION,
        "validation": FEMM_VALIDATION,
    }
    if topic == "all":
        return "\n\n".join(sections[k] for k in
                           ["overview", "matrix", "magnetics", "scalar",
                            "lamination", "validation"])
    if topic in sections:
        return sections[topic]
    raise ValueError(f"Unknown femm_parity topic '{topic}'. "
                     f"Available: {', '.join(sections)}, all")
