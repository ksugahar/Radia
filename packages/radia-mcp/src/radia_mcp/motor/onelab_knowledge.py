"""ONELAB/GetDP electric-machine knowledge for radia_mcp.motor.

Distilled from the ONELAB ElectricMachines bundle (Sabariego-Gyselinck-Geuzaine,
ULiège-ULB; COPYING.txt is copyright-only with no explicit license text —
treat as academic/research code and cite the papers).  The bundle is
the canonical open-source reference for 2D rotating-machine FEA, and is
the basis for a 2025 student thesis that benchmarks the ONELAB torque
waveforms against an independent reference.

Reading this knowledge is the prerequisite for translating an ONELAB
GetDP model (.pro + .geo) into a radia_mcp.radia_ngsolve workflow.

Section map:
  overview       — bundle layout, motor types, file pattern
  groups         — Region / FunctionSpace / Constraint vocabulary
  formulation    — A-formulation (2D), magnetodynamics_a, moving-band
  magstadyn_source — VERIFIED weak form / 3 torque methods / slip-σ trick /
                   v×B motional term / Park dq0, read from .pro (2026-06-09)
  twod_corrections — 2D→3D corrections actually IN GetDP: axial length, (anti)
                   periodicity scaling, end-effect lumped circuit, skew multi-
                   slice, lamination-loss homogenisation, shape-sensitivity
  motor_types    — IM, PMSM, SRM, WFSM, shaded-pole reference table
  ngsolve_xlate  — GetDP → NGSolve material/boundary/source mapping
  analysis_modes — Static / FreqDomain (slip) / TimeDomain (transient)
  circuit        — Coupled circuit (im_circuit.pro pattern)
  post           — Torque / iron loss / phase-current / FFT post-pro
  liu_thesis     — 2025 thesis: ONELAB torque vs an independent reference
"""

from __future__ import annotations

# -----------------------------------------------------------------
# Section content
# -----------------------------------------------------------------

OVERVIEW = """\
## ONELAB/GetDP ElectricMachines — overview

Source: the ONELAB `ElectricMachines` bundle (open-source).
Maintainers: Ruth Sabariego, Johan Gyselinck, Christophe Geuzaine
            (ULiège-ULB, WBGreen FEDO, WIST3 ONELAB).
License: see the bundle's `COPYING.txt`.

Quick start: open `main.pro` in Gmsh; the GUI picks a model and a
case from the dropdown `Flag_AnalysisType` (Static / Time-Domain /
FrequencyDomain).

### Bundle layout (per machine `M`)

| File | Role |
|------|------|
| `M_data.geo` | All geometric and electrical parameters (DefineConstant) |
| `M.geo` | Top-level Gmsh script, includes stator + rotor |
| `M_stator.geo` | Stator geometry + Physical regions |
| `M_rotor.geo` | Rotor geometry + Physical regions |
| `M.pro` | Problem definition (GetDP): Group/Function/Constraint/...|
| `M_circuit.pro` | External circuit (3-phase Y/Δ, brush, ...) |
| `M.msh` | Mesh (generated) |
| `M.pre`, `M.res` | GetDP pre-processed problem + results |
| `infos_M.png` | Geometry preview |

### Available reference machines

| Bundle | Type | Reference |
|--------|------|-----------|
| `pmsm` | 8-pole PMSM (simplified) | Ferreira da Luz et al., IEEE Trans. Mag. 2002 |
| `pmsm_cbmag` | 8-pole PMSM (full geom) | Gyselinck et al., CBMag 2002 |
| `lomonova` | 8-pole in-wheel PMSM | Lomonova et al., COMPEL 2011 |
| `wfsm_4p` | 4-pole wound-field SM | Gyselinck et al., ELECTRIMACS 2002 |
| `t30` | TEAM 30a — IM | TEAM workshop benchmark |
| `im_3kW` | 4-pole 3 kW IM | Gyselinck PhD 2000 |
| `im` | 4-pole IM, broken bar | Guérard-Gyselinck-Lecomte 2005 |
| `srm` | Switched reluctance | Gyselinck-Geuzaine-Sabariego 2011 |
| `shaded` | Shaded-pole motor | (see `shaded.pro`) |

### Common shell

All bundles share the same shell (driver_*.py, main.pro, BH.pro,
machine_magstadyn_a.pro, perturb.geo).  In particular,
`machine_magstadyn_a.pro` is the **generic** A-formulation problem
used by every bundle — the per-machine `.pro` only specializes the
Region groups, the source-Region defaults, and the post-processing.

### How to read a model

1. Start from `M_data.geo` to find the rated quantities
   (NbrPhases, NbrPolePairs, Rotor_R[], Stator_R[], slot/tooth
   counts, rated frequency, current, ...).
2. Open `M.geo` + `M_stator.geo` + `M_rotor.geo` to identify the
   Physical regions (numbers and names).  These become the `Group`
   entries in the `.pro` and the `Region` labels in GetDP solves.
3. Open `M.pro` to see the source definitions and the `Group` block:
   stator winding regions (`Stator_Ind_Ap`, ...), rotor magnets or
   bars (`Rotor_Magnets`, `RotorC`), iron, air-gap, moving band.
4. The actual physics is in `machine_magstadyn_a.pro` (shared).
"""

GROUPS = """\
## GetDP Region / FunctionSpace / Constraint vocabulary

GetDP organizes a machine model around **Region groups** that the
`.pro` aliases under semantic names.  Below is the de-facto contract
followed by every ONELAB/ElectricMachines bundle (see
`machine_magstadyn_a.pro`).

### Source-region groups

```getdp
DefineGroup[
  Stator_Inds,       // all stator inductor regions (Ap, An, Bp, Bn, Cp, Cn)
  Stator_Ind_Ap, Stator_Ind_An,
  Stator_Ind_Bp, Stator_Ind_Bn,
  Stator_Ind_Cp, Stator_Ind_Cn,
  Rotor_Inds,        // wound-field rotor (wfsm) only
  RotorC,            // squirrel-cage rotor bars (im, t30)
  Rotor_Bars,
  Stator_Magnets, Rotor_Magnets,
  RotorC_Skin, RotorC_Core,    // segregated cage parts
];
```

### Domain groups

```getdp
Group {
  Stator      = Region[STATOR_IRON];
  StatorC     = Region[{ }];                   // conducting stator (eddy)
  Rotor       = Region[{ROTOR_IRON, ROTOR_AIR_GAP, Rotor_Magnets}];
  Rotor_Bnd_MB = Region[ROTOR_BND_MOVING_BAND];
  Stator_Bnd_MB = Region[STATOR_BND_MOVING_BAND];
  MB           = MovingBand2D[Rotor_Bnd_MB, Stator_Bnd_MB];
  Air_Gap_Mvg  = Region[ MB ];
  Domain_NonMag = Region[{ AIR, Air_Gap_Mvg, Stator_Inds, Rotor_Inds }];
  Domain_Mag    = Region[{ Stator, Rotor }];
  DomainCC_Mag  = Region[ Domain_Mag ];        // magnetic but non-conducting
  DomainC_Mag   = Region[ StatorC, RotorC ];   // conducting + magnetic
}
```

### FunctionSpace `Hcurl_a_Mag_2D`

`machine_magstadyn_a.pro` declares the 2D vector-potential function
space with Whitney 1-forms (Edge basis), and adds a constant-per-region
DOF for each conductor (the so-called "constant-of-integration" or
"floating-potential" DOF that ties the integral of `-σ ∂A/∂t` to the
external circuit).  This is the standard FEEC pattern (cohomology
class 1) and translates directly to NGSolve `HCurl` with one extra
H1-piecewise-constant space per conducting region.

### Constraints

Three categories:

1. **Boundary**: A = 0 on `Surface_Inf` (Dirichlet on machine outer
   surface for 2D scalar-A models; for true open-domain see Kelvin).
2. **Anti-periodicity**: A reduced sector (1 pole / pole-pair) is the
   norm; the constraint `Antiperiodicity_A` flips the sign of A on the
   far cut.  In NGSolve this maps to a `Periodic` space with a sign
   flip in the master-slave map.
3. **Source assignment**: phase currents `I_A(t)`, `I_B(t)`, `I_C(t)`
   are imposed via `Region[Stator_Ind_Ap/An]` constraints; in
   voltage-driven mode, the constraint is on the circuit DOF instead.
"""

FORMULATION = """\
## 2D magnetodynamic A-formulation (machine_magstadyn_a.pro)

In 2D Cartesian (`(x, y, 0)`) with rotation axis along `z`, the
electromagnetic field reduces to a scalar A_z(x, y, t).  The strong
form is:

  ∇ × (ν ∇ × A) + σ ∂A/∂t = J_s + ∇ × (ν B_r)         (Ampere–Faraday)

with `J_s` the externally-imposed conductor current (`Stator_Inds`,
`Rotor_Inds`), `B_r` the permanent-magnet remanence (Stator/Rotor
magnets), `σ` the electrical conductivity (cage bars, laminated iron
with homogenization), and `ν = 1/μ(|B|)` the (possibly non-linear)
magnetic reluctivity.

### Weak form (Whitney 1-form / edge element)

For each test 1-form `a' ∈ Hcurl_a_Mag_2D`:

  ∫_Ω  ν (∇ × A) · (∇ × a')  dΩ
+ ∫_Ω  σ ∂A/∂t · a'           dΩ          (only on DomainC_Mag)
- ∫_Ω  ν B_r · (∇ × a')      dΩ          (PM remanence)
- ∫_Ω  J_s · a'              dΩ           = 0

Plus the per-conductor constant-of-integration DOFs (one per region in
`StatorC ∪ RotorC`) that pin `∫_R (-σ ∂A/∂t) dΩ` to the externally
imposed loop current (or voltage, via the circuit).

### Moving band

The air-gap is split into a fixed stator-side annulus, a fixed
rotor-side annulus, and a thin **moving band** `MovingBand2D[…]` whose
mesh is regenerated at every time step so that nodes on the
rotor-side slide with the rotation angle θ(t).  This avoids the
remeshing cost of full rotor rotation and is the de-facto standard
since Davat-Ren-Lajoie-Mangin (1985) / Antunes-Bastos-Sadowski (2008).

### Non-linear iron (BH curve)

`BH.pro` provides:

```getdp
Function {
  nu_iron[] = nu_Iron[$1] ;   // ν as a function of |B|² ($1)
  dnu_iron[] = dnu_Iron[$1] ; // dν/d|B|²  for Newton
}
```

Newton iteration uses the Jacobian:

  K_jac = ∫ (ν ∇×a' · ∇×da + 2 dν/d|B|² (∇×A·∇×da) (∇×A·∇×a'))

Frequency-domain extension uses the **harmonic-balance** method
(Gyselinck-Vandevelde-Dular-Melkebeek-Oliveira-Kuo-Peng COMPEL 2003)
where the time-periodic A is expanded in N Fourier harmonics; the
coupled (N×NbrDOF) algebraic system is solved by Newton.

### Iron loss

This bundle carries only the **classical-eddy** lamination term
`(d²/12)·σ_eq·⟨(dB/dt)²⟩` (order-0 homogenisation, OFF by default since
`sigma_fe=0`; see `twod_corrections` #6).  Hysteresis + excess (the full Bertotti
decomposition), Jiles-Atherton and dynamic Preisach are NOT included here, and
there is **no `machine_post_a.pro`** in the bundle — plug such models in
a-posteriori from the time-stepped `B(t)` if you need them.
"""

MOTOR_TYPES = """\
## Motor-type reference table

| Bundle | NbrPolePairs | Stator Slots | Rotor type | Source |
|--------|-------------|--------------|------------|--------|
| `pmsm` | 4 (8-pole) | 24 | Surface PM (NdFeB) | Stator three-phase current |
| `pmsm_cbmag` | 4 | 24 | Surface PM | Stator three-phase current |
| `lomonova` | 4 | 24 (concentr.) or 48 (distr.) | Surface PM | Stator three-phase current |
| `wfsm_4p` | 2 (4-pole) | 36 | Wound rotor (DC) + damper | Stator 3-φ + rotor DC + damper cage |
| `t30` (TEAM 30a) | 1 | 36 | Squirrel cage | Stator 3-φ, slip-frequency |
| `im_3kW` | 2 | 36 | Aluminum cage | Stator 3-φ, slip-frequency |
| `im` | 2 | 36 | Cage with 1 broken bar | Stator 3-φ + bar fault |
| `srm` | 2 (8 stator / 6 rotor) | 8 | Salient passive iron | Phase-pulsed stator |
| `shaded` | 1 | 2 + shading ring | Squirrel cage | Single-phase + shading short-circuit |

### Permanent-magnet rotor (PMSM)

In `pmsm.pro` the rotor PMs are declared as

```getdp
Group {
  Rotor_Magnet_Np = Region[ROTOR_MAG_N_POSITIVE];
  Rotor_Magnet_Sp = Region[ROTOR_MAG_S_POSITIVE];
  ...
  Rotor_Magnets   = Region[{Rotor_Magnet_Np, ..., Rotor_Magnet_Sn}];
}
Function {
  br[Rotor_Magnet_Np] = br_value * Vector[ Cos[Th[]], Sin[Th[]], 0] ;
  br[Rotor_Magnet_Sp] = -br_value * Vector[ Cos[Th[]], Sin[Th[]], 0] ;
  ...
}
```

where `Th[]` is the per-magnet orientation angle.  For NGSolve
translation, the equivalent is a piecewise-constant
`CoefficientFunction` mapping each magnet-region label to a Br vector
in machine-fixed coordinates, then transforming by the rotation
`R(θ_rotor(t))` at each time step.

### Cage rotor (IM)

In `im_3kW.pro` the squirrel-cage bars are conducting AND magnetic
(if laminated lossy iron) AND mechanically connected via end-rings:

```getdp
Group {
  Rotor_Bars       = Region[{ROTOR_BAR_1, ..., ROTOR_BAR_N}];
  RotorC           = Region[Rotor_Bars];     // conducting
  RotorCC          = Region[ROTOR_IRON];      // non-conducting iron
}
```

The end-ring resistance / inductance is handled via the
`im_3kW_circuit.pro` external circuit, not in the FE domain (2D
limitation).

### Wound-field SM (WFSM)

`wfsm_4p` has rotor field winding (DC) + damper cage:

```getdp
Group {
  Rotor_Ind_F_p = Region[ROTOR_FIELD_POSITIVE];
  Rotor_Ind_F_n = Region[ROTOR_FIELD_NEGATIVE];
  Rotor_Ind     = Region[{Rotor_Ind_F_p, Rotor_Ind_F_n}];
  RotorC        = Region[{Rotor_Damper_Bar_1, ...}];  // damper cage
}
```

### SRM

`srm.pro` is purely **doubly salient** — no PM, no rotor coil, no
cage.  Phases are excited in sequence and the torque comes entirely
from the reluctance variation as the rotor poles align with the
energized stator poles.  Highly non-linear ν(B), strongly position-
dependent inductance.
"""

NGSOLVE_XLATE = """\
## ONELAB GetDP → NGSolve translation patterns

Goal: take a GetDP `.pro` + `.geo` ONELAB model and re-express it in
the radia_mcp.radia_ngsolve workflow (NGSolve `HCurl` for A, optional
`Periodic`/Kelvin for open domain, NGSolve `TimeStepping` for
transient).

### 1. Geometry / mesh

| GetDP | NGSolve |
|-------|---------|
| Gmsh `.geo` Physical Surfaces (2D) | Netgen `OCCGeometry` 2D faces with `Material("name")` |
| Gmsh Physical Lines | Netgen `face.edges` with `name=` for `Boundaries("name")` |
| Anti-periodicity cut | `mesh.Identify(left, right, IdentificationType.PERIODIC)` with sign-flip |

For a 2D motor the **fastest path** is:
1. Run Gmsh in batch (`gmsh -2 M.geo -o M.msh`) on the ONELAB `.geo`.
2. Convert .msh → .vol via NGSolve `ReadGmsh`.
3. Re-label Physical groups with `mesh.SetMaterial(idx, "stator_iron")`
   etc. *or* edit the .msh PhysicalNames before reading.

For a 3D motor (rare in ONELAB), use Cubit or Netgen OCC directly.

### 2. Function space

```python
from ngsolve import HCurl, H1, Periodic, FESpace
mesh = Mesh("M.vol")
fes_A = HCurl(mesh, order=2, dirichlet="surface_inf")
# Per-conductor 'constant of integration' DOFs:
fes_I = NumberSpace(mesh, definedon="rotor_bar_.*")
fes = FESpace([fes_A, fes_I])
```

For anti-periodic 1-pole sector:
```python
fes_A = Periodic(HCurl(mesh, order=2, dirichlet="surface_inf"),
                 phase=-1.0)
```

### 3. Material function

```python
nu_air  = 1.0 / 4e-7 / pi
nu_iron = CoefficientFunction([1.0 / mu_table(BH_curve_file)])
nu = mesh.MaterialCF({"air": nu_air, "stator_iron": nu_iron, ...})
sigma = mesh.MaterialCF({"rotor_bar": 3.5e7, "stator_iron": 2e6, ...},
                        default=0.0)
Br = mesh.MaterialCF({"magnet_n_p": (br_val * cos(th), br_val * sin(th), 0),
                      ...})
```

For non-linear ν(|B|), use NGSolve's nonlinear iteration with a
piecewise BH-curve `CoefficientFunction` (see Wakao examples in
`radia_mcp.feec` for the standard NGSolve Newton+line-search pattern).

### 4. Moving band

NGSolve does **not** ship a 2D moving-band primitive.  Two options:

**Option A (recommended, simple)**: re-mesh at each time step in the
air-gap annulus.  For a stator-fixed mesh `Stator_Air` and a rotor-fixed
mesh `Rotor_Air`, the air-gap annulus `MB` is meshed as a thin layer
rotated by θ_rotor(t).  Use NGSolve `BuildRefinementTree` + Netgen
`Identify` to stitch the layers.  Cost: ~1 s per step on a 50k-element
mesh.

**Option B (advanced)**: sliding-interface mortar coupling.  Maintain
a fixed stator side and a fixed rotor side, glue them via a Lagrange
multiplier on the air-gap circle.  See Bouillault-Buffa-Ren 1994 or
the NGSolve `MultiscaleHCurl` example.  Faster per-step but more code.

For a first port of an ONELAB bundle, use Option A.

### 5. Source assignment

Stator winding currents:
```python
J_s = mesh.MaterialCF({
    "stator_ind_Ap": ( 0, 0, +I_A(t) * N_turns / A_slot),
    "stator_ind_An": ( 0, 0, -I_A(t) * N_turns / A_slot),
    "stator_ind_Bp": ( 0, 0, +I_B(t) * N_turns / A_slot),
    ...
})
# Right-hand side:
f += J_s[2] * v_A * dx   # v_A is the A_z test function
```

For voltage-driven mode, add the per-coil DOF + the external
RL-circuit equation as a Schur-complement (or solve the coupled system
directly with NGSolve `BlockMatrix`).

### 6. Torque computation

Maxwell stress on the air-gap circle of radius `r_mid`:

```python
torque = Integrate(
    r_mid * (B_r * B_phi / mu_0) * 2 * pi * r_mid,   # 2D scalar
    mesh, definedon=mesh.Boundaries("air_gap_circle"))
```

Cross-check with the **virtual work** method (Coulomb 1983) — see
`radia_mcp.differential_forms.forces_knowledge`.
"""

ANALYSIS_MODES = """\
## Analysis modes

In `M_data.geo` and `main.pro`:

```getdp
DefineConstant[
  Flag_AnalysisType = {0, Choices{0="Static", 1="Time-Domain",
                                  2="Frequency-Domain"},
                       Name "Input/Analysis type"}
];
```

### Static (Flag_AnalysisType = 0)

DC field at fixed rotor position θ.  Used for:
- Cogging torque map (sweep θ)
- d/q inductance map (sweep I_d, I_q with FEA)
- No-load back-EMF (sweep θ at I = 0)

NGSolve equivalent: solve `K_nl(A) A = f(θ)` with Newton; no
∂A/∂t term, no σ term.

### Frequency-Domain (Flag_AnalysisType = 2)

Single-harmonic phasor solve at slip-frequency ω_slip = s · ω_e.
**Linear** material (ν constant per region) or **harmonic-balance**
non-linear (N harmonics, Gyselinck-Dular COMPEL 2003).

For an IM at slip `s`:
- Stator currents at frequency f_e
- Rotor currents (cage bars) induced at frequency s·f_e *in
  rotor frame*, which is f_e in stator frame — both express as
  same-frequency phasors.
- σ ∂A/∂t becomes jω σ A, contributing to a complex stiffness.

NGSolve equivalent: complex `HCurl` space, `jω σ` mass term,
`PARDISO` direct solve (or `BlockMatrix` for circuit coupling).

### Time-Domain (Flag_AnalysisType = 1)

Implicit Euler (or BDF2) time stepping.  Used for:
- Transient run-up (start from standstill)
- Inverter PWM excitation (carrier ~10 kHz on top of fundamental)
- Cogging + load torque waveforms
- Rotor mechanical equation `J dω/dt = T_e - T_load`

NGSolve equivalent: `BilinearForm` for K_nl + (σ/Δt) M; right-hand
side accumulates A_prev and previous PM remanence rotation.  See
`onelab_motor_to_ngsolve.py` example (planned, see TODO).
"""

CIRCUIT = """\
## External circuit (M_circuit.pro pattern)

GetDP couples the FE domain to lumped circuit elements via a
Modified Nodal Analysis (MNA) approach.

Typical 3-phase Y-connected stator with neutral N:

```getdp
Group {
  Input_Sources = Region[{IND_A, IND_B, IND_C}];
  Capacitors   = Region[{}];
  Resistors    = Region[{R_A_endwinding, R_B_endwinding, R_C_endwinding}];
  Inductors    = Region[{L_A_endwinding, L_B_endwinding, L_C_endwinding}];
}
Function {
  V[IND_A] = V_amp * Cos[2*Pi*Freq*$Time + phi_A];
  V[IND_B] = V_amp * Cos[2*Pi*Freq*$Time + phi_B];
  V[IND_C] = V_amp * Cos[2*Pi*Freq*$Time + phi_C];
}
Resistance[Region[R_A_endwinding]] = R_end_A;
```

**Endwinding** is treated as lumped (not 2D-FE-modeled) — this is
the standard 2D simplification.

For an IM cage, each rotor bar has its own circuit node, with the
end-ring inductances and resistances modeled in `im_circuit.pro`:

```getdp
Group {
  Bars      = Region[{BAR_1, ..., BAR_N}];
  Ring_R_a  = Region[RING_A];      // upper end-ring resistance
  Ring_R_b  = Region[RING_B];      // lower end-ring resistance
}
```

### NGSolve translation

Implement as a Schur-complement block on the per-conductor constant
DOFs:

```
[ K_FE          C_circ ] [ A      ]   [ J_s  ]
[ C_circ^T  Z_circ      ] [ I_circ ] = [ V_src]
```

`Z_circ` is the lumped impedance matrix (resistors, inductors,
capacitors), `C_circ` is the per-conductor coupling (turn count /
slot area).  Solve as a 2x2 block system; PARDISO handles it.
"""

POST = """\
## Post-processing recipes (machine_post_a.pro)

### Torque (per time step)

Maxwell stress tensor integrated on a circle in the air gap:

```getdp
// VERIFIED names (machine_magstadyn_a.pro): Torque_Maxwell + PostOperation Get_Torque
PostOperation Get_Torque UsingPost MagStaDyn_a_2D {
  Print[ Torque_Maxwell[Rotor_Airgap],  OnGlobal, StoreInVariable $Trotor  ];
  Print[ Torque_Maxwell[Stator_Airgap], OnGlobal, StoreInVariable $Tstator ];
}   // rotor- and stator-side air-gap integrals must agree (sanity check)
// Torque_Maxwell = CompZ[ XYZ[] cross (T_max[{d a}]*XYZ[]) ] * 2*Pi*AxialLength/SurfaceArea[]
```

The integrand `T_max = (1/μ₀) [B (B·n) - (1/2) |B|² n]` is **volume-averaged over
the whole air-gap annulus** by the `2*Pi*AxialLength/SurfaceArea[]` factor (Arkkio
band average) — NOT a single mid-gap contour, which is what makes it mesh-robust.
This is the **Arkkio formula** (Arkkio 1987 PhD).  (The earlier `Torque_T` /
`ElectricalForce_MS` sketch was illustrative — these are the actual GetDP names;
see `magstadyn_source`.)

Sanity check: torque computed on **any** circle inside the air-gap
should give the same value (cancellation of the radial-tangential
cross terms).  In NGSolve this is the standard

```python
T = Integrate(
    r_mid * (B_x * B_y - B_y * B_x) / mu_0,   # 2D z-component
    mesh,
    definedon=mesh.Boundaries("air_gap_mid_circle"))
```

### Iron loss (accurate to this bundle)

The bundle computes only the **classical-eddy** loss, via the order-0 lamination
homogenisation term `(d²/12)·σ_eq·⟨(dB/dt)²⟩` (see `twod_corrections` #6; the iron
is non-conducting in the solve, `sigma_fe=0` by default — the σ_eq is carried only
to evaluate this loss a-posteriori).  A FULL Bertotti decomposition (the hysteresis
`k_hys·f·B_pk^α` and excess `k_excess·f^1.5·B_pk^1.5` terms) is **NOT** in this
bundle, there are no `k_hys/k_excess` in any `M_data.geo`, and there is no
`machine_post_a.pro` — add those from an external IronLoss model if required.
`JouleLosses` (σ|∂A/∂t+∇V|²) covers the real conductors (cage bars, solid parts).

### Phase current / back-EMF

For voltage-driven simulations:
```getdp
Print[ I_stator_A[Stator_Ind_Ap], OnRegion Stator_Ind_Ap,
       File StrCat[Dir_Res,"I_A.txt"], Format Table ];
```

### Field harmonics (open-circuit)

For PMSM no-load back-EMF spectrum, perform FFT of `dΦ_A/dt` over
one electrical period.  The fundamental gives `K_e` (V/rpm); harmonics
3, 5, 7 reveal slot effects and PM step.
"""

LIU_THESIS = """\
## 2025 thesis — ONELAB torque vs an independent reference

### Topic

Cross-validation of ONELAB GetDP torque waveforms against an independent
torque reference on a **synchronous reluctance motor (SynRM)**.  The
result was used to bootstrap a topology-optimization study based on
Kishi-Wakao-Murata 2025 IEEJ (autoencoder + LS method, see
`topology_opt_knowledge`).

### Key documents

| File | Content |
|------|---------|
| torque-waveform comparison | Side-by-side T(θ) waveforms for the same SynRM geometry, time domain |
| senior theses (SynRM/IPM)  | Background reading on SynRM, IPM, motor FEA |
| Wakao group papers         | Autoencoder + LS topology-opt, Darwin TD model |

### Lesson learned

ONELAB's machine_magstadyn_a.pro using **A-formulation + moving
band** reproduces the independent torque reference to within plotting
accuracy on a production SynRM design.  This validates the path: GetDP
semantics are correct, the moving-band sliding interface is numerically
clean (no spurious tangential force on the air-gap mid-circle).  Any port
to NGSolve should be benchmarked against the **same** ONELAB run, which is
open and reproducible.

### Implications for radia_mcp.motor

1. The "correct" reference for an NGSolve port is the open, reproducible
   ONELAB result.  This removes one layer of licensing /
   reproducibility friction.
2. For SynRM specifically, the dominant validation metric is
   **average torque T_ave** + **ripple T_rip = (T_max-T_min)/T_ave**.
   These are the Pareto objectives used in
   `topology_opt_knowledge` (Wakao 2025).
3. The Liu Xinyao thesis is the bridge: it converts the ONELAB
   benchmark + the Wakao topology-opt method into a single end-to-end
   workflow.  Re-implement this in NGSolve to get a full open-source
   SynRM design loop.
"""

MAGSTADYN_SOURCE = r"""## machine_magstadyn_a.pro — VERIFIED weak form & techniques (read from source 2026-06-09)

The `formulation`/`post` sections above were written from general knowledge and
carried a few idealisations.  THIS section is transcribed directly from the actual
shared template `machine_magstadyn_a.pro` and the `im.pro` induction case, so the
GetDP term and PostOperation names below are exact.

### One formulation, three regimes (Flag_AnalysisType)

A single 2D vector-potential formulation `MagStaDyn_a_2D` (a = A_z) serves all of:
  0 = Static, 1 = Time-domain (`TimeLoopTheta`, implicit Euler θ=1), 2 = Frequency
  (`Type ComplexValue; Frequency Freq`).  No re-derivation between regimes — the
  time-derivative term `DtDof[...]` is exactly what becomes `jω` in the complex solve.

### Exact Galerkin terms (Equation block)

```getdp
Galerkin { [ nu[{d a}] * Dof{d a} , {d a} ] ; In Domain ; }             // curl-curl ν
Galerkin { JacNL[ dhdb_NL[{d a}] * Dof{d a} , {d a} ] ; In DomainNL ; } // Newton (nonlin Fe)
Galerkin { [ -nu[] * br[] , {d a} ] ; In DomainM ; }                    // PM remanence (-ν·Br)
Galerkin { DtDof[ sigma[] * Dof{a} , {a} ] ; In DomainC ; }             // EDDY  σ ∂A/∂t
Galerkin { [ sigma[] * Dof{ur} , {a} ] ; In DomainC ; }                // grad V (loop voltage)
Galerkin { [ -sigma[]*(Velocity[] *^ Dof{d a}) , {a} ] ; In DomainV ; } // MOTIONAL −σ(v×B)
Galerkin { [ -js[] , {a} ] ; In DomainS ; }                            // stranded source
```

plus circuit `GlobalTerm`s (U/I, Ub/Ib, Uz/Iz) and a `GlobalEquation Network`.
Reading notes:
- PM magnets are NOT a current source — they enter as `−ν·br` tested with `{d a}` (= ∇×a').
- `DtDof[σ Dof{a}]` is the eddy term; in frequency mode it becomes `jω σ` automatically.

### Motion — TWO coexisting techniques (the key takeaway for radia-motor)

1. **Moving band** (default): `MB = MovingBand2D[MovingBand_PhysicalNb,
   Stator_Bnd_MB, Rotor_Bnd_MB, SymmetryFactor]`.  Each step:
   `ChangeOfCoordinates[NodesOf[Rotor_Moving], RotatePZ[delta_theta[]]]` then
   `MeshMovingBand2D[MB]` re-connects ONLY the thin band layer.  Works WITH eddy
   currents and in the frequency domain — this is the open-source counterpart of the
   "eddy × moving gap" capability a magnetostatic air-gap element cannot give.
2. **v×B motional term** (optional, `Term_vxb`): keep the mesh fixed and add
   `−σ(Velocity[] *^ ∇×A)` on `DomainV = RotorC`, with
   `Velocity[] = wr * XYZ[] /\ Vector[0,0,-1]`.  Cheaper (no re-mesh) but valid only
   for the steady rotating-field regime.

### Frequency-domain induction — the slip-σ trick (im.pro:165)

```getdp
sigma[Rotor_Bars] = (Flag_AnalysisType==2 ? slip : 1.) * sigma_bars ;
wr = (Flag_AnalysisType==2) ? (1-slip)*2*Pi*Freq/NbrPolePairs : rpm/60*2*Pi ;
```

In frequency domain the rotor is NOT moved; the rotor-bar conductivity is scaled by
slip s so one complex phasor solve at the stator frequency captures the
slip-frequency rotor currents.  `slip` is swept (Onelab `Loop`) to trace the
torque–slip curve.  Standard steady-state IM phasor model.

### Torque — THREE methods actually provided (exact names)

- `Torque_Maxwell` (PostOperation `Get_Torque`): Maxwell stress, Arkkio-averaged
  over the WHOLE air-gap annulus:
  `CompZ[ XYZ[] /\ (T_max[{d a}]*XYZ[]) ] * 2*Pi*AxialLength/SurfaceArea[]`,
  integrated `In Domain`.  The `/SurfaceArea[]` (= 2π·r_avg band average) is what
  makes it mesh-robust — NOT a single contour.
- `Torque_vw` (virtual work):
  `CompZ[ 0.5*nu[]*XYZ[] /\ VirtualWork[{d a}] ]*AxialLength`
  `In ElementsOf[Rotor_Airgap, OnOneSideOf Rotor_Bnd_MB]` (GetDP built-in VirtualWork).
- `Torque_Maxwell_cplx` / `_cplx_2f`: frequency-domain mean torque + the 2ω
  pulsating component.
Both rotor-side `Torque_Maxwell[Rotor_Airgap]` and stator-side `[Stator_Airgap]`
are printed and must agree (built-in sanity check).

### Park / dq0 (Flag_ParkTransformation) — the hook for Ld/Lq + MTPA

`Pinv[]`/`P[]` tensors implement abc↔dq0; `Flux_dq0[] = P[]*Vector[Flux_a,b,c]`.
Drive `Idq0=[ID,IQ,I0]`, read `Flux_d/Flux_q` ⇒ Ld=Φ_d/I_d, Lq=Φ_q/I_q.  This is the
GetDP basis for the Ld-Lq / MTPA tests already in radia-motor.

### Losses

`JouleLosses = ∫ σ |∂A/∂t + ∇V|²` (and the `−v×B` form on DomainV), split into
`Rotor`, `Rotor_Fe` ⇒ rotor cage loss + iron eddy.  Exactly the rotor eddy-current
loss the moving-band + DtDof formulation gives and a magnetostatic harmonic gap
element cannot.

### radia-ngsolve port checklist (verified mapping)

| GetDP | NGSolve |
|-------|---------|
| a = A_z (Form1P) | `H1(order=p)` scalar A_z (2D curl-curl ⇒ grad-grad on A_z) |
| `nu[{d a}]*Dof{d a}` | `nu*grad(Az)*grad(vz)*dx` |
| `DtDof[σ Dof{a}]` | implicit-Euler `σ/Δt·(Az−Az_prev)`, or `1j*ω*σ` (complex H1) |
| `−nu·br` (PM) | `−nu*(Br_perp)·grad(vz)*dx`, Br a MaterialCF rotated by θ(t) |
| moving band | re-mesh thin annulus + `mesh.Identify(PERIODIC, sign=±1)` (antiperiodic) |
| slip-IM | scale rotor-bar σ by slip, complex solve (no motion) |
| Arkkio torque | band-average over air-gap annulus (matches radia force_validation) |

Benchmark target = the ONELAB run itself (open, reproducible).  Refs:
Ferreira-da-Luz/Dular/Sadowski/Geuzaine/Bastos, IEEE Trans. Mag. 2002 (periodicity
+ moving band); Gyselinck et al., IEEE Trans. Mag. 2003 (frequency-domain modelling
of motion).
"""

TWOD_CORRECTIONS = r"""## 2D-analysis corrections actually in ONELAB/GetDP (verified from source 2026-06-09)

A 2D (x, y) machine slice ignores axial / end / skew / lamination physics.  The
ONELAB bundle carries a stack of corrections that put the 2D answer back onto the
3D reality.  Each entry below is grounded in the actual files — not general knowledge.

### 1. Axial length & per-unit-depth scaling  (machine_magstadyn_a.pro)
Everything 2D is per metre of depth.  `AxialLength` (stack length L_stk) multiplies
flux linkage, torque, loss, circuit coupling.  E.g. `Flux` =
`SymmetryFactor*AxialLength*Idir[]*NbWires[]/SurfCoil[]*CompZ[{a}]`; inductor term
`DtDof[ AxialLength*NbWires/SurfCoil*Dof{a}, {ir} ]`.
→ radia-ngsolve: carry L_stk as a scalar on every global (flux, torque, R, L).

### 2. Symmetry / (anti)periodicity scaling  (Constraint MVP_2D + SymmetryFactor)
Model only `NbrPolesInModel` of `NbrPolesTot`.  The sector cut is linked by
`RotatePZ[-NbrPolesInModel*2*Pi/NbrPolesTot]` with `Coefficient ((#poles odd)?-1:1)`
(antiperiodic = −1, periodic = +1), and globals are scaled by `SymmetryFactor`
(`Dof{Ub}/SymmetryFactor`, `Flux` ×SymmetryFactor).
→ radia-ngsolve: `mesh.Identify(PERIODIC)` + sign; multiply globals by SymmetryFactor.

### 3. End-effects = lumped circuit, NOT in the FE domain  (im_circuit.pro)
The 2D slice has no end-windings / end-rings (a genuinely 3D region).  They are added
as LUMPED circuit elements:
- stator end-winding leakage reactance: `Inductance[L1,L2,L3] = Ls` (per phase);
- squirrel-cage end-ring: per-segment `Resistance[Rers_k]=R_endring_segment`,
  `Inductance[Lers_k]=L_endring_segment`, wired bar→bar as a ring network.
→ radia-ngsolve: add to the per-conductor circuit Schur block (Z_circ), not the mesh.

### 4. 3D winding-resistance factor  (Rb[], Factor_R_3DEffects, FillFactor_Winding)
`Rb[] = Factor_R_3DEffects * AxialLength * FillFactor_Winding * NbWires^2 / SurfCoil / sigma`.
`Factor_R_3DEffects` (>1) inflates the slot DC resistance to include the end-turn
copper length the 2D model cannot see; `FillFactor_Winding` (e.g. 0.96) is packing.

### 5. Skew = multi-slice  (generate_3d.geo: Model3D = Multi-slice)
Skewed slots cannot exist in one 2D plane.  `Model3D = {0=2D, 1=Multi-slice, 2=3D}`,
`NumSlices`, `SliceZOffset`: stack `NumSlices` z-translated copies of the 2D mesh;
the skew is realised by giving each slice's rotor a different angular offset, then the
torque / EMF are averaged over slices.  This is THE skew correction (Gyselinck
multi-slice).  → radia-ngsolve: solve N rotated sector copies (cheap, embarrassingly
parallel) and average — no true 3D mesh needed.

### 6. Lamination iron eddy = homogenisation of order 0  (im_3kW_data.geo:60-68)
Laminated iron is modelled as **non-conducting in the field solve** (`sigma_fe = 0`) —
resolving each lamina in 2D is impossible.  The classical eddy loss is recovered
a-posteriori by the homogenised term, quoted verbatim from the source:

    "value of sigma with homogenization of order zero ... it only makes sense when
     adding a term in the formulation:  d^2/12 * sigma * d_t b
     so that we can account for losses in the laminated domain, which is still
     non-conducting"

i.e. `P_lam = (d^2/12) * sigma_eq * <(dB/dt)^2>`, d = lamina thickness.  This is the
CLASSICAL eddy term ONLY (= Bertotti's middle term).  Hysteresis and excess are NOT
in this bundle and there is no `machine_post_a.pro`.  → radia-ngsolve (coreloss): use
B(t) post-hoc with the same `d^2/12·σ` term; add hyst/excess from a separate model.

### 7. Frequency-domain motion folded into a coefficient = slip-σ  (see magstadyn_source)
Scale rotor-bar σ by slip; one complex solve replaces a transient run-up.  A
"correction" in the sense of folding rotor motion into a material coefficient.

### 8. Numerical relaxation startup  (im.pro: Frelax[])
`Frelax[] = 0.5*(1 − Cos[Pi*t/Trelax])` smoothly ramps a nonlinear voltage source over
the first `NbTrelax` periods to avoid an inrush transient shock — a solver-robustness
correction, not physics.  → radia-ngsolve: same smooth ramp on the source.

### 9. Shape-sensitivity = mesh velocity field  (perturb.geo)
For shape / topology optimisation: perturb one geometric parameter by ε, remesh,
difference the node coordinates → a **mesh velocity field** `(x_pert − x)/ε` that
feeds the design derivative of torque/loss (semi-analytic shape gradient).  Ties
directly to the topology_opt work.  → radia-ngsolve: a velocity `GridFunction` used
with the shape-derivative of the torque functional.

### Adoption priority for radia-motor
(a) #1, #2, #4 are trivial scalings — adopt immediately.
(b) #3 end-effect circuit — needed for any voltage-driven / induction accuracy.
(c) #5 multi-slice skew — high value, cheap in NGSolve (N rotated solves averaged).
(d) #6 lamination loss — adopt the classical `d^2/12` term in radia coreloss, and be
    honest that hysteresis/excess need a separate model.
"""

# -----------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------

SECTIONS = {
    "overview": OVERVIEW,
    "groups": GROUPS,
    "formulation": FORMULATION,
    "magstadyn_source": MAGSTADYN_SOURCE,
    "twod_corrections": TWOD_CORRECTIONS,
    "motor_types": MOTOR_TYPES,
    "ngsolve_xlate": NGSOLVE_XLATE,
    "analysis_modes": ANALYSIS_MODES,
    "circuit": CIRCUIT,
    "post": POST,
    "liu_thesis": LIU_THESIS,
}


def get_onelab_knowledge(topic: str = "overview") -> str:
    """Return ONELAB / GetDP motor-FEA knowledge.

    Args:
        topic: section name (see SECTIONS keys) or "all" for everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(SECTIONS[k] for k in SECTIONS)
    if t not in SECTIONS:
        valid = ", ".join(sorted(SECTIONS))
        return f"Unknown topic {t!r}. Valid topics: {valid}\n"
    return SECTIONS[t]
