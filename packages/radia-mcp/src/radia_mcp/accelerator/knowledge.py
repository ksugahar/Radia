"""
Accelerator magnet design — end-pole chamfer analytical theory,
multipole analysis, rotating coil measurement, Radia case studies.
"""


# Authoritative topic enum for the dispatcher tool (wired into
# `<short>_topics()` via common.register_topics_tool).
TOPICS: dict[str, str] = {
    "end_pole": "Analytical chamfer design (Delferriere)",
    "kolkata": "Radia + TOSCA validation case study",
    "rotating_coil": "Multipole measurement + field reconstruction",
    "all": "return \"\n\n\".join([",
}

END_POLE_DESIGN = """
# Analytical end-pole chamfer for accelerator magnets
(Delferriere-de Menezes-Duperrier, CEA Saclay, SOLEIL design context)

## The problem

Accelerator magnets (dipoles, quadrupoles, sextupoles) need to deliver
a precise field profile along the beam axis with tight multipole
tolerances (~10⁻⁴ relative).  In the 2D cross-section, this is done by
careful pole-tip shaping.  But in 3D, the END of the magnet introduces:

- A fringe field region with mixed-multipole content
- An effective magnetic length that depends on end geometry
- Local saturation at the pole tip corners

A simple 45° chamfer (the historical default) is INADEQUATE for
modern light sources (SOLEIL, ESRF, NSLS-II) where tolerances reach 10⁻⁵.

## The analytical model

Start from the 3D magnetic scalar potential expansion:

    V(r, θ, z) = Σ_n (A z + B_n(z)) · r^n · (a_n cos nθ + b_n sin nθ)

For the desired multipole order n (n=1 dipole, n=2 quadrupole, n=3
sextupole), set L_m = L_f (magnetic length = iron length) at the
design field.  Then the longitudinal pole profile r(z) for z in the
end region II is:

    r(z) = ∆ · (1/2 - z/L_f)^(1/n)

where ∆ is a depth parameter, L_f is the iron length, and n is the
multipole order.

For a dipole (n=1):  r(z) = ∆ · (1/2 - z/L_f)   — linear taper
For a quadrupole (n=2): r(z) = ∆ · √(1/2 - z/L_f)   — square-root
For a sextupole (n=3): r(z) = ∆ · (1/2 - z/L_f)^(1/3)   — cube-root

## The depth ∆ determines magnetic length

For the dipole with half air-gap g:
    ∆ = L_f · (g·B_dyn/(μ₀·NI) - 1)

where B_dyn is the desired uniform-field value and NI is the
ampere-turns.

For the quadrupole with bore radius r_g and gradient G_dyn:
    ∆ = L_f · (r_g²·G_dyn/(2·μ₀·NI) - 1)

## Numerical verification (TOSCA)

Tested in TOSCA 3D for both dipole and quadrupole geometries:
- Sharp end (no chamfer): integrated multipoles A_3, A_5, A_7 ~ 10⁻³
- 45° chamfer: ~10⁻⁴
- Analytical profile (10-slope approximation): ~10⁻⁵ — 3 orders of
  magnitude improvement

The first slope is NOT 45° — it depends on ∆ and the multipole order.
Smooth profile via 10 slopes is a manufacturable approximation.

## Connection to Radia

Radia's strength is rapid evaluation of magnet geometries WITHOUT
iron meshing — perfect for parametric studies of end-pole ∆.

A Radia workflow:
1. Use the analytical formula r(z) = ∆ · (1/2 - z/L_f)^(1/n) to
   generate the pole tip surface
2. Use Radia ObjHexahedron / ObjPolyhedron with magnetization
   M_s from the saturation curve
3. Sweep over ∆ values, evaluate the integrated multipoles via
   rad.Fld() along beam trajectory
4. Choose ∆ that minimizes A_3 + A_5 + A_7 (or the relevant
   multipole budget)
5. Validate with TOSCA / NGSolve (using the radia-mcp.radia_ngsolve
   FEM pipeline)

This is exactly the workflow used at SOLEIL, ESRF, and Kindai's
accelerator partners.
"""

KOLKATA_CYCLOTRON = """
# Radia validation: Kolkata Superconducting Cyclotron case study
(Pradhan et al, VECC Kolkata, Cyclotrons 2007)

## The problem

Kolkata SCC has 2 superconducting coils (α-coil + β-coil) excited
to various current combinations.  Internal field can only be measured
up to r = 673 mm (extraction radius).  Beyond that, simulation must
extrapolate.

## Why Radia + TOSCA

The paper compares two simulation codes:

| Aspect | TOSCA (FEM) | Radia |
|---|---|---|
| Solver | Finite element | Boundary integral |
| Mesh | 450k elements needed | None |
| Geometry input | Complex (FEM mesh) | Simple (Mathematica) |
| CPU time | Long (1-2 hours) | Fast (few minutes) |
| Mathematica integration | None | NATIVE (RADIA exports to MATHEMATICA) |
| Saturation handling | Native | Needs trick (see below) |

The paper explicitly says: "RADIA requires geometry creation using
simple MATHEMATICA code, requires less CPU time of the solver."

This validates the radia-mcp ecosystem's Mathematica+Radia coupling —
it's not a new invention, it's the natural way RADIA has always been
used at major accelerator labs.

## The "saturated pole" trick (when FEM mesh is too expensive)

For complex SC cyclotron pole tips (shims, trim-coil inserts,
splittings):

1. Assume pole-tip surface is FULLY SATURATED at design field
2. Compute the field from saturated surface B_sat(r, θ) from ONE
   measured data set (calibration)
3. Compute the coil field B_coils(r, θ, I_α, I_β) analytically
4. The yoke contribution B_iron(r, θ, I_α, I_β) is computed via
   Radia (just the iron, ignoring detailed pole-tip features)
5. Total:  B_total = B_sat + B_coils + B_iron

This DECOMPOSITION trick avoids meshing the detailed pole-tip
geometry while preserving the dominant saturation physics.

## Results

Average iron field: |measured - simulated| < 0.5% across all coil
current combinations (extrapolation to higher currents than
measured).

3rd, 6th, 9th harmonics: agreement < 0.2%.

This validates Radia as a TRUSTWORTHY tool for SC cyclotron design,
including extrapolation beyond measured operating points.

## Connection to radia-mcp

The Kolkata workflow is exactly what radia-mcp.electromagnet supports:
- Build the coil with CoilBuilder
- Build the iron yoke with Radia primitives
- Apply BH curve with MatSatIsoTab
- Compute field via rad.Fld()
- Verify via Mathematica (radia-mcp.mathematica) symbolic identities
- Use TOSCA/NGSolve for cross-validation when needed

The Pradhan 2007 paper is THE blueprint for using Radia at an
accelerator facility.
"""

ROTATING_COIL_MEASUREMENT = """
# Rotating coil measurement and field reconstruction

The standard technique for measuring multipole content of accelerator
magnets.  Closely tied to Radia simulation for validation.

## The principle

A coil rotates around the magnet axis.  As it rotates, the magnetic
flux through it varies sinusoidally with frequency components:
    n = 1 → dipole
    n = 2 → quadrupole
    n = 3 → sextupole
    ...

Fourier-decompose the induced voltage to extract the multipole
amplitudes.

## What you get

For each multipole order n, you measure:
- The complex amplitude C_n (magnitude + phase)
- The "main" component (the design multipole)
- The "allowed" multipoles (with same symmetry as main)
- The "forbidden" multipoles (manufacturing errors)

## Local field reconstruction

A single rotating-coil measurement gives INTEGRATED multipoles
along the entire magnet length.  But for storage-ring optics, the
LOCAL field profile along z matters.

The reconstruction problem:
- Inputs: rotating-coil C_n(z) at multiple z positions OR
  a single integrated measurement
- Outputs: local B_x(x, y, z), B_y(x, y, z) field map

Approach 1: many rotating coils at different z (expensive)
Approach 2: 1 rotating coil + 3D Radia model fit to measurement
Approach 3: integrated + analytical end-pole model (Delferriere)

## Connection to radia-mcp

For Sugahara lab accelerator partner work:
1. Measure: rotating coil → integrated C_n
2. Simulate: Radia 3D model with proposed end-pole ∆ parameter
3. Compare: integrated C_n_sim should match C_n_measured
4. Iterate: adjust ∆ until match
5. Local field: use the validated 3D model to predict B(x, y, z)
   along beam orbit
6. Cross-validate: use radia-mcp.mathematica for symbolic
   verification of any analytical end-pole formulae

This is the canonical magnet validation cycle, supported entirely
by the radia-mcp tool stack.
"""


def get_accelerator_documentation(topic: str = "all") -> str:
    """Dispatch by topic.

    Topics:
      "all"
      "end_pole"         - Analytical chamfer design (Delferriere)
      "kolkata"          - Radia + TOSCA validation case study
      "rotating_coil"    - Multipole measurement + field reconstruction
    """
    topic = topic.lower().strip()
    if topic in ("end_pole", "chamfer", "delferriere"):
        return END_POLE_DESIGN
    if topic in ("kolkata", "cyclotron", "validation"):
        return KOLKATA_CYCLOTRON
    if topic in ("rotating_coil", "multipole", "measurement"):
        return ROTATING_COIL_MEASUREMENT
    if topic == "all":
        return "\n\n".join([
            END_POLE_DESIGN, KOLKATA_CYCLOTRON, ROTATING_COIL_MEASUREMENT,
        ])
    return (
        f"Unknown topic '{topic}'. Available: "
        "all, end_pole, kolkata, rotating_coil."
    )
