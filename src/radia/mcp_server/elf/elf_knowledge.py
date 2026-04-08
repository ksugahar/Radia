"""
ELF600 electromagnetic field analysis knowledge base.

Covers file formats (.mai, .mei, .meg), three solvers (MAGIC, ELFIN, BEAM),
element types, mesh scripting, B-H curves, and analysis workflows.
"""

# ============================================================
# Overview
# ============================================================

ELF_OVERVIEW = """\
# ELF600 Electromagnetic Field Analysis Suite

ELF600 is a **BEM (Boundary Element Method)** based electromagnetic analysis
suite consisting of three solver modules:

| Module | Purpose | Key Application |
|--------|---------|-----------------|
| **MAGIC** | Magnetostatic field analysis | Transformers, actuators, motors, magnets |
| **ELFIN** | Eddy current / low-frequency EM | Induction devices, shielding, NDE |
| **BEAM** | Charged particle beam tracking | Electron optics, accelerators |

## File Types

| Extension | Description | Role |
|-----------|-------------|------|
| `.mei` | Mesh script | Procedural mesh generation language |
| `.meg` | Compiled mesh data | Node/element data (output of IEmesh) |
| `.mai` | Analysis input | Solver control, material data, SOL blocks |
| `.mag` | Result data | Solver output (magnetization, fields) |
| `.model` | Motor config | Key-value config for IPM motor analysis |
| `.props` | Motor properties | Additional motor parameters |

## Workflow
```
.mei (mesh script)  -->  IEmesh processor  -->  .meg (compiled mesh)
.mai (analysis) + .meg (mesh)  -->  MAGIC/ELFIN/BEAM solver  -->  results
```

## Symmetry Types
- **T** = Full 3D
- **K** = 2D axisymmetric (revolution about Z axis)
- **R** = 2D cross-section (infinitely long in Z)
"""

# ============================================================
# .mai format
# ============================================================

MAI_FORMAT = """\
# .mai File Format (Analysis Input)

The `.mai` file controls solver execution: material properties, excitation,
solver commands, and output options.

## Structure

```
USE  <SOLVER>  3.50        # Header: MAGIC, ELFIN, or BEAM
PRE                        # Preprocessing block
  <material/excitation commands>
END
SOL <TYPE>                 # Solver execution block (repeatable)
  <solver options>
END
SOL <TYPE>                 # Additional solver blocks
  ...
END
```

## PRE Block Keywords

### MAGIC PRE Keywords
| Keyword | Description | Example |
|---------|-------------|---------|
| `HBUN` | B-H curve unit (OE/G or A/M/T) | `HBUN 1 OE/G` |
| `HBSC` | B-H curve scale factors | `HBSC 1 1.0 1.0` |
| `HBCU` | B-H curve data point (H, B) | `HBCU 1 0.0 0.0` |
| `HBRM` | Remanent magnetization (permanent magnet) | `HBRM 2 1 1.35 0.0 0.0 1.0` |
| `HBCN` | B-H curve assignment to material | `HBCN 1 1` |
| `COI1` | Coil definition (turns, layers) | `COI1 1 0 100 1` |
| `AMP1` | Coil current amplitude (real) | `AMP1 1 0 1.0` |
| `AMP1I` | Coil current amplitude (imaginary) | `AMP1I 1 0 0.5` |
| `CMU1` / `CMU1I` | Complex permeability (real/imag) | `CMU1 1 0 1000.0` |
| `OHM2` / `OHM3` | Resistance values | `OHM3 1 0 0.01` |
| `ORI1` | Rotation axis origin (for IPM motors) | `ORI1 0.0 0.0 0.0 0.0 0.0 1.0` |
| `MOV1` | Movement/rotation per timestep | `MOV1 1 7.5` |
| `UNI2` | Uniform external field | `UNI2 1 0 0.0 0.0 1.0` |

### ELFIN PRE Keywords
| Keyword | Description | Example |
|---------|-------------|---------|
| `VOL1` / `VOL1I` | Voltage source (real/imag) | `VOL1 1 0 1.0` |
| `EDSC` | Conductivity scale | `EDSC 1 1.0` |
| `EDCU` | Conductivity curve data | `EDCU 1 5.8e7` |
| `CHA2` | Charge definition | `CHA2 1 0 1.0` |

### BEAM PRE Keywords
| Keyword | Description | Example |
|---------|-------------|---------|
| `CHAR` | Particle charge (electron) | `CHAR 1.602E-19` |
| `MASS` | Particle mass | `MASS 9.109E-31` |
| `VOLB` | Beam voltage | `VOLB 1000.0` |
| `FILE` | Import field from MAGIC/ELFIN | `FILE MAGIC BX01.meg` |
"""

# ============================================================
# SOL commands
# ============================================================

SOL_COMMANDS = """\
# SOL Block Types and Sub-Keywords

SOL blocks define solver execution. Multiple SOL blocks can appear in one .mai.

## SOL Types

### MAGIC SOL Types
| SOL Type | Description |
|----------|-------------|
| `SOL MOME` | Moment method (BEM) — main magnetostatic solver |
| `SOL MOMC` | Moment method complex — AC/frequency domain |
| `SOL FIEL` | Field computation at evaluation points |
| `SOL FIXA` | Flux linkage computation |
| `SOL FIXB` | Fixed B-field computation |
| `SOL FORC` | Force computation (Maxwell stress tensor) |
| `SOL FORT` | Torque computation |
| `SOL MAS` | Mass/magnetization output |

### ELFIN SOL Types
| SOL Type | Description |
|----------|-------------|
| `SOL MOME` | Moment method (BEM) — main eddy current solver |
| `SOL FIEL` | Field computation at evaluation points |

### BEAM SOL Types
| SOL Type | Description |
|----------|-------------|
| `SOL BEAM` | Charged particle beam tracking |

## Common SOL Sub-Keywords

| Keyword | Description | Example |
|---------|-------------|---------|
| `NCPU 0` | Auto CPU count (0 = all available) | `NCPU 0` |
| `BLAS 8` | BLAS block size | `BLAS 8` |
| `TIME` | Time stepping (step, dt, tmax) | `TIME 1 0.001 0.1` |
| `NONL` | Nonlinear iteration (maxiter, tol) | `NONL 50 0.001` |
| `DMEG` | Output results to .meg file | `DMEG` |
| `PASS GENE` | Skip (pass) this SOL block | `PASS GENE` |
| `PASS STEP ALL` | Pass all time steps | `PASS STEP ALL` |
| `NOGO` | Dry run (no execution) | `NOGO` |
| `FREQ` | Frequency for AC analysis (MOMC) | `FREQ 50.0` |
| `EMFM` | EMF computation for coils | `EMFM 1` |
| `FLUM` | Flux linkage for coil | `FLUM 1` |
| `COIM` | Coil impedance | `COIM 1` |
| `SELM` | Select elements for force integration | `SELM 1 3` |
| `DELM` | Deselect elements | `DELM 1 5` |
| `STAR` | Star connection (motor coils) | `STAR` |
| `STEP` | Beam tracking step count | `STEP 1000` |
| `STOP` | Beam tracking bounds | `STOP 0.0 0.1` |
| `RELA` | Relativistic correction (beam) | `RELA` |
| `BMEG` | Beam trajectory output to .meg | `BMEG` |
"""

# ============================================================
# .mei format
# ============================================================

MEI_FORMAT = """\
# .mei File Format (Mesh Script)

The `.mei` file is a **procedural mesh scripting language** for building
BEM surface/line models. It is NOT raw mesh data — it is processed by
IEmesh to generate .meg files.

## Structure

```
MODEL <SOLVER> <SYMMETRY>    # Header
GSC <scale>                  # Global scale (e.g., 0.001 for mm->m)
IMA X                        # Image symmetry plane(s)
PUT A=10.0                   # Variable assignment

<source region commands>     # Coils, magnets, conductors

SPACE                        # Separator

<evaluation region commands> # Field evaluation surfaces

STOP                         # End of script
```

## Header
- `MODEL MAGIC T` — 3D magnetostatic
- `MODEL MAGIC K` — axisymmetric magnetostatic
- `MODEL ELFIN R` — 2D cross-section eddy current
- `MODEL BEAM T` — 3D beam tracking

## Variable Assignment
```
PUT A=10.0
PUT B=A*2.0
PUT C=_COS(30.0)*A         # Built-in functions: _COS, _SIN, _TAN, _SQRT, _ABS
```

## Two Regions
1. **Before SPACE**: Source geometry (coils, iron, magnets, conductors)
2. **After SPACE**: Evaluation geometry (field computation surfaces/points)
"""

# ============================================================
# .mei commands
# ============================================================

MEI_COMMANDS = """\
# Mesh Scripting Commands (.mei)

## Auto-Mesh Commands (AA)
Generate mesh elements automatically for standard shapes.

| Command | Shape | Symmetry |
|---------|-------|----------|
| `AA(RECT,x1,y1,x2,y2,nx,ny)` | Rectangle | K, R |
| `AA(BLOCK,x1,y1,z1,x2,y2,z2,nx,ny,nz)` | 3D block | T |
| `AA(ARCH2,r1,r2,a1,a2,nr,na)` | Arch/annulus (2D) | K, R |
| `AA(ARCH3,r1,r2,a1,a2,h1,h2,nr,na,nh)` | Torus sector (3D) | T |
| `AA(RING1,r1,r2,a1,a2,nr,na)` | Ring surface | T |
| `AA(RINGB,r1,r2,a1,a2,h,nr,na,nh)` | Ring block | T |
| `AA(SQRING,...)` | Square cross-section ring | T |
| `AA(BOX,...)` | Hollow box | T |
| `AA(LINEX,x1,y1,x2,y2,n)` | Line elements | K, R |
| `AA(O,x,y,z)` | Set origin for next AA | all |
| `AA(P,XZ)` / `AA(P,YZ)` | Set plane orientation | all |

## Node Commands (AN)
| Command | Description |
|---------|-------------|
| `AN(XYZ,x,y,z)` | Add node at (x,y,z) |
| `AN(RT2,r,theta)` | Add node in cylindrical coords |
| `AN(FL,n1,n2,count)` | Interpolate nodes between n1 and n2 |
| `AN(F3,n1,n2)` | Add midpoint between n1 and n2 |

## Element Commands (AE)
| Command | Description |
|---------|-------------|
| `AE(N2,n1,n2)` | Add 2-node line element |
| `AE(N4,n1,n2,n3,n4)` | Add 4-node quad element |
| `AF(...)` | Add face |

## Material/Name (NB)
| Command | Description |
|---------|-------------|
| `NB(MID,n)` | Set material ID for subsequent elements |
| `NB(NAME,xxx)` | Set element label (e.g., MWL, MCL, MMB) |
| `NB(CHAN,MM,old,new)` | Change material ID |
| `NB(CHAN,MN,id,name)` | Change element name |

## Selection Commands
| Command | Description |
|---------|-------------|
| `SE(ADD,ALL)` | Select all elements |
| `SE(ADD,n1,n2)` | Select element range |
| `SE(AND,...)` | Intersect selection |
| `SN(ADD,...)` | Select nodes |

## Mesh Operations (ME)
| Command | Description |
|---------|-------------|
| `ME(CPL,dx,dy,dz,n)` | Copy elements (linear translation) |
| `ME(CPM,plane)` | Copy elements (mirror) |
| `ME(EXL,dx,dy,dz,n)` | Extrude elements (linear) |
| `ME(EXC,ax,ay,az,angle,n)` | Extrude elements (circular/revolve) |
| `ME(MVL,dx,dy,dz)` | Move elements (linear) |
| `ME(MVC,ax,ay,az,angle)` | Move elements (circular) |
| `ME(DEL)` | Delete selected elements |

## Other Commands
| Command | Description |
|---------|-------------|
| `MN(TAX,...)` | Transform node coordinates (axis swap) |
| `BE(L4,...)` | Bisect elements |
| `TE(TRI,...)` | Convert quads to triangles |
| `TE(REV)` | Reverse element normals |
| `RB(MERGE,tol)` | Merge coincident nodes |
| `CL(LONE)` | Clean up orphan nodes |
| `VEC(X,dx,dy,dz)` | Define vector for copy/move |
| `R/<count>` | Repeat next command count times |
| `LOOP` / `LOOPEND` | Loop block |
| `SUBR` / `CALL` / `RETURN` | Subroutines |
"""

# ============================================================
# .meg format
# ============================================================

MEG_FORMAT = """\
# .meg File Format (Compiled Mesh Data)

The `.meg` file is fixed-format text output from IEmesh, containing node
coordinates and element connectivity.

## Structure

```
BOOK  MEP  3.50                    # Header
* ELF/MESH VERSION 5.3.7          # Version comment
* SOLVER = MAGIC                   # Solver type
MGSC 0.001                        # Global scale (M=MAGIC, E=ELFIN, B=BEAM)
MIMA -X                           # Image symmetry

<source nodes and elements>        # GR1 group
<evaluation nodes and elements>    # GR2 group

BOOK  END                          # Footer
```

## Node Records
```
MGR1  <id>  <timestep>  <x>  <y>  <z>      # Source node (MAGIC)
MGR2  <id>  <timestep>  <x>  <y>  <z>      # Evaluation node (MAGIC)
EGR1  <id>  <timestep>  <x>  <y>  <z>      # Source node (ELFIN)
BGR1  <id>  <timestep>  <x>  <y>  <z>      # Source node (BEAM)
```

## Element Records
```
<TYPE>  <id>  <timestep>  <material_id>  <n1>  <n2>  [<n3>  <n4>]
```

Element TYPE follows naming convention: **PREFIX(3) + NODECOUNT(1) + DIM(1)**

Example: `MWL4K  1  0  1  1  2  3  4` = MAGIC Wall 4-node Axisymmetric

## Prefix Codes in Global Settings
- `MGSC` / `EGSC` / `BGSC` — scale factor (MAGIC/ELFIN/BEAM)
- `MIMA` / `EIMA` / `BIMA` — image symmetry
"""

# ============================================================
# Element types
# ============================================================

ELEMENT_TYPES = """\
# ELF Element Type Naming Convention

Element names follow: **PREFIX(3 chars) + NODE_COUNT(1 digit) + DIM(1 char)**

## DIM suffix
| Suffix | Symmetry | Description |
|--------|----------|-------------|
| `T` | 3D | Full three-dimensional |
| `K` | Axisymmetric | Revolution about Z axis |
| `R` | 2D cross-section | Infinitely long in Z |

## MAGIC Element Prefixes (Magnetostatic)
| Prefix | Full Name | Description |
|--------|-----------|-------------|
| `MCL` | MAGIC Coil | Current-carrying coil |
| `MMB` | MAGIC Magnetic Body | Permeable body (linear) |
| `MWL` | MAGIC Wall | Non-linear iron surface |
| `MBB` | MAGIC Bulk Body | Volume element |
| `MCO` | MAGIC Contour | Field evaluation surface |
| `MAB` | MAGIC Air Body | Air region body |
| `MAT` | MAGIC Air Triangle | Air region triangle |
| `MMS` | MAGIC Magnetic Surface | 3D-only magnetic surface |
| `MMP` | MAGIC Magnetic Permanent | Permanent magnet surface |
| `MWV` | MAGIC Wall V | Wall variant |
| `MCM` | MAGIC Complex Material | Complex permeability surface |

## ELFIN Element Prefixes (Eddy Current)
| Prefix | Full Name | Description |
|--------|-----------|-------------|
| `ESC` | ELFIN Source Conductor | Current source conductor |
| `EMB` | ELFIN Magnetic Body | Magnetic/conductive body |
| `ECB` | ELFIN Conductive Body | Conductive body |
| `ECO` | ELFIN Contour | Field evaluation surface |

## BEAM Element Prefixes
| Prefix | Full Name | Description |
|--------|-----------|-------------|
| `BGR` | BEAM Grid | Grid nodes (source/evaluation) |

## Examples
| Element Code | Meaning |
|-------------|---------|
| `MCL2K` | MAGIC Coil, 2-node line, axisymmetric |
| `MWL4K` | MAGIC Wall, 4-node quad, axisymmetric |
| `MMB4T` | MAGIC Magnetic Body, 4-node quad, 3D |
| `MCO4K` | MAGIC Contour, 4-node quad, axisymmetric |
| `ESC2K` | ELFIN Source Conductor, 2-node line, axisymmetric |
| `EMB4R` | ELFIN Magnetic Body, 4-node quad, 2D cross-section |

## Available Labels by DIM (for MEG export)

### T (3D)
MMB, MMS, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB

### K (Axisymmetric) and R (2D)
MMB, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB (same minus MMS)
"""

# ============================================================
# B-H Curves
# ============================================================

BH_CURVES = """\
# B-H Curve Specification in ELF

B-H curves define the nonlinear magnetic permeability of iron/steel materials.

## Unit Declaration
```
HBUN  <curve_id>  <unit_system>
```
- `OE/G` — Oersted / Gauss (CGS)
- `A/M/T` — Ampere/meter / Tesla (SI)

## Scale Factors
```
HBSC  <curve_id>  <H_scale>  <B_scale>
```

## Curve Data Points
```
HBCU  <curve_id>  <H_value>  <B_value>
```
Multiple HBCU lines define the piecewise-linear B-H curve.

## Curve-to-Material Assignment
```
HBCN  <material_id>  <curve_id>
```

## Permanent Magnet (Remanence)
```
HBRM  <material_id>  <curve_id>  <Br>  <dir_x>  <dir_y>  <dir_z>
```
- `<Br>` = remanent flux density
- `<dir_x/y/z>` = magnetization direction unit vector

## Example (B-H curve for silicon steel in SI)
```
PRE
  HBUN  1  A/M/T
  HBSC  1  1.0  1.0
  HBCU  1      0.0    0.0
  HBCU  1    200.0    1.0
  HBCU  1    500.0    1.4
  HBCU  1   1000.0    1.6
  HBCU  1   5000.0    1.8
  HBCU  1  50000.0    2.0
  HBCN  1  1
END
```

## Example (NdFeB permanent magnet, Br = 1.35 T, Z-direction)
```
PRE
  HBRM  2  1  1.35  0.0  0.0  1.0
END
```
"""

# ============================================================
# MAGIC solver
# ============================================================

MAGIC_DOCS = """\
# MAGIC Solver (Magnetostatic BEM)

MAGIC solves magnetostatic field problems using the Boundary Element Method.

## Applications
- Transformers and reactors
- Electromagnets and solenoids
- Permanent magnet devices
- Electric motors (IPM, SPM)
- Magnetic actuators and relays

## Typical .mai Structure
```
USE  MAGIC  3.50
PRE
  HBUN  1  A/M/T             # B-H curve unit
  HBCU  1  0.0  0.0          # B-H data points
  HBCU  1  200.0  1.0
  ...
  HBCN  1  1                  # Assign curve to material
  COI1  1  0  100  1          # Coil: 100 turns
  AMP1  1  0  5.0             # Current: 5A
END
SOL MOME                      # Main BEM solve
  NCPU 0                      # Auto CPU
  BLAS 8                      # BLAS block size
  NONL 50 0.001               # Nonlinear: 50 iter, 0.1% tol
  DMEG                        # Write results to .meg
END
SOL FIEL                      # Field computation
  DMEG
END
SOL FORC                      # Force computation
  SELM 1 3                    # Elements for Maxwell stress
  DMEG
END
```

## Solver Chain
1. `SOL MOME` — assemble and solve BEM matrix (mandatory first step)
2. `SOL FIEL` — compute B, H at evaluation points (after MOME)
3. `SOL FORC` — force on selected surfaces (after MOME)
4. `SOL FORT` — torque on selected surfaces (after MOME)
5. `SOL FIXA` — flux linkage through coils (after MOME)
6. `SOL MAS` — magnetization output (after MOME)

## Time-Stepping (Transient/Motor)
```
SOL MOME
  TIME 1 0.001 0.1           # Start step, dt, end time
  NONL 50 0.001
  DMEG
END
```

## Example Categories
- **MK/** — Axisymmetric (K-symmetry): solenoids, pot magnets
- **MR/** — 2D cross-section (R-symmetry): C-core, relay
- **MT/** — Full 3D: complex geometry, 3D magnets
- **IPM/** — Interior permanent magnet motors
- **LscLl/** — Inductance computation
- **MOMC/** — AC complex analysis
- **magne/** — Magnetization/demagnetization
"""

# ============================================================
# ELFIN solver
# ============================================================

ELFIN_DOCS = """\
# ELFIN Solver (Eddy Current BEM)

ELFIN solves low-frequency electromagnetic (eddy current) problems using BEM.

## Applications
- Eddy current NDE (non-destructive evaluation)
- Induction heating (coil + workpiece)
- Electromagnetic shielding
- Transient electromagnetic analysis

## Typical .mai Structure
```
USE  ELFIN  3.50
PRE
  VOL1  1  0  1.0            # Voltage source: 1V
  EDSC  1  1.0               # Conductivity scale
  EDCU  1  5.8e7             # Conductivity (copper: 5.8e7 S/m)
END
SOL MOME
  NCPU 0
  BLAS 8
  TIME 1 0.0001 0.01         # Time stepping
  DMEG
END
SOL FIEL
  DMEG
END
```

## Element Types
- **ESC** (Source Conductor) — carries imposed current/voltage
- **EMB** (Magnetic Body) — permeable + conductive material
- **ECB** (Conductive Body) — non-magnetic conductor
- **ECO** (Contour) — evaluation surface (after SPACE in .mei)

## Example Categories
- **EK/** — Axisymmetric eddy current (induction heating, NDE probes)
- **ER/** — 2D cross-section eddy current
- **ET/** — Full 3D eddy current
"""

# ============================================================
# BEAM solver
# ============================================================

BEAM_DOCS = """\
# BEAM Solver (Charged Particle Tracking)

BEAM tracks charged particles through electromagnetic fields computed by
MAGIC or ELFIN.

## Typical .mai Structure
```
USE  BEAM  3.50
PRE
  CHAR  1.602E-19             # Electron charge
  MASS  9.109E-31             # Electron mass
  VOLB  10000.0               # Beam voltage: 10 kV
  FILE  MAGIC  BX01.meg       # Import field from MAGIC result
END
SOL BEAM
  STEP 1000                   # Number of tracking steps
  STOP 0.0 0.1               # Tracking bounds (z_min, z_max)
  RELA                        # Relativistic correction
  BMEG                        # Output trajectory to .meg
END
```

## Features
- Import EM fields from MAGIC or ELFIN via `FILE` command
- Relativistic correction with `RELA`
- Space charge effects
- Multiple beam sources

## Example Categories
- **BX01-BX08** — Various electron optics configurations
- **BECB** — Beam in ELFIN-computed eddy current fields
"""

# ============================================================
# IPM Motor
# ============================================================

IPM_MOTOR = """\
# IPM Motor Analysis with MAGIC

Interior Permanent Magnet (IPM) motor analysis uses MAGIC with time-stepping
rotation.

## Required Files
1. `.mei` — Motor mesh (rotor + stator surfaces)
2. `.meg` — Compiled mesh
3. `.mai` — Analysis input with rotation
4. `.model` — Motor configuration
5. `.props` — Motor properties

## .model File Format
```
PartsM1 {
  RotAxis    = 0.0, 0.0, 0.0, 3     # Origin + axis (3=Z)
  Params_Rot = 1, 0.0, 60.0, 0.0    # Rotation parameters
}
Calculations {
  TimeSeries = 180, 0.5, 180.0      # Steps, dt (deg), total angle
}
```

## .mai for Motor Analysis
```
USE  MAGIC  3.50
PRE
  HBUN  1  A/M/T
  HBCU  1  ...                       # Rotor/stator iron B-H curves
  HBRM  2  1  1.35  0.0  0.0  1.0   # NdFeB magnets
  COI1  1  0  10  1                  # Phase A coil
  COI1  2  0  10  1                  # Phase B coil
  COI1  3  0  10  1                  # Phase C coil
  AMP1  1  0  5.0                    # Phase currents
  AMP1  2  0  5.0
  AMP1  3  0  5.0
  ORI1  0.0  0.0  0.0  0.0  0.0  1.0  # Rotation axis (Z)
  MOV1  1  7.5                        # 7.5 deg/step
END
SOL MOME
  TIME 1 1.0 48.0                    # 48 time steps
  NONL 50 0.001
  DMEG
END
SOL FORT                             # Torque
  SELM 1 3
  DMEG
END
SOL FIXA                             # Flux linkage (for back-EMF)
  FLUM 1
  FLUM 2
  FLUM 3
  STAR                               # Star connection
  DMEG
END
```

## Example Files
- `magic/IPM/Motor1.mai` / `Motor1.meg` / `Motor1.model`
- `magic/IPM/Motor2.mai` / `Motor2.meg` / `Motor2.model`
"""

# ============================================================
# Inductance computation
# ============================================================

INDUCTANCE_DOCS = """\
# Inductance Computation (LscLl)

ELF can compute short-circuit inductance (Lsc) and leakage inductance (Ll)
for transformers using the MAGIC solver.

## Workflow
1. Model transformer windings as MCL (coil) elements
2. Model core as MWL (wall) elements with B-H curve
3. Run SOL MOME for field solution
4. Run SOL FIXA with FLUM for flux linkage
5. Compute Lsc/Ll from flux linkage and current

## Example Files
Located in `magic/LscLl/`:
- `Sample1` through `Sample6` — various transformer configurations
- Each includes `.mai`, `.meg`, `.mei`, and `_DATA処理.xls` for post-processing

## Key .mai Pattern for Inductance
```
SOL FIXA
  FLUM 1                     # Flux linkage for coil 1
  FLUM 2                     # Flux linkage for coil 2
  EMFM 1                     # EMF for coil 1
  COIM 1                     # Impedance for coil 1
  DMEG
END
```
"""

# ============================================================
# Magnetization / demagnetization
# ============================================================

MAGNETIZATION_DOCS = """\
# Magnetization / Demagnetization Analysis

ELF MAGIC can simulate magnetization and demagnetization processes for
permanent magnets.

## Directory Structure
- `magic/magne/magnetize/` — Magnetization analysis (5 cases + sample)
- `magic/magne/demagne/` — Demagnetization analysis (4 samples)

## Magnetization Workflow (CASE files)
Each case has:
- `CASE*-A.mei` — Analysis mesh definition
- `CASE*-C.mai` + `CASE*-C.meg` — Computation input
- `CASE*-PARA.txt` — Step parameters (multi-step for CASE4/5)
- `ReadMe*.txt` — Case description
- `result/` — Pre-computed results

## Multi-Step Magnetization (CASE4, CASE5)
Cases 4 and 5 use multiple parameter files for multi-step magnetization:
```
CASE4-PARA-1.txt    # Step 1 parameters
CASE4-PARA-2.txt    # Step 2 parameters
CASE4-PARA-3.txt    # Step 3 parameters
```

## Demagnetization Samples
- `Sample01` through `Sample04` — various demagnetization scenarios
- See `demag.pdf` for documentation
"""

# ============================================================
# MEG Export
# ============================================================

MEG_EXPORT = """\
# MEG Export from Cubit

The Cubit mesh export tool (`radia_export meg`) generates ELF-compatible .meg files.

## Command
```
radia_export meg "output.meg" [threed|twod|axisymmetric] [labels "1:MMB,2:MWL,..."] [overwrite]
```

## DIM Selection
| Keyword | DIM suffix | Symmetry |
|---------|-----------|----------|
| `threed` | T | Full 3D |
| `twod` | R | 2D cross-section |
| `axisymmetric` | K | Axisymmetric |

## Label Assignment
Each Cubit block is assigned an ELF element label:
```
labels "1:MMB,2:MWL,3:MCL,4:MCO"
```

## Element Name = PREFIX(3) + NODE_COUNT(1) + DIM(1)
- The node count is determined from the mesh element type (quad=4, tri=3, line=2)
- The DIM suffix is determined from the symmetry selection

## Available Labels by DIM

### T (3D)
MMB, MMS, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB

### K (Axisymmetric) / R (2D)
MMB, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB

## Example
```python
# In Cubit Python console:
cubit.cmd('radia_export meg "model.meg" axisymmetric labels "1:MWL,2:MCL,3:MCO" overwrite')
```
"""

# ============================================================
# Examples catalog
# ============================================================

EXAMPLES_CATALOG = """\
# ELF600 Example File Catalog

Located in `C:\\ELF600\\examples\\`.

## beam/ — Charged Particle Beam Tracking
| Example | Description |
|---------|-------------|
| `BASIC/ECB` | Beam in ELFIN eddy current field |
| `BASIC/BECB` | Beam + eddy current combined |
| `BX01` - `BX08` | Various electron optics (standard + variants) |
| `BX01m` - `BX08m` | Modified versions |
| `BX05e` - `BX08e` | Extended versions |

## elfin/ — Eddy Current Analysis
| Example | Description |
|---------|-------------|
| `BASIC/EMB` | Basic magnetic body with eddy currents |
| `BASIC/EQL` | Q-L formulation |
| `BASIC/ESC*` | Source conductor variants (F, Q, V) |
| `EK/EK01` - `EK06` | Axisymmetric eddy current problems |
| `ER/ER01` - `ER05` | 2D cross-section eddy current |
| `ET/ET01` - `ET06` | 3D eddy current problems |
| `WorkBook/EWork*` | Workbook tutorials |

## magic/ — Magnetostatic Analysis
| Example | Description |
|---------|-------------|
| `BASIC/MBCL` | Magnetic body + coil (basic) |
| `BASIC/MBCLV` | With volume elements |
| `BASIC/MBCLXYZ` | With XYZ output |
| `BASIC/HCOIL` | Helmholtz coil |
| `BASIC/POLY` | Polygonal geometry |
| `BASIC/WVF/WVH/WVU` | Wave function examples |
| `BASIC/FLUMD/FLUMU` | Flux density examples |
| `BASIC/EMFU` | EMF computation |
| `MK/MK01` - `MK06` | Axisymmetric magnetostatic |
| `MR/MR01` - `MR05` | 2D cross-section magnetostatic |
| `MT/MT01` - `MT09` | 3D magnetostatic |
| `IPM/Motor1,Motor2` | IPM motor analysis |
| `LscLl/Sample1-6` | Inductance computation |
| `MOMC/MOMCFJ` | AC complex analysis |
| `magne/magnetize/` | 5 magnetization cases |
| `magne/demagne/` | 4 demagnetization samples |
| `V6Conv/` | Version 6 conversion tool |
| `WorkBook/Work*` | Workbook tutorials |

## File Naming Convention
- Base name without suffix → standard version
- `m` suffix → modified version
- `e` suffix → extended version
- `_01`, `_02` suffixes → sub-variants (e.g., different symmetry)
"""

# ============================================================
# Router function
# ============================================================

_TOPICS = {
    "overview": ELF_OVERVIEW,
    "mai_format": MAI_FORMAT,
    "mei_format": MEI_FORMAT,
    "meg_format": MEG_FORMAT,
    "magic": MAGIC_DOCS,
    "elfin": ELFIN_DOCS,
    "beam": BEAM_DOCS,
    "element_types": ELEMENT_TYPES,
    "bh_curves": BH_CURVES,
    "sol_commands": SOL_COMMANDS,
    "mei_commands": MEI_COMMANDS,
    "ipm_motor": IPM_MOTOR,
    "inductance": INDUCTANCE_DOCS,
    "magnetization": MAGNETIZATION_DOCS,
    "examples": EXAMPLES_CATALOG,
    "meg_export": MEG_EXPORT,
}


def get_elf_documentation(topic: str = "all") -> str:
    """Return ELF documentation for the given topic."""
    if topic == "all":
        return "\n\n".join(_TOPICS.values())

    result = _TOPICS.get(topic)
    if result is None:
        available = ", ".join(sorted(_TOPICS.keys()))
        return f"Unknown topic: '{topic}'. Available topics: {available}"
    return result
