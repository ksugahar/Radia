"""
ELF600 electromagnetic field analysis knowledge base.

Covers file formats (.mai, .mei, .meg), three solvers (MAGIC, ELFIN, BEAM),
element types, mesh scripting, B-H curves, analysis workflows, tools,
error messages, and advanced topics (anisotropy, lamination, demagnetization,
sinusoidal response, IPM motors, inductance, magnetization).

Organized into 30+ topics derived from:
- C:/ELF600/examples/ (all beam/elfin/magic examples)
- C:/ELF600/help/m_rf1/ (MAGIC reference manual)
- C:/ELF600/help/e_rf1/ (ELFIN reference manual)
- C:/ELF600/help/b_rf1/ (BEAM reference manual)
- C:/ELF600/help/m_treasure/ (quick-reference tables)
- C:/ELF600/help/d_ken/ (training materials and technical docs)
- C:/ELF600/help/u_support/ (error messages)
- C:/ELF600/help/t_mesh/ (IEmesh mesh commands)
- C:/ELF600/help/t_iemesh/ (IEmesh tool overview)
"""

# ============================================================
# Overview
# ============================================================

ELF_OVERVIEW = """\
# ELF600 Electromagnetic Field Analysis Suite

ELF600 is a **BEM (Boundary Element Method / Integral Element Method)**
electromagnetic analysis suite with three solver modules:

| Module | Purpose | Key Application |
|--------|---------|-----------------|
| **MAGIC** | Magnetostatic field (static, transient, AC) | Transformers, motors, magnets |
| **ELFIN** | Electrostatic field analysis | Capacitors, insulators, electrodes |
| **BEAM** | Charged particle beam tracking | Electron optics, accelerators |

**Note:** Despite the name "ELFIN", the ELFIN solver handles **electrostatic**
(not eddy current) analysis. Eddy current analysis in ELF is done by MAGIC
using MAB/MAT/MBB elements with time-stepping.

## File Types

| Extension | Description | Role |
|-----------|-------------|------|
| `.mei` | Mesh script | Procedural mesh generation language |
| `.meg` | Compiled mesh data | Node/element data (output of IEmesh) |
| `.mai` | Analysis input | Solver control, material data, SOL blocks |
| `.mao` | Computation log | Convergence info, errors, timing |
| `.mag` | Result data | Fields, forces, flux linkage |
| `.mas` | Charge/source file | For BEAM trajectory input |
| `.mat` | Matrix file | LU decomposition cache |
| `.mac` | Work file | Intermediate timestep data |
| `.model` | Motor config | Key-value config for IPM analysis |

## Workflow
```
.mei (mesh script) --> IEmesh (mesh750.exe) --> .meg (compiled mesh)
.mai (analysis) + .meg (mesh) --> MAGIC/ELFIN/BEAM --> .mag/.mao results
```

## Symmetry Types
- **T** = Full 3D
- **K** = 2D (infinite in Z, cross-section in XY plane)
- **R** = Axisymmetric (revolution about Z axis, cross-section in XZ plane)

## Tools
- **IEmesh** = Mesh script editor + mesh generator (mesh750.exe) + viewer (Wmap3)
- **MaiEditor3** = GUI editor for .mai control files
- **Wmap3** = 3D visualization of .meg and .mag files
- **MagFilter2** = Result extraction to CSV
- **ELF/Bench** = Problem/project manager with batch processing
- **ComplexTransfer2** = Convert MOMC complex results to pseudo-transient for animation
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
PRE  <MAXMID> <MAXSTEP>   # Preprocessing block
  <material/excitation commands>
END
SOL <TYPE>                 # Solver execution block (repeatable)
  <solver options>
END
```

- Comment lines: first character = `*` (leading spaces/tabs OK)
- `*` before `SOL` skips that entire block
- **Max line length: 120 half-width characters**

## MAGIC PRE Keywords

### B-H Curves
| Keyword | Description | Example |
|---------|-------------|---------|
| `HBUN` | B-H curve unit (OE/G or A/M/T) | `HBUN 1 A/M T` |
| `HBSC` | B-H curve scale factors | `HBSC 1 1.0 1.0` |
| `HBCU` | B-H curve data point (H, B) | `HBCU 1 0.0 0.0` |
| `HBCN` | B-H curve assignment per step (demagnetization) | `HBCN 1 0 1` |
| `HBRM` | Recoil permeability + max remanence (demagnetization) | `HBRM 1 1.05 1.35` |
| `HBEQ` | Copy B-H curve (dest, source) | `HBEQ 2 1` |
| `HBCP` | Copy B-H curve (source, dest1, dest2, ...) | `HBCP 1 2 3 4` |

### Coil / Current
| Keyword | Description | Example |
|---------|-------------|---------|
| `COI1` | Coil definition (MID, 0, turns, NB1, NB2, A/V) | `COI1 1 0 100 8 4 A` |
| `AMP1` / `AMP1I` | Current amplitude real/imag | `AMP1 1 0 5.0` |
| `VOL1` / `VOL1I` | Voltage source real/imag | `VOL1 1 0 1.0` |
| `OHM1` | Coil resistance (normal calc) | `OHM1 1 0 0.1` |
| `OHM2` | Volume resistivity (eddy current elements) | `OHM2 3 0 2.7e-8` |
| `OHM3` / `OHM3I` | Impedance real/imag (sinusoidal calc) | `OHM3 1 0 0.1` |

### Anisotropy / Lamination
| Keyword | Description | Example |
|---------|-------------|---------|
| `HBA1` | Susceptibility anisotropy ratio (0 < k < 1) | `HBA1 1 0.5` |
| `HBA2` | Lamination stacking factor (0 < p < 1) | `HBA2 1 0.95` |
| `VEC1` | Direction vector per element (for HBA1/HBA2/MWV) | `VEC1 1 0 0.0 0.0 1.0` |
| `VEC3` | Direction vector per material | `VEC3 1 0 0.0 0.0 1.0` |

**HBA1 and HBA2 are mutually exclusive** -- cannot use both on same material.

### Thin Plate / Wire Diameter
| Keyword | Description | Example |
|---------|-------------|---------|
| `THIN` | Plate thickness [m] or wire diameter [m] | `THIN 1 0.001` |

Two uses: (1) MMT*/MAT* thin plate thickness, (2) MCL2T/MCL1* wire diameter.

### Motion
| Keyword | Description | Example |
|---------|-------------|---------|
| `MOV1` | Displacement + rotation per step | `MOV1 1 1 0.001 0 0 7.5` |
| `ORI1` | Rotation axis (origin + axis 1=X,2=Y,3=Z) | `ORI1 1 0 0 0 0 3` |
| `MVEQ` / `MVCP` | Copy motion conditions | `MVEQ 2 1` |
| `STED` | Steady-state eddy current motion | `STED 3 0 0 0 5 10` |

### Uniform Field
| Keyword | Description | Example |
|---------|-------------|---------|
| `UNI1` / `UNI1I` | Uniform field A/m (real/imag) | `UNI1 0 0 0.0 0.0 100.0` |
| `UNI2` / `UNI2I` | Uniform field Oe (real/imag) | `UNI2 0 0 0.0 0.0 1.0` |

### Complex Permeability (SOL MOMC only)
| Keyword | Description | Example |
|---------|-------------|---------|
| `CMU1` | Complex permeability real (>= 1) | `CMU1 1 0 1000.0` |
| `CMU1I` | Complex permeability imag (<= 0) | `CMU1I 1 0 -10.0` |

Replaces B-H curves for sinusoidal analysis. **Nonlinear B-H NOT supported in MOMC.**

### Output Control
| Keyword | Description | Example |
|---------|-------------|---------|
| `MAH` | Result file precision (+digit=space, -digit=comma) | `MAH 6` |
| `MAF` | Force output precision (FORT/FIXB per material) | `MAF 6` |

## ELFIN PRE Keywords
| Keyword | Description |
|---------|-------------|
| `EDSC` / `EDCU` | D-E curve scale and data points |
| `EDEQ` / `EDCP` | Copy D-E curve |
| `VOL1` | Electrode potential [V] |
| `CHA1` | Charge density [C/m^2] for EQL elements |
| `CHA2` | Total charge [C] on floating electrode (unknown potential) |
| `VAL1` | Boundary condition value (ESH: ratio, ESS/ESN: E-field) |
| `VEC2` / `VEC4` | Contour-point vector for cylindrical electrodes |

## BEAM PRE Keywords
| Keyword | Description | Example |
|---------|-------------|---------|
| `CHAR` | Particle charge [C] | `CHAR 1 1.602E-19` |
| `MASS` | Particle mass [kg] | `MASS 1 9.109E-31` |
| `VOLB` | Acceleration voltage [V] | `VOLB 1 10000.0` |
| `FILE` | Import field from MAGIC/ELFIN .mas | `FILE MAGIC problem_M` |
| `AMPB` | Beam current [A] for space-charge | `AMPB 1 0.001` |
"""

# ============================================================
# SOL commands
# ============================================================

SOL_COMMANDS = """\
# SOL Block Types and Sub-Keywords

SOL blocks define solver execution. Multiple SOL blocks can appear in one .mai.

## MAGIC SOL Types

| SOL Type | Description | After MOME | After MOMC |
|----------|-------------|------------|------------|
| `SOL MOME` | Moment method (BEM) -- main solver | -- | -- |
| `SOL MOMC` | Complex moment (sinusoidal/AC) | -- | -- |
| `SOL FIEL` | Field at evaluation points | Yes | Yes |
| `SOL FIXA` | Vector potential / flux linkage | Yes | Yes |
| `SOL FORC` | Element surface force | Yes | No |
| `SOL FORT` | Maxwell stress tensor force | Yes | No |
| `SOL FIXB` | Lorentz force on coils | Yes | No |
| `SOL MAS` | Charge file output (for BEAM) | Yes | No |
| `SOL FMAG` | Center field (MJH elements) | Yes | No |

**Key constraint:** SOL MOME or SOL MOMC must execute first. Choose one (not both).
After SOL MOMC, only SOL FIEL and SOL FIXA are available.

## ELFIN SOL Types
Same structure: SOL MOME (mandatory), then SOL FIEL, SOL FORC, SOL FORT, SOL MAS.

## BEAM SOL Type
Only `SOL BEAM` (trajectory computation).

## Common Sub-Keywords

### SOL MOME / SOL MOMC
| Keyword | Description | Example |
|---------|-------------|---------|
| `TIME` | Time stepping (MAXSTEP, a, b; dt=a/b) | `TIME 20 0.001 1.0` |
| `NONL` | Nonlinear convergence (N, eps, old) | `NONL -3 0.01` |
| `DMEG` | Load geometry + write results to .meg | `DMEG` |
| `BLAS` | Matrix precision (8=double, 4=single) | `BLAS 8` |
| `NCPU` | Parallel threads (0=auto) | `NCPU 0` |
| `NOGO` | Dry run (check input only) | `NOGO` |
| `EMFM` | Induced current (MID, NB3) | `EMFM 1 4` |
| `STAR` | Star/neutral connection | `STAR 1 2 3` |
| `HOLE` | Zero resistance near axis | `HOLE 0.001 3` |
| `PASS GENE` | Skip matrix assembly (reuse .mat) | `PASS GENE` |
| `PASS STEP ALL` | Skip MOME, compute post-blocks only | `PASS STEP ALL` |
| `PASS SW` | Speed optimization (0/1/2) | `PASS 2` |
| `HBGO` | Ignore B-H curve defects | `HBGO` |
| `EMGO` | Force induced current with current input | `EMGO` |
| `CLGO` | Force unbalanced coil subdivision | `CLGO` |
| `MACH` | DC superposition (2-pass) | `MACH problem1` |
| `FREQ` | Frequency [Hz] (MOMC only) | `FREQ 0 1000` |

### NONL Convergence Strategy
```
NONL  -3  0.01        # Newton-Raphson: 3 iterations (fast approach)
NONL   3  0.01        # Iterative correction: 3 iter (true solution)
NONL  -3  0.01        # NR again
NONL  10  0.01        # Iterative: 10 iter (final convergence)
```
- N > 0: successive substitution (if converges = true solution)
- N < 0: Newton-Raphson (fewer iterations but approximate)
- **Recommended:** alternate NR and iterative, finish with iterative
- `old` parameter (0 to <1): stabilization for non-converging iterative steps

### PASS Speed Optimization
| PASS | Matrix | LU | Use Case |
|------|--------|----|----------|
| (auto) | auto | auto | Default -- use normally |
| 0 | rebuild | rebuild | Full rebuild (baseline) |
| 1 | skip | rebuild | Nonlinear B-H with motion |
| 2 | skip | skip | Linear B-H (fastest) |

Validate PASS 1/2 results against PASS 0 for a few steps.

### SOL FIEL
| Keyword | Description |
|---------|-------------|
| `TIME` | Step selection (STEP1, STEP2, STEP3); supports `LAST`, `LAST-N` |
| `DMEG` | Read space nodes/MCO elements |
| `ELEM` | Compute surface field on poly elements |

### SOL FORC / SOL FORT
| Keyword | Description |
|---------|-------------|
| `SEL1` | Select single element face |
| `SELE` | Select by element range |
| `SELM` | Select by material range |
| `SELP` | Select poly elements |
| `DELE/DELM/DELP` | Deselect duplicate faces |
| `PART` | Face subdivision for accuracy |

### SOL FIXB (Lorentz Force on Coils)
| Keyword | Description |
|---------|-------------|
| `COIE` | Select coil by element (EID, NB3) |
| `COIM` | Select coil by material (MID, NB3) |

### SOL FIXA (Flux Linkage)
| Keyword | Description |
|---------|-------------|
| `FLUM` | Flux linkage target (MIDT, NB3, MIDS) |
| `DMEG` | Read flux line drawing elements |

FLUM with MIDT=MIDS: self-inductance. Different: mutual inductance.

### SOL BEAM
| Keyword | Description |
|---------|-------------|
| `STEP` | Max tracking steps |
| `TIME` | Time step size [s] |
| `STOP XYZ` | Bounding box [m] |
| `STOP ELEM` | Beam stop surface elements |
| `BMEG` | Read beam stop meg file |
| `DMEG` | Read particle file |
| `OMAG` | Output data thinning |
| `RELA` | Relativistic correction |
| `RUNG` | Solver method (1=Euler, 4=RK4) |
| `BINT` | Beam-beam interaction (space charge) |
| `GRAV` | Gravitational acceleration |
"""

# ============================================================
# .mei format
# ============================================================

MEI_FORMAT = """\
# .mei File Format (Mesh Script)

The `.mei` file is a **procedural mesh scripting language** processed by
IEmesh (mesh750.exe) to generate .meg files. It is NOT raw mesh data.

## Structure

```
MODEL <SOLVER> <SYMMETRY>    # Header (MAGIC/ELFIN/BEAM + T/K/R)
GSC <scale>                  # Global scale (e.g., 0.001 for mm->m)
IMA X                        # Image symmetry plane(s)
PUT A=10.0                   # Variable assignment

<source region commands>     # Coils, magnets, conductors, dielectrics

SPACE                        # Separator

<evaluation region commands> # Field evaluation surfaces (MCO/ECO)

STOP                         # End of script
```

## Variable Assignment
```
PUT A=10.0
PUT B=A*2.0
PUT C=_COS(30.0)*A         # Functions: _COS, _SIN, _TAN, _SQRT, _ABS, ...
```

## Two Regions
1. **Before SPACE**: Source geometry (material nodes MGR1/EGR1 + elements)
2. **After SPACE**: Evaluation geometry (space nodes MGR2/EGR2 + MCO/ECO elements)

## Comments
- `*` at line start
- `#` or `!` anywhere to end of line
- `{` / `}` block comments
- `==` at line start = end of data

## Control Flow
- `LOOP var = start end step` / `LOOPEND` (nest up to 10 levels)
- `IF (condition)` / `ELSE` / `ENDIF`
- `CASEIF var` / `CASE v1` / `CASEELSE` / `CASEEND`
- `GOTO label$` / `GOAL label$`

## Subroutines
```
SUBR name$(var1 var2)
  ...
RETURN

CALL name$(v1 v2)
```

## Configuration
- `INIT ELEM MAX` -- max node/element number (default 200000)
- `DEF EPS eps` -- merge tolerance (default 1E-10)
"""

# ============================================================
# .mei commands
# ============================================================

MEI_COMMANDS = """\
# Mesh Scripting Commands (.mei)

## AA -- Direct Element Creation (Auto-Mesh)

| Command | Shape | Description |
|---------|-------|-------------|
| `AA(O,x,y,z)` | -- | Set origin for next AA |
| `AA(N,N1)` | -- | Set origin by node |
| `AA(P,AX$)` | -- | Change coordinate axes (XY,XZ,YZ,...) |
| `AA(B,M1,M2,M3)` | -- | Set default division counts |
| `AA(LINEX,lx,Mx)` | Line | X-direction line |
| `AA(LINEY,ly,My)` | Line | Y-direction line |
| `AA(RECT,lx,ly,Mx,My)` | Quad | Rectangle |
| `AA(BLOCK,lx,ly,lz,Mx,My,Mz)` | Hex | Rectangular solid |
| `AA(BOX,lx,ly,lz,Mx,My,Mz)` | Quad | Box surface (6 faces) |
| `AA(BOX1,lx,ly,Mx,My)` | Line | Rectangle outline |
| `AA(ARCH3,rs,rl,t1,t2,lz,Mr,Mt,Mz)` | Hex | Arch/torus sector |
| `AA(ARCH2,rs,rl,t1,t2,Mr,Mt)` | Quad | Arch surface |
| `AA(ARCH1,r,t1,t2,Mt)` | Line | Arch line |
| `AA(RING1,r,Mt)` | Line | Ring (line) |
| `AA(RINGB,rs,rl,lz,Mt)` | Hex | Ring solid |
| `AA(SQRING,lx,ly,lz,lt,r,Mr)` | Hex | Square ring (coil) |
| `AA(SPHSUR,r,M)` | Tri | Sphere surface (20*M*M triangles) |
| `AA(CYLSUR,r,lz,Mr,Mt,Mz)` | Quad | Cylinder surface |
| `AA(CYLSOL,r,lz,Mr,Mt,Mz)` | Hex | Cylinder solid |
| `AA(DISK,r,Mr,Mt)` | Quad | Disk surface |
| `AA(HELICAL,r,Mt,zp,K,L)` | Line | Helix (K turns, L layers) |
| `AA(HOLE,En,r,t,Mr,M1,M2)` | Quad | Punch hole in quad |

## AN -- Node Creation

| Command | Description |
|---------|-------------|
| `AN(XYZ,x,y,z)` | Cartesian |
| `AN(RT2,r,t)` | Polar (2D) |
| `AN(RT3,r,t,z)` | Cylindrical |
| `AN(F3,N1,...,N2)` | Equal-division on line |
| `AN(FL,N1,N2,ND,p)` | Ratio-division on line |
| `AN(FC,NC,N1,N2,ND,p)` | Ratio-division on arc |
| `AN(PL,N1,N2,M,p)` | Equal-ratio by count |
| `AN(PC,NC,N1,N2,M,p)` | Arc by count |
| `AN(PB,N1,N2,N3,N4,M1,M2)` | Grid on quadrilateral |
| `AN(LLX,N1,N2,N3,N4)` | Line-line intersection |
| `AN(CNW,N1,N2,N3)` | Circle through 3 points |

## AE -- Element Creation by Node

| Command | Description |
|---------|-------------|
| `AE(N1,N1)` | Point element |
| `AE(N2,N1,N2)` | Line (2 nodes) |
| `AE(N3,N1,N2,N3)` | Triangle |
| `AE(N4,N1,N2,N3,N4)` | Quad |
| `AE(N6,N1,...,N6)` | Triangular prism |
| `AE(N8,N1,...,N8)` | Hexahedron |
| `AE(NP,SW,N1,...,Nn)` | 2D poly element |

## NB -- Material/Name Assignment

| Command | Description |
|---------|-------------|
| `NB(NAME,xxx)` | Set element label (MWL, MCL, MMB, etc.) |
| `NB(MID,n)` | Set material ID |
| `NB(PID,n)` | Set poly ID |
| `NB(CHAN,MM,old,new)` | Change material ID |
| `NB(CHAN,NM,name,mid)` | Rename by name |

## SE/SN -- Selection

| Command | Description |
|---------|-------------|
| `SE(ADD,ALL)` | Select all elements |
| `SE(ADD,E1,E2)` | Select element range |
| `SE(ADD,M,M1,M2)` | Select by material |
| `SN(ADD,ALL)` | Select all nodes |

## ME -- Mesh Operations

| Command | Description |
|---------|-------------|
| `ME(CPL,d,K)` | Copy elements (linear, K copies) |
| `ME(CPC,t,K)` | Copy elements (rotational) |
| `ME(CPM,SW)` | Copy elements (mirror) |
| `ME(EXL,d,M)` | Extrude linear (surface -> solid) |
| `ME(EXC,t,M)` | Extrude rotational |
| `ME(MVL,d)` | Move linear |
| `ME(MVC,t)` | Move rotational |
| `ME(DEL)` | Delete selected |
| `ME(B2S,Pn,Mn)` | Solid to surface conversion |

## BE -- Element Division

| Command | Description |
|---------|-------------|
| `BE(L2,M,p)` | Divide line |
| `BE(L4,M1,M2,p1,p2)` | Divide quad |
| `BE(L8,M1,M2,M3)` | Divide hexahedron |
| `BE(H4,SW)` | Quad to 2 triangles |

## TE -- Element Direction

| Command | Description |
|---------|-------------|
| `TE(REV)` | Reverse element direction |
| `TE(WLV)` | Align magnet to vector |
| `TE(WLR,SW)` | Align magnet radially |
| `TE(CLV)` | Align coil around vector (right-hand rule) |
| `TE(CME,En)` | Align normals by reference element |
| `TE(TRI)` | Convert quad to triangle |

## Other

| Command | Description |
|---------|-------------|
| `RB(MERGE,eps)` | Merge coincident nodes |
| `RB(EN)` | Renumber elements and nodes |
| `RB(STORE)` | Store current data |
| `CL(LONE)` | Clean orphan nodes |
| `VECX x1 y1 z1 x2 y2 z2` | Define direction vector |
| `VECN N1 N2` | Direction by nodes |
| `GE(N2,P)` | Geometric division ratio |
| `DMEG file.meg 1` | Import meg data |
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
MGSC 0.001                        # Scale (M=MAGIC, E=ELFIN, B=BEAM)
MIMA -X                           # Symmetry

MGR1  1  0  x  y  z               # Source nodes
<element records>                  # Source elements

MGR2  1  0  x  y  z               # Evaluation nodes
<element records>                  # Evaluation elements

BOOK  END                          # Footer
```

**Max line length: 120 half-width characters.**

## Scale Headers
- `MGSC <scale>` -- MAGIC coordinate multiplier to meters
- `EGSC <scale>` -- ELFIN
- `BGSC <scale>` -- BEAM

## Symmetry Headers
- `MIMA X` / `MIMA -X` -- same/opposite-pole mirror about YZ plane
- `MIMA ZC <N>` -- N-fold periodic rotational symmetry
- `MIMA -ZC <N>` -- opposite-pole periodic (N must be even)
- `EIMA` / `BIMA` -- same for ELFIN/BEAM

### Symmetry Availability
| Type | X | Y | Z | ZC |
|------|---|---|---|----|
| T (3D) | Yes | Yes | Yes | Yes |
| R (Axisym) | No | No | Yes | No |
| K (2D) | Yes | Yes | No | No |

### Force/Torque Multipliers for Symmetry
| MIMA | Fx | Fy | Fz | Tx | Ty | Tz |
|------|----|----|----|----|----|----|
| +/-X | 0 | 2 | 2 | 2 | 0 | 0 |
| +/-Y | 2 | 0 | 2 | 0 | 2 | 0 |
| +/-Z | 2 | 2 | 0 | 0 | 0 | 2 |
| +/-ZC n | 0 | 0 | n | 0 | 0 | n |

Multiply across all active symmetries. Zero means that component cancels.

## Periodic Boundary Conduction
`MLIN <eps>` -- eddy currents conduct across periodic boundary (normally insulated).

## Node Records
```
MGR1  <id>  0  <x>  <y>  <z>     # Source node
MGR2  <id>  0  <x>  <y>  <z>     # Evaluation node
```

## Element Records
```
<TYPE>  <id>  0|<PID>  <MID>  <n1>  <n2>  [<n3> ... <n8>]
```

Element TYPE = **PREFIX(2-3 chars) + NODECOUNT(1 digit) + DIM(1 char)**
"""

# ============================================================
# Element types
# ============================================================

ELEMENT_TYPES = """\
# ELF Element Type Naming Convention

Element names: **1st char (solver) + 2nd-3rd chars (type) + 4th char (nodes) + 5th char (DIM)**

## DIM Suffix
| Suffix | Symmetry | Description |
|--------|----------|-------------|
| `T` | 3D | Full three-dimensional |
| `K` | 2D | Infinite in Z, XY cross-section |
| `R` | Axisymmetric | Revolution about Z, XZ cross-section |

## Node Count and Shape
| Nodes | Shape |
|-------|-------|
| 8 | Hexahedron (solid) |
| 6 | Triangular prism (solid) |
| 4 | Quadrilateral (surface) |
| 3 | Triangle (surface) |
| 2 | Line segment |
| 1 | Point |

## MAGIC Element Prefixes

### Magnetic Body (soft magnetic, B-H curve in 1st quadrant)
| Prefix | Description | DOF (3D) | Notes |
|--------|-------------|----------|-------|
| `MMB` | General magnetic body | 6(hex),5(prism) | **Recommended** |
| `MMS` | Thick plate (no top/bottom charge) | 4,3 | MMB preferred |
| `MMT` | Thin plate (edge charge) | 1 | MMB preferred (v4.00+) |
| `MMP` | Poly element (closed polyhedra) | 1 | For arbitrary shapes |

### Magnet (permanent magnet, B-H in 2nd/3rd quadrant)
| Prefix | Description | DOF | Notes |
|--------|-------------|-----|-------|
| `MWL` | Direction by element (bottom=S, top=N) | 1 | |
| `MWV` | Direction by VEC1/VEC3 vector | 1 | **Recommended** |

### Coil (current-carrying)
| Prefix | Description | DOF | Notes |
|--------|-------------|-----|-------|
| `MCL` | Coil element | 0 | Direction = node order |

### Eddy Current
| Prefix | Description | DOF (3D) | Notes |
|--------|-------------|----------|-------|
| `MAB` | Non-magnetic conductor | 6,5 | Adjacent = conducts |
| `MAT` | Non-magnetic thin plate | 1 | |
| `MBB` | Magnetic conductor | 12,10 | B-H + eddy current |

**Holes in eddy current bodies:** Fill with high-resistivity MAB (target ~1 ohm*m).

### Evaluation / Other
| Prefix | Description | DOF | Notes |
|--------|-------------|-----|-------|
| `MCO` | Contour (field evaluation) | 0 | Space nodes, SOL FIEL |
| `MCM` | Maxwell stress evaluation | 0 | Space nodes, SOL FORT |
| `MJH` | Center field | 0 | Material nodes, SOL FMAG |
| `MCB` | Beam stop surface | 0 | 3D only |
| `M..` | Shape (display only) | 0 | No computation |

## ELFIN Element Prefixes

| Prefix | Description | DOF |
|--------|-------------|-----|
| `EMB` | Dielectric body | 6,5,4,3 |
| `EMS` | Dielectric thick plate | 4,3 |
| `EMT` | Dielectric thin plate | 1 |
| `EMP` | Dielectric poly | 1 |
| `ESC` | Electrode (known potential) | 1 |
| `EQL` | Fixed charge | 0 |
| `ESH` | Boundary (permittivity ratio) | 1 |
| `ESS` | Boundary (normal E both sides) | 1 |
| `ESN` | Boundary (normal E front side) | 1 |
| `ECO` | Space field evaluation | 0 |
| `ECM` | Maxwell stress evaluation | 0 |
| `ECB` | Beam stop surface | 0 |

## BEAM Element Prefixes
| Prefix | Description |
|--------|-------------|
| `BGR1` | Particle initial position |
| `BGR2` | Particle velocity direction |

## Label Availability by DIM

### T (3D): MMB MMS MMT MMP MWL MWV MCL MAB MAT MBB MCO MCM MJH MCB
### K (2D): MMB MMT MMP MWL MWV MCL MAB MAT MBB MCO MCM MJH
### R (Axisym): MMB MMT MMP MWL MWV MCL MAB MAT MBB MCO MCM

(3D adds MMS, MCB; 3D adds MJH for K; R has no MJH)
"""

# ============================================================
# B-H Curves
# ============================================================

BH_CURVES = """\
# B-H Curve Specification in ELF

## Unit Declaration
```
HBUN  <curve_id>  <unit_system>
```
- `OE/G` or just `OE` -- Oersted/Gauss (CGS)
- `A/M/T` or just `A/M` -- A/m and Tesla (SI)

## Scale Factors
```
HBSC  <curve_id>  <H_scale>  <B_scale>
```

## Curve Data Points
```
HBCU  <curve_id>  <H_value>  <B_value>
```
- Minimum 2 points, max 100 per curve
- Must be monotonically increasing in H
- **Magnetic bodies:** start from origin (0,0), 1st quadrant only
- **Magnets:** start from negative H (2nd/3rd quadrant), can extend to 1st

**Extrapolation warning:** If solution H exceeds max input H, extrapolates
from last 2 points. If slope >> mu_0, moment becomes excessive. Use
MaiEditor3 auto-extension to add safe extrapolation points.

## Curve Copy
```
HBEQ  <dest_id>  <source_id>        # Note: dest first
HBCP  <source_id>  <dest1> <dest2>  # Note: source first
```

## Curve-to-Material Assignment (per timestep)
```
HBCN  <MID>  <STEP>  <curve_id>
```
Different B-H per step enables demagnetization analysis.

## Recoil Permeability (Demagnetization)
```
HBRM  <curve_id>  <recoil_slope>  <Bmax>
```
- Bmax = max residual magnetization for irreversible process
- Knee point = intersection of recoil line through Bmax with B-H curve
- Below knee point: irreversible demagnetization occurs

## Example (Silicon steel in SI)
```
HBUN  1  A/M  T
HBCU  1      0.0    0.0
HBCU  1    200.0    1.0
HBCU  1    500.0    1.4
HBCU  1   5000.0    1.8
HBCU  1  50000.0    2.0
```

## Example (NdFeB magnet, Br=1.35T, Z-direction)
```
HBUN  1  A/M  T
HBCU  1  -700000.0  0.0
HBCU  1       0.0   0.90
```
Direction set by VEC1/VEC3 for MWV elements, or by element node order for MWL.
"""

# ============================================================
# MAGIC solver
# ============================================================

MAGIC_DOCS = """\
# MAGIC Solver (Magnetostatic BEM)

MAGIC solves magnetostatic field problems using the integral element method
(moment method / BEM). Supports static, transient, and sinusoidal analysis.

## Analysis Types

| Type | TIME | Elements | NONL | Description |
|------|------|----------|------|-------------|
| Single static | STEP=0 | MB,WL,CL | Yes | Basic magnetostatic |
| Multiple static | STEP>0 | MB,WL,CL (no AB/BB) | Yes | Batch parameter sweep |
| Transient | STEP>0 | MB,WL,CL,AB,AT,BB | Yes | Eddy current + motion |
| Sinusoidal (MOMC) | STEP>=0 | All | No | AC/frequency domain |

## Typical .mai Structure
```
USE  MAGIC  3.50
PRE
  HBUN 1 A/M T
  HBCU 1 0.0 0.0
  HBCU 1 200.0 1.0
  ...
  COI1 1 0 100 8 4 A
  AMP1 1 0 5.0
END
SOL MOME
  NCPU 0
  TIME 0 1 1              # Static (STEP=0)
  NONL -3 0.01
  NONL 10 0.01
  DMEG
END
SOL FIEL
  TIME 0
  DMEG
END
```

## Solver Chain
1. `SOL MOME` -- BEM matrix solve (mandatory first)
2. `SOL FIEL` -- B,H at evaluation points
3. `SOL FORC` -- Force on element surfaces (approximate)
4. `SOL FORT` -- Maxwell stress force (**more accurate** than FORC)
5. `SOL FIXB` -- Lorentz force on coils
6. `SOL FIXA` -- Flux linkage, vector potential
7. `SOL MAS` -- Charge file for BEAM
8. `SOL FMAG` -- Center field at MJH elements

## Eddy Current in MAGIC
- Use MAB/MAT (non-magnetic) or MBB (magnetic) elements
- Eddy currents are **automatic** when field varies with TIME steps
- No special header needed -- just include TIME steps and eddy current elements
- Adjacent elements with shared nodes: current conducts
- Double nodes or gaps: current insulated

## Example Categories
| Directory | DIM | Description |
|-----------|-----|-------------|
| MK/ | K (axisym) | Solenoids, pot magnets |
| MR/ | R (2D) | C-core, relay, actuator |
| MT/ | T (3D) | Complex 3D geometry |
| IPM/ | T | Interior permanent magnet motors |
| LscLl/ | K,T | Inductance computation |
| MOMC/ | T | AC complex analysis |
| magne/ | -- | Magnetization/demagnetization |
"""

# ============================================================
# ELFIN solver (electrostatic)
# ============================================================

ELFIN_DOCS = """\
# ELFIN Solver (Electrostatic BEM)

ELFIN solves **electrostatic** field problems using the integral element method.
It is the electric-field counterpart to MAGIC (magnetic field).

## Key Differences from MAGIC

| Feature | ELFIN | MAGIC |
|---------|-------|-------|
| Physical quantity | E, V, D | H, B |
| Material curve | D-E curve (EDSC/EDCU) | B-H curve (HBUN/HBCU) |
| Potential | VOL1 (electrode V) | -- |
| Charge | CHA1 (density), CHA2 (total/floating) | -- |
| Boundary | ESH (eps ratio), ESS/ESN (E normal) | -- |
| Element prefix | E (EMB, ESC, EQL, ...) | M (MMB, MCL, MWL, ...) |
| USE declaration | `USE ELFIN 3.50` | `USE MAGIC 3.50` |

## Applications
- Capacitor design
- High-voltage insulator analysis
- Electrode potential distribution
- Electrostatic force computation
- Electric field shielding

## Element Types
- **EMB** -- Dielectric body (D-E curve, 1st quadrant)
- **ESC** -- Electrode with known potential (VOL1)
- **EQL** -- Fixed charge element (CHA1)
- **ESH** -- Permittivity ratio boundary
- **ESS/ESN** -- Normal E-field boundary
- **ECO** -- Space field evaluation
- **ECM** -- Maxwell stress evaluation

## Typical .mai Structure
```
USE  ELFIN  3.50
PRE
  EDSC 1 1.0 1.0            # D-E curve scale
  EDCU 1 0.0 0.0            # D-E data points
  EDCU 1 1000.0 5000.0
  VOL1 1 0 100.0            # Electrode 1: +100V
  VOL1 2 0 -100.0           # Electrode 2: -100V
END
SOL MOME
  NCPU 0
  NONL -3 0.01
  NONL 10 0.01
  DMEG
END
SOL FIEL
  TIME 0
  DMEG
END
```

## D-E Curve Notes
- D* = epsilon* x E (dimensionless permittivity ratio)
- Monotonically increasing, minimum 2 points
- **CHA2** creates floating electrode: solver finds unknown potential and charge distribution

## Example Categories
| Directory | DIM | Description |
|-----------|-----|-------------|
| EK/ | K (axisym) | Axisymmetric electrostatic |
| ER/ | R (2D) | 2D cross-section |
| ET/ | T (3D) | Full 3D electrostatic |
"""

# ============================================================
# BEAM solver
# ============================================================

BEAM_DOCS = """\
# BEAM Solver (Charged Particle Tracking)

BEAM tracks charged particles through electromagnetic fields computed by
MAGIC and/or ELFIN via .mas charge files.

## Workflow
1. Run MAGIC with SOL MAS -> magnetic .mas file
2. Run ELFIN with SOL MAS -> electric .mas file (optional)
3. Create particle .meg file (BGR1 positions + BGR2 directions)
4. Create .mai with FILE, particle properties, SOL BEAM
5. Run BEAM -> trajectory .mag file

**2D model .mas files cannot be used** -- only 3D and axisymmetric.

## Particle File (.meg)
```
BOOK MEP 3.50
BGSC 0.001
BGR1  1  0  x  y  z       # Particle 1 position
BGR2  1  0  vx vy vz      # Particle 1 velocity direction
BGR1  2  0  x  y  z       # Particle 2 position
BGR2  2  0  vx vy vz
BOOK END
```

## .mai Structure
```
USE  BEAM  3.50
PRE
  FILE  MAGIC  problem_M   # Import magnetic field
  FILE  ELFIN  problem_E   # Import electric field (optional)
  CHAR  1  1.602E-19       # Electron charge [C]
  MASS  1  9.109E-31       # Electron mass [kg]
  VOLB  1  10000.0         # Acceleration voltage [V]
END
SOL BEAM
  STEP  1000               # Max tracking steps
  TIME  1.0E-10            # Time step [s]
  STOP XYZ -0.1 -0.1 0.0 0.1 0.1 0.5   # Bounding box [m]
  RELA                     # Relativistic correction
  DMEG
END
```

## Advanced Features
- **RELA** -- Special relativity (longitudinal + transverse mass)
- **BINT** -- Beam-beam interaction (space charge via line charges + iteration)
- **STOP ELEM** -- Beam stop surfaces (ECB4T/ECB3T) with reflection
- **GRAV** -- Gravitational acceleration
- **RUNG 4** -- 4th-order Runge-Kutta (default: proprietary high-speed method)
- **BIMA** -- Symmetry for beam-beam interaction trajectories

## Example Categories
| Example | Description |
|---------|-------------|
| BX01-BX08 | Various electron optics |
| BX*m | Modified versions |
| BX*e | Extended versions |
| BECB | Beam in ELFIN field |
"""

# ============================================================
# Treasure Box (quick reference tables)
# ============================================================

TREASURE_BOX = """\
# ELF Treasure Box (Quick Reference Tables)

## 1. Element-to-PRE Property Mapping

Which PRE headers apply to which element types:

| Property | Header | MB | MS | MT | MP | WL | WV | CL | AB | AT | BB |
|----------|--------|----|----|----|----|----|----|----|----|----|----|
| B-H curve | HBUN/HBCU | Req | Req | Req | Req | Req | Req | -- | -- | -- | Req |
| B-H copy | HBEQ | Opt | Opt | Opt | Opt | Opt | Opt | -- | -- | -- | Opt |
| Anisotropy | HBA1 | Opt | Opt | Opt | -- | -- | -- | -- | -- | -- | Opt |
| Lamination | HBA2 | Opt | Opt | Opt | -- | -- | -- | -- | -- | -- | Opt |
| Direction | VEC1/VEC3 | Opt | Opt | Opt | -- | -- | Req | -- | -- | -- | Opt |
| Complex mu | CMU1/I | Req* | Req* | Req* | Req* | Req* | Req* | -- | -- | -- | Req* |
| Turns | COI1 | -- | -- | -- | -- | -- | -- | Req | -- | -- | -- |
| Current | AMP1 | -- | -- | -- | -- | -- | -- | Sel | -- | -- | -- |
| Voltage | VOL1 | -- | -- | -- | -- | -- | -- | Sel | -- | -- | -- |
| Resistance | OHM1 | -- | -- | -- | -- | -- | -- | Sel | -- | -- | -- |
| Resistivity | OHM2 | -- | -- | -- | -- | -- | -- | -- | Req | Req | Req |
| Thickness | THIN | -- | -- | Req | -- | -- | -- | Opt | -- | Req | -- |
| Motion | MOV1/ORI1 | Opt | Opt | Opt | Opt | Opt | Opt | Opt | Opt | Opt | Opt |

Req = Required, Opt = Optional, Sel = Selectively required, * = MOMC only

## 2. Calculation Mode Compatibility

| Mode | MOME/MOMC | TIME | MB | WL | CL | AB/AT | BB | EMFM | NONL |
|------|-----------|------|----|----|----|-------|----|----|------|
| Single static | MOME | n=0 | OK | OK | OK | -- | tri | -- | Yes |
| Multiple static | MOME | n>0 | OK | OK | OK | X | X | X | Yes |
| Transient | MOME | n>0 | OK | OK | OK | OK | OK | Opt | Yes |
| Single sinusoidal | MOMC | n=0 | OK | OK | OK | OK | OK | Opt | No |
| Multiple sinusoidal | MOMC | n>0 | OK | OK | OK | OK | OK | Opt | No |

X = cannot use, tri = usable but pointless

## 3. SOL Block Execution Order

| After | FIEL | FORC | FORT | FIXB | FIXA | MAS | FMAG |
|-------|------|------|------|------|------|-----|------|
| MOME | Yes | Yes | Yes | Yes | Yes | Yes | Yes |
| MOMC | Yes | No | No | No | Yes | No | No |

## 4. Common Mistakes (Knotty Sample)

| # | Task | Required Elements | Required Settings |
|---|------|------------------|------------------|
| 1 | Eddy current | AB/AT/BB elements | TIME steps (automatic) |
| 2 | Joule heat | AB/AT/BB elements | TIME steps (automatic) |
| 3 | Induced current | CL elements | EMFM + OHM1 + TIME |
| 4 | Flux linkage | CL elements | SOL FIXA + FLUM |
| 5 | Vector potential | Space nodes | SOL FIXA + DMEG |
| 6 | Force on coils | CL elements | SOL FIXB + COIE/COIM |
| 7 | Force on iron | Non-CL elements | SOL FORC + SEL* |
| 8 | Maxwell stress | MCM elements | SOL FORT + SEL* |
| 9 | Plate thickness | MMT/MAT elements | THIN in PRE |

## 5. Current-Direction Subdivision (NB3)

EMFM, FLUM, COIE/COIM all have NB3 parameter (default 4) for subdividing
the coil element along the current direction. Increase when B varies
significantly along a single coil element's length. This is separate from
cross-section subdivision (NB1 x NB2).
"""

# ============================================================
# Sinusoidal Response (SOL MOMC)
# ============================================================

SINUSOIDAL_RESPONSE = """\
# Sinusoidal Response Analysis (SOL MOMC)

SOL MOMC replaces d/dt with jw for frequency-domain (complex linear) analysis.
Results are complex numbers convertible to amplitude and phase.

## Key Points
- **Nonlinear B-H NOT supported** -- use complex permeability (CMU1/CMU1I) instead
- Current/voltage/resistance become complex (AMP1I, VOL1I, OHM3I)
- After MOMC, only SOL FIEL and SOL FIXA available (no FORC/FORT/FIXB)
- Multiple frequencies can be swept (FREQ per STEP)

## Force and Joule Heat in AC Fields

Forces and Joule heat are proportional to the **square** of sinusoidal
quantities, so they oscillate at **2x the driving frequency**:

```
F(t) = A * cos(2*omega*t + P) + B
```

MagFilter2 extracts amplitude A, phase P, and center value B.

## ComplexTransfer2 Tool
Converts complex .mag file to pseudo-transient .mag (1 period / N steps)
for Wmap3 animation:
```
ComplexTransfer.EXE problem.mag 36
```

## Frequency Sweep Example
```
USE MAGIC 3.50
PRE
  CMU1  1 0 100.0           # Real permeability
  CMU1I 1 0 -10.0           # Imaginary permeability
  OHM2  2 0 1.0e-7          # Resistivity
  COI1  3 0 100 8 4 A
  AMP1  3 0 5.0             # Current amplitude (real)
  AMP1I 3 0 0.0             # Current amplitude (imag)
END
SOL MOMC
  TIME 3 1 1                # 4 frequencies (steps 0-3)
  FREQ 0 50
  FREQ 1 100
  FREQ 2 500
  FREQ 3 1000
  EMFM 3 4                  # Induced current
  DMEG
END
SOL FIEL
  TIME 0 3
  DMEG
END
```
"""

# ============================================================
# Anisotropy and Lamination
# ============================================================

ANISOTROPY_LAMINATION = """\
# Anisotropy (HBA1) and Lamination (HBA2)

## HBA1 -- Magnetic Susceptibility Anisotropy
```
HBA1  <MID>  <k>     # 0 < k < 1
```
- k = ratio of hard-direction to easy-direction susceptibility
- Easy direction set by VEC1/VEC3
- Hard-direction susceptibility = k x easy-direction susceptibility
- Applies to: MMB*, MBB* (not MMP*). Also MMS*, MMT* for easy-axis only.

## HBA2 -- Lamination Stacking Factor
```
HBA2  <MID>  <p>     # 0 < p < 1
```
- p = stacking factor (lamination fill factor)
- Models air gap magnetic reluctance in laminated steel
- Stacking direction set by VEC1/VEC3
- Applies to: MMB*, MBB* (not MMP*)

## Mutual Exclusion
**HBA1 and HBA2 cannot be used together on the same material.**

## Direction Setting
```
VEC1  <EID>  0  <vx>  <vy>  <vz>   # Per element
VEC3  <MID>  0  <vx>  <vy>  <vz>   # Per material
```

Used for:
1. MWV* magnet magnetization direction
2. HBA1 easy-magnetization direction
3. HBA2 lamination stacking direction

## For MOMC (Sinusoidal)
Same HBA1/HBA2 syntax, but base permeability from CMU1/CMU1I instead of B-H curve.

## Property Directions Summary

| Element | Direction from meg | HBA1 direction | HBA2 direction | Magnet direction |
|---------|-------------------|----------------|----------------|-----------------|
| MMB | -- | VEC1/VEC3 | VEC1/VEC3 | -- |
| MMS | thickness | VEC (in-plane) | NOT recommended | -- |
| MMT | thickness | VEC (in-plane) | NOT recommended | -- |
| MWL | bottom=S, top=N | -- | -- | node order |
| MWV | -- | -- | -- | VEC1/VEC3 |
| MCL | node order | -- | -- | -- |
| MBB | -- | VEC1/VEC3 | VEC1/VEC3 | -- |
"""

# ============================================================
# Steady-State Eddy Current (STED)
# ============================================================

STEADY_STATE_EDDY = """\
# Steady-State Eddy Current Motion (STED)

For rotating metal disks or translating metal plates in static magnetic fields.

```
STED  <MID>  <dx>  <dy>  <dz>  <dt>  <N>
```
- dx, dy, dz = displacement per step [m]
- dt = rotation angle per step [degrees]
- N = element number difference after one step overlap

## Mesh Requirements
1. Regular mesh division in motion direction
2. Element numbers consecutive perpendicular to motion
3. After one step, mesh overlaps exactly with original
4. Element number difference N must be uniform
5. Periodic rotational symmetry compatible
6. For linear motion: mesh must extend to where eddy current decays

## Example (3000 rpm disk, 5-degree mesh pitch)
```
STED 3 0 0 0 5 10
SOL MOME
  TIME 32 5 3000/60*360
  ...
END
```
"""

# ============================================================
# IPM Motor
# ============================================================

IPM_MOTOR = """\
# IPM Motor Analysis with MAGIC

## Ld/Lq Calculation Flow

1. **Run 1** (Motor1.mai): No coil current, rotate rotor -> flux linkage phi1 (magnet only)
2. **Run 2** (Motor2.mai): With sinusoidal coil current, rotate rotor -> phi2 (magnet + current)
3. **Post-processing:** phi3 = phi2 - phi1 (current-only contribution)
4. Fourier transform to extract fundamental d/q components
5. **Ld = phi3d / Id, Lq = phi3q / Iq**

## .model File Format
```
PartsM1 {
  RotAxis    = 0.0, 0.0, 0.0, 3     # Origin + axis (3=Z)
  Params_Rot = 1, 0.0, 60.0, 0.0
}
Calculations {
  TimeSeries = 180, 0.5, 180.0      # Steps, dt (deg), total angle
}
```

## Key Settings
- 3-phase Y-connection: use STAR for neutral point
- Magnets: HBRM with remanence + direction vectors
- Rotor/stator: MMS elements (laminated steel)
- Motion: ORI1 (Z-axis) + MOV1 (angle per step)
- Flux linkage: SOL FIXA with FLUM for each phase
- Torque: SOL FORT with MCM elements around rotor

## Phase Current Setup (MaiEditor3)
- U-phase: AMP1 amplitude, initial phase 225 deg (sin convention)
- V-phase: 105 deg, W-phase: 345 deg (120-degree offsets)

## Alternative LdLq Method
Three-phase-to-two-phase transform (Clarke) -> rotating coordinate transform
(alpha-beta -> d-q). Requires all 3 phase flux linkages.

## Example Files
- `magic/IPM/Motor1.mai` + `Motor1.meg` + `Motor1.model`
- `magic/IPM/Motor2.mai` + `Motor2.meg` + `Motor2.model`
- Excel: `LdLqの計算.xlsx` (Fourier), `LdLqの計算_座標変換.xlsx` (transform)
"""

# ============================================================
# Inductance computation
# ============================================================

INDUCTANCE_DOCS = """\
# Inductance Computation

## Self-Inductance (Basic)
1. Set linear permeability, 1A current, N turns
2. Run SOL MOME + SOL FIXA with FLUM
3. L = flux_linkage / current

## Two Definitions of Leakage Inductance

### Lsc (JIS Standard) -- Short-Circuit Inductance
Measured by short-circuiting secondary, computing primary impedance:
```
Lsc = X / (2*pi*f)
```
where X = Im(Z), Z = V/I (complex impedance at frequency f).

With zero resistance: Lsc = L1 * (1 - M^2 / (L1 * L2))

### Ll (IEEJ Definition) -- Leakage Inductance
Based on leakage flux energy:
```
Ll = L1 - 2*(n1/n2)*M + (n1/n2)^2 * L2
```
where Le1 = L1 - (n1/n2)*M, Le2 = L2 - (n2/n1)*M, Ll = Le1 + (n1/n2)^2*Le2

Energy method: Set I2 = -(n1/n2)*I1 to cancel main flux, then Ll = 2*W/I1^2

### Comparison
- Ll >= Lsc always; they converge as coupling approaches 1
- Lsc is easily measurable; Ll is theoretical

## 6 Calculation Samples (magic/LscLl/)

| Sample | Method | Description |
|--------|--------|-------------|
| 1 | Lsc | Basic: secondary short-circuit at 50Hz |
| 2 | Lsc | With current distribution (9 parallel sub-coils) |
| 3 | Lsc | Three-phase AC (Y-primary, delta-secondary) |
| 4 | Ll | Field energy method (cancel main flux) |
| 5 | Ll | From L1, L2, M via formula |
| 6 | Lsc | From L1, L2, M (resistance=0) |

## Key Formulas
| Quantity | Formula |
|----------|---------|
| Self-inductance | L = phi / I |
| Mutual inductance | M = phi_12 / I_1 |
| Lsc (JIS) | X / (2*pi*f) |
| Lsc (R=0) | L1 * (1 - M^2/(L1*L2)) |
| Ll (IEEJ) | L1 - 2*(n1/n2)*M + (n1/n2)^2*L2 |
| Ll (energy) | 2*W/I1^2 |
"""

# ============================================================
# Magnetization / demagnetization
# ============================================================

MAGNETIZATION_DOCS = """\
# Magnetization and Demagnetization Analysis

## Demagnetization by Opposing Field (3 steps)
- Step 0: No opposing field (operating point on B-H curve)
- Step 1: Opposing field applied (operating point moves left)
- Step 2: Opposing field removed (returns on recoil line)

If operating point passes knee point: **irreversible demagnetization**.
Recoil remanence ratio = L2/L1.

### ELF Implementation
- `HBRM <curve_id> <recoil_slope> <Bmax>` -- recoil line definition
- `HBCN <MID> <STEP> <curve_id>` -- different B-H per timestep

### Temperature Demagnetization (3 steps)
- Step 0: Normal temperature
- Step 1: Low/high temperature (shifted B-H curve)
- Step 2: Return to normal (inherit remanence ratio)

## Magnetization Analysis (3-Step Workflow)

### Step 1: Magnetizer Calculation
- Model magnet material as MMB (magnetic body) in magnetizer with coil
- Apply strong current -> compute B distribution
- Use initial magnetization curve (1st quadrant B-H)

### Step 2: Extract Magnet Model (MAGNE2.exe)
- Determines residual magnetization from Step 1 results
- Magnetization ratio k from k-H curve
- Converts MMB elements to MWV (magnet) elements
- Creates M.meg (shape) and M.mai (properties)
- Configuration: Step1-2.txt with problem names, CMID, 2nd-quadrant B-H, k-H curve

### Step 3: Magnetic Circuit Evaluation (VEC1C2.exe)
- Import M.meg into magnetic circuit via IEmesh DMEG command
- VEC1C2.exe corrects A.mai with magnetization vectors from M.mai
- Run final ELF/MAGIC calculation

### Multi-Step Extraction (MAGNE2.exe modes)
| Mode | Description |
|------|-------------|
| 1 | Per-element maximum timestep (recommended) |
| 2 | Global maximum timestep |
| 3 | Average maximum timestep |
| 4 | Directly specified timestep |
| 5 | Specified element's maximum timestep |

## Example Files
- `magic/magne/magnetize/CASE1-5` + `Samlpe/`
- `magic/magne/demagne/Sample01-04`
"""

# ============================================================
# Meshing Guidelines
# ============================================================

MESHING_GUIDELINES = """\
# Model Creation and Meshing Guidelines

## Element Shape Quality
- Hexahedra closer to cubes give better accuracy
- Since v4.00, thin/elongated hexahedra have improved accuracy
- **Obliquely crushed shapes degrade accuracy** -- lines connecting opposite
  face centers should be mutually perpendicular
- **Quads better than triangles** -- use triangles only to resolve diagonal issues

## Element Connectivity
- **Matching boundaries** between adjacent elements improve accuracy
- **Unmatched boundary** with face exposed to air: parasitic charges degrade accuracy
- **Gaps** between elements: increased magnetic reluctance -> underestimated moment
- **Overlapping** elements: negative reluctance -> overestimated moment
- **T-junctions** (intermediate nodes): cause insulation for eddy currents
  (acceptable for non-eddy-current analysis)

## Gap Region Meshing
- Face-to-face offset angle theta: **< 45 degrees**
- Element edge length at gap face: **<= gap distance**
- For moving models: mesh alignment at each step matters

## Coil Subdivision
- Solid coil replaced by NB1 x NB2 line currents
- Near-coil spatial field: increase NB1/NB2 so line current pitch <= half
  distance to space node
- Unbalanced NB1:NB2 ratio causes warning (ELF-W 901-903); use CLGO to override

## Eddy Current Element Holes
- Holes must be filled with high-resistivity MAB (target ~1 ohm*m)
- Fill creates equal/opposite currents that cancel, leaving edge current
- Bay-shaped openings connecting to exterior do not need filling

## Poly Elements (MMP)
- One poly = multiple pieces with same MID and PID
- Adjacent poly elements need boundary pieces from each
- Poly center must be inside the closed region; subdivide if not
- Ideal shape: sphere
"""

# ============================================================
# Nonlinear Convergence
# ============================================================

NONLINEAR_CONVERGENCE = """\
# Nonlinear Convergence

## Log File Columns
| Column | Meaning |
|--------|---------|
| ITERATION | Iteration count |
| CONVERGENCE | Maximum element error |
| MU*(NOW) | Relative permeability this iteration |
| MU*(NEXT) | Relative permeability for next |
| KIND | Element type with max error |
| ELEMENT | Element number with max error |
| MID | Material number with max error |

Iterative method lines: `++`, Newton-Raphson lines: `--`

## Convergence Criteria
- **Iterative:** |MU*(NOW) - MU*(NEXT)| / MU*(NOW) <= eps for ALL elements
- **Newton-Raphson:** |Bk - Bk-1| / Bk <= eps

## Troubleshooting Non-Convergence
1. **Increase iterations** -- if error is decreasing, more iterations may suffice
2. **Check mesh quality** -- extreme aspect ratios at max-error element
3. **Smooth B-H curve** -- check XT (differential susceptibility) in log:
   - RXT > 1.0: piecewise-linear is convex downward
   - RXT < 1.0: convex upward
   - Oscillating RXT: convergence is difficult; smooth data points
4. **Use stabilization coefficient** -- `NONL N eps old` with 0 < old < 1
5. **Extend B-H curve** -- `!` in moment output means H exceeded input range;
   extend curve until slope approaches mu_0

## Without NONL
B-H curve treated as **linear** using first 2 points only, even with 3+ points.
"""

# ============================================================
# Force Calculation Methods
# ============================================================

FORCE_METHODS = """\
# Force Calculation Methods

## Comparison

| Method | SOL Block | Elements | Accuracy | Use Case |
|--------|-----------|----------|----------|----------|
| Element surface force | FORC | Material elements | Lower | Quick estimate |
| Maxwell stress | FORT | MCM elements | **Higher** | Accurate force/torque |
| Lorentz force | FIXB | MCL coil elements | Good | Force on coils |

## SOL FORC (Element Surface Force)
- Computes force on selected faces of magnetic body and magnet elements
- **SOL FORT is more accurate** -- use FORC only for quick estimates
- Face subdivision with `PART <NB>` improves accuracy

## SOL FORT (Maxwell Stress)
- Requires MCM evaluation elements **in air** surrounding the target
- MCM elements must have consistent outward-pointing normals
- Finer mesh in regions of large/varying field
- More accurate than FORC because it evaluates in air (not on material surface)

## SOL FIXB (Lorentz Force on Coils)
- Computes F = IL x B at coil element positions
- NB3 subdivisions in current direction for accuracy
- For MCL2*/MCL1*: specify wire diameter via THIN

## MOMC Limitation
**After SOL MOMC, only FIEL and FIXA are available.**
For sinusoidal force, use the A/P/B oscillation parameters from MOMC results:
```
F(t) = A * cos(2*omega*t + P) + B
```
"""

# ============================================================
# Error messages
# ============================================================

ERROR_MESSAGES = """\
# ELF Error Messages

## Error Types
- **ELF-W(nnn):** Warning -- computation continues
- **ELF-E(nnn):** Error -- computation terminates mid-run
- **ELF-Q(nnn):** Error at end -- computation terminates

## Key Errors by Category

### Matrix Errors (301-323)
| Code | Message | Remedy |
|------|---------|--------|
| 311-319 | MATRIX = 0 / SINGULAR | Check B-H/D-E curve or element shape |
| 322-323 | DGETRS/SGETRS ERROR | PASS LU matrix is abnormal |

### File Errors (401-442)
| Code | Message | Remedy |
|------|---------|--------|
| 402,414-416 | MAT-FILE ERROR OR NO-KEY | Check license key connection |
| 408,413 | MAT-FILE NOT FOUND | Check .mat file in working folder |
| 409-410 | CANNOT WRITE IN DISK | Disk full -- free space or reduce model |
| 441 | NO LICENSE | Wrong HASP key for this software |

### Input Errors (501-599)
| Code | Message | Remedy |
|------|---------|--------|
| 501 | MAI FILE NOT FOUND | Check .mai file |
| 502 | USE LINE NOT FOUND | Fix USE block in .mai |
| 504-505 | MAI-FILE IS NOT MAGIC/ELFIN | Wrong solver for this .mai |
| 508 | SOL MOME/MOMC NOT FOUND | Add mandatory SOL block |
| 509 | MEMORY OVER | Reduce number of elements |
| 581 | UNDEFINED GRID | Check node coordinates in .meg |
| 582-584 | UNDEFINED VECT | Check VEC1 for element |
| 585-586 | UNDEFINED THIN | Add THIN for thin plate element |
| 587 | UNDEFINED AMPE | Add AMP1 for coil |
| 588 | UNDEFINED OHM2 | Add OHM2 for eddy current element |
| 589 | UNDEFINED TURN | Add COI1 for coil |

### Validation Errors (601-695)
| Code | Message | Remedy |
|------|---------|--------|
| 607 | COI1 DIVISION > 100 | Reduce NB1/NB2 to <= 100 |
| 673-674 | ONLY ONE DATA LINE | Add more B-H/D-E data points |
| 677-678 | H/E MUST INCREASE | Fix monotonicity |
| 679 | XT SHOULD BE PLUS | Fix B-H/D-E curve slope |
| 682 | HBCU SHOULD BEGIN (0,0) | Magnetic body curve must start at origin |
| 692-693 | SHAPE OF ELEMENT | Fix element shape |
| 743-744 | DO NOT USE HBA1 AND HBA2 TOGETHER | Remove one |

### Circuit/Limit Errors (701-744)
| Code | Message | Remedy |
|------|---------|--------|
| 701 | CURRENT MID > 1000 | Keep EMFM materials <= 1000 |
| 703 | STAR NEED EMFM | Add EMFM for star-connected coils |
| 716 | FREQ < 1E-10 | Increase frequency |
| 721 | CANNOT USE PASS GENE WITH MOTION | Remove PASS GENE |

### Warnings (901-906)
| Code | Message | Remedy |
|------|---------|--------|
| 901-903 | CHECK DIVISION (NB1,NB2) | Coil cross-section ratio unbalanced |
| 904-905 | XT SHOULD BE PLUS | Smooth B-H/D-E curve |
| 906 | DID NOT CONVERGE | Increase NONL iterations, refine mesh, smooth curve |

### BEAM Messages
| Message | Meaning |
|---------|---------|
| VELOCITY EXCEEDS LIGHT SPEED | Particle too fast -- check VOLB |
| TOO BIG SPEED CHANGE | TIME step too large |
"""

# ============================================================
# MEG Export from Cubit
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
Each Cubit block gets an ELF element label:
```
labels "1:MMB,2:MWL,3:MCL,4:MCO"
```

## Element Name = PREFIX(3) + NODE_COUNT(1) + DIM(1)
Node count from mesh type (quad=4, tri=3, line=2), DIM from symmetry selection.

## Available Labels by DIM

### T (3D)
MMB, MMS, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB

### K (Axisymmetric) / R (2D)
MMB, MMT, MMP, MWL, MWV, MCL, MCO, MAB, MAT, MBB
"""

# ============================================================
# Examples catalog
# ============================================================

EXAMPLES_CATALOG = """\
# ELF600 Example File Catalog

Located in `C:\\ELF600\\examples\\`.

## beam/ -- Charged Particle Beam Tracking
| Example | Description |
|---------|-------------|
| BASIC/ECB, BECB | Beam in ELFIN field |
| BX01-BX08 | Electron optics (+ m/e variants) |

## elfin/ -- Electrostatic Analysis
| Example | Description |
|---------|-------------|
| BASIC/EMB,EQL,ESC* | Basic electrostatic problems |
| EK/EK01-EK06 | Axisymmetric electrostatic |
| ER/ER01-ER05 | 2D cross-section |
| ET/ET01-ET06 | 3D electrostatic |
| WorkBook/EWork* | Tutorials |

## magic/ -- Magnetostatic Analysis
| Example | Description |
|---------|-------------|
| BASIC/MBCL,HCOIL,... | Simple magnet + coil models |
| MK/MK01-MK06 | Axisymmetric magnetostatic |
| MR/MR01-MR05 | 2D cross-section |
| MT/MT01-MT09 | 3D magnetostatic |
| IPM/Motor1,Motor2 | IPM motor Ld/Lq analysis |
| LscLl/Sample1-6 | Leakage inductance (6 methods) |
| MOMC/MOMCFJ | AC complex force/Joule heat |
| magne/magnetize/ | 5 magnetization cases |
| magne/demagne/ | 4 demagnetization samples |
| WorkBook/Work* | Tutorials |

## No TEAM problem benchmarks are included in the examples.
"""

# ============================================================
# IEmesh Tool
# ============================================================

IEMESH_TOOL = """\
# IEmesh -- Mesh Preprocessor

IEmesh is the integrated mesh creation tool for ELF. It consists of:

## Architecture
- **IEmesh.exe** -- Script editor with Command Guide, Command Filter, Hint View
- **mesh750.exe** -- Script-to-mesh compiler (.mei -> .meg). Requires HASP key.
- **Wmap3.exe** -- 3D visualization of .meg files

## Workflow
1. Create/edit `.mei` script in IEmesh editor
2. mesh750.exe compiles to `.meg`
3. Wmap3 visualizes for verification
4. Iterate

## Script Advantages
- Model creation procedure is saved (reproducible)
- Variable changes update dimensions globally
- Work can be divided among team members

## Command Categories (~18 two-letter families)
| Family | Purpose |
|--------|---------|
| AA | Direct element creation (auto-mesh shapes) |
| AN | Node creation |
| AE | Element creation by nodes |
| AF | Figure (geometry line/circle) creation |
| FN | Figure utilization (intersections, tangents) |
| NB | Names and material IDs |
| SE | Element selection |
| SN | Node selection |
| ME | Element editing (copy, move, extrude, delete) |
| MN | Node editing (copy, move, transform) |
| BE | Element division |
| TE | Element direction |
| GE | Division ratio |
| CL | Clear/delete |
| RB | Renumber |
| VEC | Vector definition |
| PUT | Variable assignment |
| DIM | Array declaration |

## Mathematical Functions
_ABS, _SQRT, _EXP, _LOG, _LOG10, _MAX, _MIN, _INT, _SIGN, _MOD,
_SIN, _COS, _TAN, _ASIN, _ACOS, _ATAN, _ATAN2

## Node/Element Information Functions
_XNODE, _YNODE, _ZNODE, _DIST, _RNODE, _TNODE, _NELEM, _KELEM,
_MNODE, _MELEM, _MMID, _KMID, _MID, _PID, _SEC, _SNC
"""

# ============================================================
# Wmap3 and MagFilter2
# ============================================================

TOOLS_DOCS = """\
# ELF Visualization and Post-Processing Tools

## Wmap3 -- 3D Visualization
- Opens .meg (geometry) and .mag (results) files
- Multiple views with independent rotation/zoom/pan
- Display modes: wireframe, shrink, solid
- Vector plots (field arrows with level colors)
- Contour plots (scalar fields with custom ranges)
- Animation across time steps
- Object selection by node, element, material, coordinates
- Data list for tabular numerical values

## MagFilter2 -- Result Extraction
- Extracts .mag results to CSV for spreadsheet analysis
- Physics tabs: magnetic flux, force, Joule heat, flux linkage, etc.
- Configurable output range, items, units
- Settings saveable for reuse

## MaiEditor3 -- Control File GUI Editor
- Creates/edits .mai files with mode-specific interfaces
- MAGIC mode: B-H curves, permeability, current/voltage, motion, force selection
- ELFIN mode: D-E curves, potential, charge, boundary conditions
- BEAM mode: particle properties, beam info, relativistic effects

## ELF/Bench -- Problem Manager
- Groups linked files as "problems" (.meg, .mai, .mag, .mao, .mas)
- Batch processing of multiple calculations
- Quick-launch of all ELF tools
- Drag-and-drop file registration

## ComplexTransfer2 -- MOMC Animator
- Converts complex .mag to pseudo-transient .mag (N steps per period)
- Enables Wmap3 animation of sinusoidal results
- Usage: `ComplexTransfer.EXE problem.mag 36`
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
    "treasure_box": TREASURE_BOX,
    "sinusoidal": SINUSOIDAL_RESPONSE,
    "anisotropy": ANISOTROPY_LAMINATION,
    "sted": STEADY_STATE_EDDY,
    "meshing": MESHING_GUIDELINES,
    "convergence": NONLINEAR_CONVERGENCE,
    "force_methods": FORCE_METHODS,
    "errors": ERROR_MESSAGES,
    "iemesh": IEMESH_TOOL,
    "tools": TOOLS_DOCS,
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
