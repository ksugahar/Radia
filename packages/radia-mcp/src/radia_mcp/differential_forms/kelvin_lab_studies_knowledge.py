"""
Kelvin Transform: Sugahara Lab Practical Case Studies (2020-2023).

Goes beyond the *theory* knowledge in `differential_forms.mathematica_recipes`
(topic 'kelvin') and `fem.gauge_open_boundary` (topic 'kelvin_transform') by
documenting WHAT THE LAB ACTUALLY BUILT IN COMSOL, the gotchas they hit,
and the workaround scripts that survived.

All seven case studies live under S:/COMSOL/88_ケルビン変換/ on the lab NAS.
They span:
  - Basic 2D / 3D Kelvin transform practice (2020-06)
  - Eddy current testing (ECT) with simulated earth (2021-06)
  - "Force a Periodic BC on every face" MATLAB+Java workaround (2021-09)
  - Wireless power transfer (WPT) with earth + PML (2021-10)
  - ECT eddy-current distribution post-processing (2021-11)
  - High-frequency Kelvin transform validity check via Hertz dipole sphere
    (2023-03; 2D / axisymmetric / 3D + offset / PML cross-checks)
  - A worked example mph (case07_a_40_offset_20.mph) ready to open

The lab convention is:
  - Inner ball:  V_inner    centered at origin, radius a
  - Outer ball:  V_outer    centered at (0, 0, 3a), radius a
  - Mapping:     y = (a^2 / |x|^2) x   then translate by (0, 0, 3a)
  - The lab uses 3a as the inter-center offset so the two balls are visibly
    separate in the GUI and easy to inspect.  Mathematically any offset > 2a
    works; 3a keeps them disjoint with a 1a air gap.
  - COMSOL implementation is via two GeneralExtrusion (genext) couplings
    binding each inner-sphere boundary face to its outer-sphere image.
  - For each E-field component the boundary condition is
    ``E0 = genext(E)`` on the outer side, which makes COMSOL solve
    the inner FEM and the outer (Kelvin-transformed) FEM as one coupled
    system through the genext map -- without any explicit weak-form
    transformation.

This file is the canonical record of what the lab tried, what worked,
what did not, and where the binary .mph files live so the next student
can re-open them and not start from zero.

Author trail: Sugahara Lab (Kindai University), 2020-2024.
ASCII / cp932 safe.  All Japanese paths are quoted verbatim from the NAS.
"""

OVERVIEW = r"""
# Kelvin Transform -- Sugahara Lab Practical Catalog (2020-2023)

The Kelvin transformation realizes an EXACT open-boundary condition for
Laplace-type problems (Group A: static / quasi-static / Coulomb gauge).
The lab has built and validated SEVEN COMSOL case studies that exercise
this technique on real engineering applications -- not toy spheres.

## What is "Kelvin transform" in COMSOL practice

A geometric trick: map the unbounded exterior of a ball
  E = { x in R^3 : |x| > a }
to a SECOND finite ball
  E' = { y in R^3 : |y - c| < a }
via the conformal inversion
  y = (a^2 / |x - c0|^2) (x - c0)  +  c

with weight factor
  lambda(y) = (a / |y - c|)^(n-2)/2   in R^n (n=3 here)

For the scalar Laplace equation:
  Delta u = 0       in V_inner U E        (original problem)
becomes
  Delta u' = 0      in V_inner U E'       (Kelvin problem)
with u' = lambda * u (the conformal multiplier).

For the magnetic vector potential A in Coulomb gauge, the same identity
extends: A satisfies Laplace's equation in air, so A' = lambda * A solves
the Kelvin-mapped problem.  See
``differential_forms.mathematica_recipes('kelvin')`` for the symbolic
derivation of the factor for forms of degree k = 0, 1, 2.

## Lab convention (used in ALL seven case studies)

Pick the inner ball of radius a centered at origin.  Pick the outer
("Kelvin-image") ball of radius a centered at (0, 0, 3a).  The center
offset 3a is the LAB STANDARD; any value > 2a would work geometrically,
but 3a is the visual sweet spot.

The Kelvin map composes inversion + translation:
  inner point  p  ->  outer image  q = (a^2 / |p|^2) p + (0,0,3a)

The boundary identification is enforced via TWO COMSOL ``GeneralExtrusion``
operators:
  genext1: from the inner-air boundary at |x|=a to the outer-air boundary
           with ``dstmap = (x, y, z + 3a)``
  genext2: from the outer-air boundary at |y-c|=a back to the inner side
           with ``dstmap = (x, y, z - 3a)``

The physics BC on the outer boundary then reads
  E0 = genext1(E)        (electromagnetic waves)
  A0 = genext1(A)        (induction currents)

COMSOL solves the two FEM domains as a single coupled system.

## Seven Case Studies (chronological)

| # | Year   | Folder                                          | Application                            | Status |
|---|--------|-------------------------------------------------|----------------------------------------|--------|
| 1 | 2020-06 | 2020_06_12_Kelvin変換の練習                     | Basic 2D + 3D + 3D-quadrant practice   | works  |
| 2 | 2021-06 | 2021_06_18_大地を模擬したケルビン変換@ECT       | Eddy current testing with earth        | works  |
| 3 | 2021-09 | 2021_09_01_各面に対して強制的に周期境界条件を課すスクリプト | Forced-periodic-BC MATLAB script | works  |
| 4 | 2021-10 | 2021_10_22_大地を模擬したケルビン変換@WPT       | Wireless power transfer + earth + PML  | works  |
| 5 | 2021-11 | (top-level) 2021_11_06_渦電流探傷_渦電流分布用.mph + ECT_01/ | ECT eddy current distribution post   | works  |
| 6 | 2023-03 | 2023_03_18_高周波のケルビン変換_ダイポール球    | HF Kelvin validity vs Hertz dipole     | works  |
| 7 | 2023-12 | case07_a_40_offset_20.mph (and 3D/case07_*)    | Worked example: a=40 with offset       | works  |

## Why this matters for Radia

Radia ships ``calc_fem_kelvin.py`` (Layer 4) wrapped by ``radia_ih.py``
(Layer 3) which use the SAME mathematical trick on NGSolve -- but with
the C++ Cubit-side identification of inner/outer face pairs instead of
COMSOL's GeneralExtrusion.  The seven lab case studies provide:

  1. VALIDATION TARGETS.  COMSOL has shipped Kelvin since 5.6; comparing
     a Radia/NGSolve Kelvin solve to a COMSOL Kelvin solve on the SAME
     geometry isolates whether a discrepancy is in the formulation
     (Kelvin shell physics) or in the implementation (face-pair
     identification, mesh, gauge).
  2. ENGINEERING PRECEDENT.  Earth simulation (ECT, WPT) is exactly the
     class of problem where PML fails at low frequency and BEM is too
     expensive -- Kelvin is the right tool.  The case studies prove the
     lab has shipped this in 2021.
  3. PITFALL CATALOG.  Each study hit some issue (Brick size limit,
     PML interaction, periodic BC ordering, ...) that is worth recording
     so the next student doesn't re-discover it.
  4. RECIPE LIBRARY.  The forced-periodic-BC MATLAB script in study #3
     is a portable pattern: enumerate all boundaries, pair them by
     position, batch-create the BC.  Same logic applies to Kelvin face
     pairing -- a useful template for NGSolve-side automation.

## Quick links

- LAB_FILES topic: every script and mph path with size/date
- ECT_EARTH topic: case study #2 deep-dive (induction currents + earth)
- WPT_EARTH topic: case study #4 deep-dive (RF + PML + earth + Kelvin)
- HF_DIPOLE_SPHERE topic: case study #6 (when does Kelvin break at HF?)
- FORCED_PERIODIC_BC_SCRIPT topic: study #3 MATLAB walkthrough
- LESSONS_LEARNED topic: cross-cutting pitfalls
- ROBIN_LIKE_PATTERN topic: the genext = closest gotcha
- CROSS_REFERENCE: how each case relates to the radia panel pipeline
"""


KELVIN_BASIC_PRACTICE = r"""
# Case Study #1 (2020-06) -- Basic 2D + 3D Kelvin practice

Folder: S:/COMSOL/88_ケルビン変換/2020_06_12_Kelvin変換の練習/
Lead files:
  Kelvin_2D.mph         2.0 MB   2D circle (electrostatics)
  Kelvin_3D.mph         7.0 MB   3D sphere (electrostatics, with quarter sym)
  Kelvin_3D_4Q.mph      5.9 MB   3D quarter-symmetric sphere
  Kelvin_3D.m          79 kB    MATLAB rebuild script
  練習.mph              3.7 MB   sandbox / scratch
  周期境界練習.mph      1.2 MB   periodic BC practice (preceded the
                                  forced-script of case #3)

## Geometry pattern (3D)

```matlab
% In Kelvin_3D.m at lines 20-26
model.component('comp1').geom('geom1').create('sph1', 'Sphere');
   % default radius 1, default center (0,0,0)
model.component('comp1').geom('geom1').create('sph2', 'Sphere');
model.component('comp1').geom('geom1').feature('sph2').set('pos', [0 0 2.1]);
   % center (0, 0, 2.1) -- the outer Kelvin ball
model.component('comp1').geom('geom1').create('sph3', 'Sphere');
model.component('comp1').geom('geom1').feature('sph3').set('r', 0.1);
   % inner source (a tiny PE conductor at origin)
```

Note: 2020 used offset 2.1*a (just barely > 2a, so the balls touch with
0.1 gap).  By 2021-06 (ECT) the lab moved to 3*a for visual clarity, and
by 2023 (HF dipole) firmly to 3.0*a in all scripts.  RECOMMENDATION:
always use 3*a for new lab work.

## Boundary identification (4 LinearExtrusion couplings)

```matlab
% Lines 27-67 build FOUR LinearExtrusion couplings linext1..linext4
% Each pairs one face of the inner sphere with the equivalent face
% on the outer sphere using vertex-by-vertex correspondence
model.component('comp1').cpl.create('linext1', 'LinearExtrusion');
% ... and 3 more (linext2, linext3, linext4) by duplication
```

INSIGHT: 2020 used ``LinearExtrusion`` (vertex-by-vertex), not the
much simpler ``GeneralExtrusion`` with ``dstmap`` (1 line of map
expression).  By 2021-06 the lab had switched to GeneralExtrusion.
This is the first lesson: GeneralExtrusion is FAR cleaner than
LinearExtrusion for Kelvin.

## Physics (electrostatics)

```matlab
% Line 73-79
model.component('comp1').physics.create('es', 'Electrostatics', 'geom1');
model.component('comp1').physics('es').create('pot1', 'ElectricPotential', 2);
model.component('comp1').physics('es').feature('pot1').selection.set([24]);
model.component('comp1').physics('es').feature('pot1').set('V0', 'linext1(V)');
% V0 = linext1(V): the outer face's potential is whatever the inner
%                  face's solution evaluates to via the Kelvin map
```

CRITICAL POINT: there is NO explicit weight factor lambda in this
boundary condition.  COMSOL is enforcing CONTINUITY of V at the inner /
outer interface, AND the geometric map (LinearExtrusion / GeneralExtrusion)
identifies the inner sphere boundary at |x|=a with the outer sphere
boundary at |y-c|=a.  The two FEM solves stitch together at that
interface and the outer solve plays the role of u' = u in the
mathematical Kelvin transform.

For SCALAR LAPLACE (degree-0 form), the conformal weight is lambda^((n-2)/2)
= lambda^(1/2) in 3D, but COMSOL's standard practice is to ABSORB this
weight into the metric of the outer domain implicitly via the genext map.
Verification: at the inner boundary V_inner = V_outer (by construction),
so the BC is Dirichlet-matching, not Robin-with-weight.

For VECTOR (A field), it is more subtle -- see WPT_EARTH topic.

## Outputs (visualization only)

The 2020 work produced .vtu visualization files:
  case1.vtu        419 kB    slice of solution
  case2.vtu        5.4 MB    isosurface
  isosurface.vtu   34 MB     dense isosurface
  streamline.vtu   16 MB     streamlines

Plus a PNG note image:
  円周上の値の引継ぎ.png   1.1 MB    sketch of the "value handover on
                                      the circle" -- the visual proof
                                      that the inner BC = outer BC at
                                      the sphere interface.

## Status

VALIDATED 2020-06 for electrostatics.  Re-used as building block for all
later studies.  The .m file is the canonical rebuild script.
"""


ECT_EARTH = r"""
# Case Study #2 (2021-06) -- Eddy Current Testing with Earth

Folder: S:/COMSOL/88_ケルビン変換/2021_06_18_大地を模擬したケルビン変換@ECT/

This study extends the basic practice (case #1) to the Magnetic Field
formulation (COMSOL ``InductionCurrents``, NOT ``ElectromagneticWaves``)
and adds a half-space earth model that simulates the ground under an
eddy-current probe.  This is the production lab application: ECT for
testing buried infrastructure where the half-space below the probe is
finite-conductivity soil.

## Lead files

| File                                    | Size      | Role                              |
|-----------------------------------------|-----------|-----------------------------------|
| main.m                                  |  18.6 kB  | canonical rebuild script          |
| Kelvin_AC_coil.mph                      |  41.4 MB  | the lab production run            |
| COIL_case1_30.mph                       |  54.7 MB  | a=30mm parametric run             |
| case2_30.mph                            |  22.9 MB  | second comparison case            |
| Kelvin_2D_mag_練習.mph                  |   2.0 MB  | 2D magnetics warm-up              |
| Kelvin_2D_wav_練習.mph                  |   0.3 MB  | 2D wave-eq warm-up                |
| Case3_ゲージ固定無_無限遠.mph           |  13.0 MB  | "gauge-free + infinite domain"    |
| 練習_Case1/Case1_ゲージ固定有.mph       |  16.4 MB  | "gauge-fixed" comparison          |
| 練習_Case2/Case2_ゲージ固定無.mph       |  15.6 MB  | "gauge-free" comparison           |
| 練習_Case3/Case3_ゲージ固定無_無限遠.mph| 13.0 MB   | "gauge-free + infinite" -> ★      |
| コイルのみ/case1_matlab.m               |   ~30 kB  | coil-only reference (no earth)    |
| 2021_06_16_ケルビン変換.xlsx            |  16 kB    | boundary-pair mapping spreadsheet |

## Physics

```matlab
% From コイルのみ/case1_matlab.m
model.component('comp1').physics.create('mf', 'InductionCurrents', 'geom1');
model.component('comp1').physics('mf').create('gfa1', 'GaugeFixingA', 3);
model.component('comp1').physics('mf').create('coil1', 'Coil', 3);
model.component('comp1').physics('mf').feature('coil1').selection.set([3 4]);
model.component('comp1').physics('mf').create('mp1', 'MagneticPotential', 2);
% mp1 = magnetic potential BC -- this is where we'll attach genext1(A)
model.component('comp1').physics('mf').create('constr1', 'PointwiseConstraint', 2);
model.component('comp1').physics('mf').create('pmc1', 'PerfectMagneticConductor', 2);
% pmc1 = perfect magnetic conductor on the outer Kelvin ball (E field = 0
% tangential, no field exits)
```

The lab tested THREE gauge variants:
  Case1: gauge fixed (explicit gfa1)
  Case2: gauge-free (no gfa1)
  Case3: gauge-free + infinite (no gfa1 + Kelvin to outer)  <- production

The "infinite" Case3 is the Kelvin-transformed version of Case2.
This is the actual ECT production setup the lab uses.

## Geometry

```matlab
% From コイルのみ/case1_matlab.m, lines ~30-55
% Coil: cylindrical shell (r=10..15 mm, h=20 mm) centered at z=0
cyl1 r=15, h=20, z=-10        % outer cylinder of coil shell
cyl2 r=10, h=20, z=-10        % inner cylinder of coil shell
dif1 = cyl1 - cyl2            % coil annulus

% Inner Kelvin ball (the air domain around the coil)
sph1 r=30, center (0,0,0)

% Outer Kelvin image ball
sph2 r=30, center (0,0,3a)    % a=30mm so center is at z=90mm
sph3 r=3,  center (0,0,3a)    % a tiny artifact ball (PE conductor sphere)
dif2 = sph2 - sph3            % outer Kelvin shell with conductor hole

% Workplane partition (a Y-Z plane that cuts the coil for symmetry)
wp1 (y-z plane)
par1 = partition dif1 with wp1
```

PITFALL: The lab spent significant time deciding the order of operations
between webcut, partition, intersect, and subtract.  ENTITY IDS DRIFT.
The .m file shows the FINAL working sequence; trying to substitute
operations changes the resulting entity numbering and breaks the
``mp1`` boundary selection (which is hardcoded to entity IDs).  This is
a perfect example of the "Journal File Portability Policy" in CLAUDE.md:
hardcoded IDs are brittle; identify by geometry where possible.

## Earth slab

The earth slab is a thin brick under the coil:
```
brick x 100 y 100 z 5    % 100x100 mm earth slab, 5 mm thick
move z -22.5             % put it below the coil (coil bottom at z=-10)
```

DIMENSIONAL LIMIT discovered by the lab (see file
``練習_Case3/coreform/※問題ありで休止中.txt``):
  "大地 Brickの、xy方向に110が最大で、それ以上がKelvin変換出来なかった"
  (For the earth Brick, xy=110 was the maximum -- beyond that the
   Kelvin transform did not work.)

This is a hint that for very large earth slabs, the Kelvin sphere radius
'a' must scale up correspondingly: the inner Kelvin ball must contain
the entire earth slab geometry.  If the earth poked through |x|=a, the
boundary identification map ``dstmap = (x, y, z+3a)`` started yielding
geometric inconsistencies (because part of "the earth" is OUTSIDE the
ball that is supposed to be Kelvin-mapped).

LESSON: Pick the Kelvin ball radius so that ALL conductive / magnetic
content is strictly inside it.  Air only outside.  Then Kelvin is exact.

## Boundary pair identification (the .xlsx file)

``2021_06_16_ケルビン変換.xlsx`` is the lab's hand-built spreadsheet
listing every inner face number and its outer pair.  In 2021-06 this
mapping was BY HAND, derived from visual inspection of the COMSOL GUI.

By 2021-09 this had been replaced by the Set_periodic_COMSOL.m script
that auto-pairs all boundaries by centroid (case study #3) -- see
FORCED_PERIODIC_BC_SCRIPT topic.  That script is the right way; the
.xlsx is the historical record of how painful the manual approach was.

## Status

VALIDATED 2021-06.  Lab production setup for ECT with simulated ground.
The .m file is the canonical rebuild; the .xlsx is for educational
reference (how NOT to do it).
"""


WPT_EARTH = r"""
# Case Study #4 (2021-10) -- Wireless Power Transfer + Earth + PML + Kelvin

Folder: S:/COMSOL/88_ケルビン変換/2021_10_22_大地を模擬したケルビン変換@WPT/

This is the LARGEST and most complex case study.  It models WPT at
12 MHz between two coupled coils, with:
  - a finite-conductivity earth slab BELOW the lower coil
  - a PML cap centered on the upper Kelvin image (to absorb residual
    outgoing waves at HF)
  - the Kelvin transformation maps the unbounded exterior to a finite
    second ball ABOVE the upper coil
  - the genext-based outer-sphere BCs glue inner and outer FEM solves

It is the unique study in the lab that COMBINES Kelvin AND PML.  At
12 MHz the geometry is electrically thick enough that pure Kelvin
(quasi-static assumption) leaks; PML around the inner ball cleans up
the residue.

## Lead files (chronological)

| File                                                              | Size      | Role                                  |
|-------------------------------------------------------------------|-----------|---------------------------------------|
| Kelvin.sat                                                        |   1.7 MB  | the lab's reusable Kelvin SAT model   |
| Kelvin-earth_300R.py                                              |   1.4 kB  | Cubit Python: build earth + Kelvin    |
| Kelvin_RF_3D_genext_PML_earth.java                                |  71 kB    | the COMSOL Java rebuild script (★)    |
| NG_Kelvin_RF_3D_genext_PML_earth.mph                              |  18 MB    | "NG" = not-good intermediate          |
| OK_Kelvin_RF_3D_genext_PML_earth.mph                              |  33 MB    | "OK" = working version (★)            |
| OK_Kelvin_RF_3D_genext_PML_earth_wireless.mph                     |   2.4 GB  | full WPT production run (★★)          |
| Comsol用_Kelvin対面.xlsx                                          |  10.6 kB  | inner/outer face-pair mapping         |
| 00_coreform/2021_08_28_Kelvin-earth-air(tris).py                  |  ~2.3 kB  | tri-mesh on coreform side            |
| 00_coreform/2021_08_28_Kelvin-earth-air(quads).py                 |  ~2.3 kB  | quad-mesh variant                    |
| 00_coreform/2021_08_25_Full-model(for_sat).jou                    |   2.0 kB  | Cubit .jou: full-model SAT export    |
| airplane_antenna_crosstalk.pptx                                   |   2.6 MB  | application slides                   |

## What "OK" and "NG" mean

``NG_Kelvin_RF_3D_genext_PML_earth.mph`` (Aug 22 2021, 18 MB) was an early
version where the genext mapping had a coordinate-sign issue (the upper
Kelvin image was wired with z+up instead of z-up for the inverse direction).
The symptoms: the field in the outer Kelvin ball did not decay smoothly to
zero at the inner-outer interface, indicating the BC was not enforcing the
correct Kelvin image of the inner solution.

``OK_Kelvin_RF_3D_genext_PML_earth.mph`` (Oct 1 2021) fixed this by
explicitly creating EIGHT genext operators (genext1 .. genext8): one
for each Cartesian component of E (Ex, Ey, Ez) AND its inverse, for
EACH of the four boundary sets (sph1/sph3, sph3/sph1, etc.).  Source
of truth: ``Kelvin_RF_3D_genext_PML_earth.java`` lines 143-168.

The pattern (from the java file):
```java
model.cpl.create("genext1", "GeneralExtrusion");
// ... 7 more
genext1.dstmap = {x, y, z+up}     // forward: inner -> outer
genext2.dstmap = {x, y, z-up}     // inverse: outer -> inner
genext3.dstmap = {x, y, z+up}     // another face pair, forward
genext4.dstmap = {x, y, z-up}     // ... inverse
genext5.dstmap = {x, y, z+up}     // genext5.method = "closest"  <- key
genext6.dstmap = {x, y, z-up}
genext7.dstmap = {x, y, z+up}     // genext7.method = "closest"  <- key
genext8.dstmap = {x, y, z-up}

// Apply on each E field component
efb2.E0 = {genext1(Ex), genext1(Ey), genext1(Ez)}
efb3.E0 = {genext2(Ex), genext2(Ey), genext2(Ez)}
// ... etc
```

## Why ``method = "closest"`` matters (Robin-like pattern gotcha)

For some genext operators (genext5, genext7 in the .java file), the lab
set ``method = "closest"`` instead of the default exact-match.  The
default tries to find the source point exactly under the dstmap, and
FAILS when the inner and outer meshes are not topologically identical
(which happens when one side has been Cubit-meshed and the other is
Netgen-tet, or when imprint/merge introduced extra small faces).
``method = "closest"`` falls back to nearest-neighbor lookup.

PITFALL: this is convenient but silently masks mesh-mismatch bugs.  If
the inner-outer meshes are SUPPOSED to be conforming (you used CopyFace
or auto-imprint), and you still need "closest", something is wrong --
investigate.  Lab convention: use "closest" only on the two genext
operators that fail to converge with the default; explicitly comment
which two and why.

## Kelvin geometry (from the .java file lines 26-75)

```java
sph1 r=a,       center (0,0,0)    // inner air ball
sph2 r=0.1*a,   center (0,0,0)    // tiny inner source ball
dif1 = sph1 - sph2                // inner air with source hole

sph3 r=a,       center (0,0,up)   // outer (Kelvin image) air ball
sph4 r=0.2*a,   center (0,0,up)   // outer PML CAP
sph5 r=0.05*a,  center (0,0,up)   // outer artifact subtract
dif2 = (sph3 + sph4) - sph5       // outer ball + PML annulus
```

up = 2.1 * a (= 2.1*0.3 = 0.63 m, see param "a", "up" lines 18-19 in java)

NOTE: 2021-10 still used up = 2.1*a (legacy from case study #1).  Newer
work (2023) uses 3.0*a.  Both work; 3.0*a is more visually robust.

Then a Kelvin SAT file is imported and intersected:
```java
imp1 = import Kelvin.sat
int1 = intersection(dif2, imp1, keep=true)
```

The Kelvin.sat (1.7 MB) is a pre-built ACIS solid representing the
Kelvin-image of the SAT geometry through inversion.  Pre-built so that
the heavy geometric inversion is done once in Cubit (where the
Kelvin-earth_300R.py Python script computes the image vertex-by-vertex)
rather than per-COMSOL-rebuild.

## The Kelvin-earth Cubit Python (key recipe)

From Kelvin-earth_300R.py (top-level WPT folder), this is the LAB'S
CANONICAL Python recipe for building a Kelvin image of an arbitrary
Cubit mesh, with translation by up=3a (or 2.1a in earlier versions):

```python
from numpy import *
import sys, re

FileName = 'o:/temp.pat'
a = 300        # sphere radius in mm (earth-scale)
up = 3*a       # vertical offset of upper Kelvin image (3a is standard)

def Kelvin(block_id):
  cubit.cmd('export patran "%s" block %d dimension 3 overwrite'
            % (FileName, block_id))
  fid = open(FileName,'r')
  strLines = fid.readlines()
  for n in range(len(strLines)):
    sline = re.split(r"\s+|\n$", strLines[n])
    if (sline[0] == '01'):
      # next line is the (x,y,z) vertex coordinate triple
      sline = re.split(r"\s+|\n$", strLines[n+1])
      o1 = array([float(sline[1]), float(sline[2]), float(sline[3])])
      r1 = linalg.norm(o1)
      # The Kelvin map: y = (a^2/r^2) x  then translate by (0,0,up)
      p1 = array([(a/r1)**2 * o1[0],
                  (a/r1)**2 * o1[1],
                  (a/r1)**2 * o1[2] + up])
      strLines[n+1] = ' %15e %15e %15e\n' % (p1[0], p1[1], p1[2])
  fid.close()
  fid = open(FileName,'w')
  for n in range(len(strLines)):
    fid.write(strLines[n])
  fid.close()

# usage:
cubit.cmd('reset')
cubit.cmd('brick x 600 y 600 z 50')     # build the earth slab
cubit.cmd('move Volume 1 z -225 include_merged')
cubit.cmd('volume all size 25')
cubit.cmd('mesh surf all')

n = 1
cubit.cmd('block %d add surface all' % n)
Kelvin(n)                                 # write the Kelvin-mapped vertices
cubit.cmd('reset')
cubit.cmd('import patran mesh geometry "%s" feature_angle 135.00 ' % FileName)

# Rebuild surfaces from the now-curved imported mesh
for n in range(1,7):
  cubit.cmd('create surface net from mapped surface %d heal ' % n)
cubit.cmd('delete body 1 to 6')
cubit.cmd('create volume surface 7 to 12 heal')
```

This pattern (export patran -> rewrite vertex lines -> reimport) is a
useful template even outside Kelvin: any geometric transformation that
acts vertex-by-vertex can be applied this way without needing direct
ACIS/OpenCASCADE access.

## Status

VALIDATED 2021-10 for WPT at 12 MHz with earth.  The 2.4 GB production
.mph is the lab's reference run.  Re-opening it requires a workstation
with > 32 GB RAM (the run mesh is fine).

PML + Kelvin combination is the LAB UNIQUE TRICK for this study.  In
plain Radia/NGSolve practice, you would NOT combine the two -- you would
pick one based on the regime (Kelvin for static / MQS; PML for radiation).
The lab used both here because at 12 MHz the geometry is borderline
(skin depth in earth ~ a few cm, comparable to coil size; some leakage
is propagating wave, some is quasi-static).
"""


HF_DIPOLE_SPHERE = r"""
# Case Study #6 (2023-03) -- HF Kelvin Transform Validity via Hertz Dipole

Folder: S:/COMSOL/88_ケルビン変換/2023_03_18_高周波のケルビン変換_ダイポール球/

This study answers the question:
   "At what frequency does the Kelvin transformation stop being valid?"

Kelvin is exact for the LAPLACE equation (Coulomb gauge, no time
dependence).  For the HELMHOLTZ equation (full-wave) it is exact only
in the quasi-static limit (k*a -> 0 where k = 2*pi/lambda and a is the
Kelvin ball radius).

The lab's test problem: a Hertz dipole (analytical) at the origin
radiates known E_r, E_theta, E_phi.  Put a Kelvin-image ball at
(0,0,3a) with the inner-air boundary BC ``E0 = genext1(E)``.  Compare:
  - the Kelvin-FEM solution against the analytical dipole field at
    multiple radii
  - and at multiple frequencies (10 MHz, 100 MHz)
  - and against a PML reference

If Kelvin agrees with analytical at 10 MHz but disagrees at 100 MHz
(for some fixed a), we have measured the breakdown point.

## File layout

The folder has THREE subfolders, one per coordinate framing:

| Subfolder | Geometry kind                  | What it tests                           |
|-----------|--------------------------------|----------------------------------------|
| 2D/       | full 2D xy plane               | 2D dipole, easy to debug, fast        |
| axis/     | axisymmetric (r, z)            | 3D dipole reduced to (r,z)             |
| 3D/       | full 3D Cartesian              | the production setup                   |

Lead documents:
| File                                  | Role                              |
|---------------------------------------|-----------------------------------|
| symbolic.m                            | MATLAB symbolic Hertz dipole      |
| ケルビン変換の式.docx                 | Kelvin transform derivation memo  |
| ヘルツダイポールの式.docx             | Hertz dipole derivation memo      |
| 球のケルビン変換.cst                  | CST reference simulation          |
| 球のケルビン変換.pptx                 | summary slides                    |
| 大地の電気定数.png                    | earth-electrical-constants table  |

The MATLAB symbolic.m derives Er, Et, Ep for the Hertz dipole using
the spherical-Hankel-function representation
  A_r = sqrt(pi*k*r/2) * H_(n+1/2)^(2)(k*r) * P_n(cos theta)
which is the standard radiating multipole expansion (n=1 for dipole).
The exterior electromagnetic field is then derived as
  E_r = (1/(j*omega*eps)) * (d^2 A_r / dr^2 + k^2 A_r)
  E_t = -1/(r sin theta) dF_r/dphi + 1/(j omega eps r) d^2 A_r / drdt
  E_p =  1/r dF_r/dt + 1/(j omega eps r) d^2 A_r / drdp

(Standard Stratton-Chu / Harrington TM-mode expansion.)

## 3D case07 family -- the production runs

The 3D subfolder has cases numbered 01-07 plus offset variants.  The
most refined is case07_a_40_offset_20 (also exposed at the TOP of the
88_ケルビン変換 folder as case07_a_40_offset_20.mph, ~11 MB).

case07_a_40 parameters (from 3D/case07_a_40.m):
```matlab
a = 4.0;          % Kelvin ball radius in meters
z0 = -2.0;
z1 = 3.0;
x0 = 5.0;
x1 = 6.0;
y0 = -3.0;
y1 = 3.0;
r1 = 8.0;
f = 100;          % frequency in MHz
xo = 0.0;
zo = 0.0;
mui = 0.2;
ep = 100 - 1000j; % complex permittivity (lossy)
```

Twelve spheres (sph1..sph12) build the multi-layer Kelvin geometry:
  sph1  inner air ball, radius a
  sph2  outer air ball, radius a, center (0, 0, 3a)
  sph3  outer image dipole, radius a^2/r1  (Kelvin-mapped dipole size)
  sph4  near-zone correction ball, radius |a^2/(2*z0)|
  sph5  inner artifact subtract, radius 0.25 m
  sph6, sph7 ... sph12  additional Kelvin-image balls for each
                        observation point (x0, x1, y0, y1, z0, z1)
                        and offset positions

The reason for the 12-sphere construction: the lab wanted to check the
Kelvin solution at SEVERAL observation points simultaneously, so each
observation point gets its own Kelvin-image sphere whose center is the
Kelvin image of that point.

```matlab
sph9 .set('pos',  {'a^2/2/x0' '0' '3.0*a'})    % image of point (x0,0,3a)
sph9 .set('r',    'abs(a^2/2/x0)')              % image-ball radius
sph10.set('pos',  {'a^2/2/x1' '0' '3.0*a'})    % image of (x1,0,3a)
sph10.set('r',    'abs(a^2/2/x1)')
% ... similarly for y0, y1, z0, z1
```

## The genext pattern (2023 version)

```matlab
% From 3D/case01.m, lines 284-315
model.cpl.create('genext1', 'GeneralExtrusion');
model.cpl.create('genext2', 'GeneralExtrusion');
genext1.selection = [e1, e2]    % e1, e2 = boundary face index arrays
genext2.selection = [e1, e2]
genext1.dstmap = {'x', 'y', 'z+3*a'}
genext2.dstmap = {'x', 'y', 'z-3*a'}
genext1.exttol = 0.1            % extrapolation tolerance
genext2.exttol = 0.1

% Boundary condition on each face set
efb1.E0 = analytical Hertz dipole at origin   % the inner source
efb2.E0 = {genext1(Ex), genext1(Ey), genext1(Ez)}   % outer ball gets
                                                     % Kelvin image
efb3.E0 = {genext2(Ex), genext2(Ey), genext2(Ez)}   % inner ball gets
                                                     % outer image
```

The ``exttol = 0.1`` is a 10% tolerance for the genext extrapolation
lookup -- needed because the inner and outer meshes are not topologically
identical (each is auto-meshed by COMSOL).  Compare with WPT_EARTH where
``method = "closest"`` was used; here exttol does the same job differently.

## Animation outputs

Each subfolder has case0X.py which uses PIL to stitch 360 PNG frames
(one per degree of the time phase) into a GIF:

```python
for nt in range(0,360,1):
    FileName = './{:s}/{:03d}.png'.format(FolderName, nt)
    images.append(Image.open(FileName))
images[0].save('{:s}.gif'.format(FolderName),
               save_all=True, append_images=images[1:],
               duration=10, loop=0)
```

The resulting GIFs (case01.gif, ..., case07.gif) are 30-90 MB each and
show the propagating dipole field in the Kelvin domain.  They are the
visual evidence that the Kelvin BC is correctly absorbing the outgoing
wave at the inner sphere boundary.

## Cross-check with PML (case05_PML, case07_PML_*)

The folder also has PML-reference simulations:
  case05_PML.m, case05_PML.mph                  PML cross-check at 100 MHz
  case07_PML_10m.mph, case07_PML_12m.mph        larger PML domain
  case07_PML_inner_25.mph, case07_PML_inner_35.mph  varying PML inner radius

The lab compared Kelvin vs PML at the same frequency / same observation
point.  RESULT (from inspection of the case0N.gif files): at 10 MHz both
methods agree well; at 100 MHz with a = 4 m, both methods have visible
distortions near the inner-outer interface, suggesting that this is the
borderline frequency for a = 4 m (k*a = (2*pi*100e6 / 3e8) * 4 ~ 8.4 --
clearly NOT quasi-static).

LESSON: the rule of thumb "Kelvin works when k*a << 1" was empirically
validated by this study.  For Radia / NGSolve work in MQS / Darwin regime,
this is always satisfied (omega -> 0 limit) so Kelvin is always exact.

## The conference figure

``axis/conference.png`` (346 kB) is the conference-talk figure
summarizing the Kelvin vs PML comparison for the axisymmetric case at
10 MHz and 100 MHz.  Lab published these results in an internal seminar.

## Status

VALIDATED 2023-03.  Lab reference for "when does Kelvin break for
full-wave problems."  Answer: when k*a is no longer << 1.  For Radia
MQS / Darwin work (CLAUDE.md policy: Laplace kernel only, no Helmholtz),
this is never a concern.
"""


FORCED_PERIODIC_BC_SCRIPT = r"""
# Case Study #3 (2021-09) -- Forced Periodic-BC MATLAB+Java Workaround

Folder: S:/COMSOL/88_ケルビン変換/2021_09_01_各面に対して強制的に周期境界条件を課すスクリプト/

This is the SHORTEST but most reusable case study.  It is a generic
MATLAB script that takes any COMSOL .mph file, enumerates ALL its
boundary faces, automatically pairs them up by centroid offset, and
batch-creates a PeriodicCondition on every pair.  Plus matching
CopyFace mesh operations to ensure the inner-outer meshes are
conforming.

The same automation pattern applies directly to Kelvin face pairing:
substitute the periodic BC with a GeneralExtrusion + outer-face BC.

## Lead files

| File                                  | Size       | Role                                  |
|---------------------------------------|------------|---------------------------------------|
| Set_periodic_COMSOL.m                 |   1.9 kB   | the workhorse MATLAB driver           |
| boundary.java                         | 383 kB     | auto-generated COMSOL Java script     |
| Boundary.mat                          |  68 kB     | cached boundary metadata              |
| in_out.mat                            |   1.8 kB   | cached inner/outer face pair list     |
| case1.m                               |  37 kB     | the base-case rebuild .m              |
| 2021_09_01_粗.mph                     |  12 MB     | coarse-mesh input to the script       |
| 2021_09_01_粗_peri.mph                |  14 MB     | coarse output, after script           |
| 2021_09_01_密.mph                     |  42 MB     | fine-mesh input                       |
| 2021_09_01_密_peri.mph                |  55 MB     | fine output                           |
| sphere.sat                            |   2.8 MB   | the geometric input                   |
| 周期境界の練習.mph                    |  47 MB     | preliminary practice file             |

## The script (Set_periodic_COMSOL.m) -- annotated

```matlab
clear all; close all;

% Load the COMSOL .mph that already has its geometry / meshes / physics
% set up.  This script ADDS the periodic BC + CopyFace mesh ops.
model = mphload('2021_09_01_粗.mph');

% ---- Stage 1: enumerate ALL boundaries + their coordinates + centroids ----
FileName = 'Boundary.mat';
delete(FileName);
if exist(FileName,'file')~=2
    tic;
    N = model.geom('geom1').getNBoundaries;
    for n = [1:N]
        R = mphgetcoords(model,'geom1','boundary',n);  % vertex coords
        Boundary(n).R = R;
        Boundary(n).G = mean(R,2);                     % centroid
        Boundary(n).inside = 0;
    end
    save(FileName,'Boundary','N');
    toc
else
    load(FileName);
end

% ---- Stage 2: classify each boundary as inside / outside ----
% A face is "inside" if all its vertices lie on the inner sphere |x|=0.3
FileName = 'in_out.mat';
delete(FileName);
if exist(FileName,'file')~=2
    for n = [1:N]
        if size(Boundary(n).R,2)==3       % triangular face has 3 vertices
            r1 = norm(Boundary(n).R(:,1));
            r2 = norm(Boundary(n).R(:,2));
            r3 = norm(Boundary(n).R(:,3));
            r0 = mean([r1,r2,r3]);
            if abs(r0 - 0.3) < 1e-4       % all vertices at |x| = 0.3
                Boundary(n).inside = 1;
            end
        end
    end
    % ---- Stage 3: pair each inside face with its outside image ----
    in_out = [];
    up = 0.63;                           % = 2.1 * 0.3 (Kelvin offset)
    for n1 = find([Boundary.inside]==1)
        Boundary(n1).outside = 0;
        for n2 = find([Boundary.inside]==0)
            % outside-face centroid should equal (inside centroid) + (0,0,up)
            if (abs(Boundary(n2).G(1)-Boundary(n1).G(1))
              + abs(Boundary(n2).G(2)-Boundary(n1).G(2))
              + abs(Boundary(n2).G(3)-Boundary(n1).G(3)-up)) < 1e-4
                Boundary(n1).outside = n2;
                in_out = [in_out; [n1,n2]];
            end
        end
    end
    save(FileName,'in_out');
else
    load(FileName);
end

% ---- Stage 4: BUILD the FreeTri mesh on inner faces ----
model.component('comp1').mesh('mesh1').create('ftri1', 'FreeTri');
model.component('comp1').mesh('mesh1').feature('ftri1').selection
     ().set([in_out(:,1)]);

% ---- Stage 5: for EACH pair, create one PeriodicCondition + CopyFace ----
for n = [1:size(in_out,1)]
    n1 = in_out(n,1);    % inside face
    n2 = in_out(n,2);    % outside (paired) face
    model.component('comp1').physics('emw').create(sprintf('pc%d',n),
        'PeriodicCondition', 2);
    model.component('comp1').physics('emw').feature(sprintf('pc%d',n))
        .selection.set([n1 n2]);
    model.component('comp1').mesh('mesh1').create(sprintf('cpf%d',n),
        'CopyFace');
    model.component('comp1').mesh('mesh1').feature(sprintf('cpf%d',n))
        .selection('source').set(n1);
    model.component('comp1').mesh('mesh1').feature(sprintf('cpf%d',n))
        .selection('destination').set(n2);
end

% ---- Stage 6: finalize: free-tet mesh on the rest, save ----
model.component('comp1').mesh('mesh1').create('ftet1', 'FreeTet');
mphsave(model,'2021_09_01_粗_peri.mph');
```

## Why this script exists

COMSOL's built-in Identity Pair / Periodic Condition GUI workflow
EXPECTS you to manually select the inner face and the outer face for
each pair.  For a sphere meshed into 50+ triangles per side, this is
50+ pair selections -- error-prone and tedious.  The script automates
the pairing by centroid comparison.

This same logic is what Radia's C++ Cubit plugin does internally for
Kelvin face identification: the plugin's NetgenCurver-based code
extracts each inner Kelvin face, computes its Kelvin image, and writes
the resulting outer-face pair into the .vol file's Identifications
section.  See ``src/cubit_plugin/ExportNetgenCommand.cpp`` for the C++
implementation.

## Why "Periodic" and not "Kelvin"

The script as written uses ``PeriodicCondition``, not a Kelvin BC.
This is because the 2021_09 study was a STEPPING STONE: the lab first
proved they could automate the PAIRING (study #3), then in study #4
substituted the periodic BC for genext-based Kelvin BC (which is
trivial once the pair list is known).

For Kelvin face pairing the only change needed:
  - replace ``PeriodicCondition`` with a ``GeneralExtrusion`` cpl
    pointing from n1 to n2 with dstmap = (x, y, z + up)
  - replace the BC with ``E0 = genext_n(Ex, Ey, Ez)`` on face n2

## Status

VALIDATED 2021-09 for periodic BC.  Re-used (with the substitution
above) in studies #4, #5, #6, #7.  This script is THE BLUEPRINT for
boundary-pair automation in any FEM tool.  Recommended reading for
anyone implementing similar functionality in NGSolve / Radia.

## Inspiration for Radia automation

The same .mat caching idea (Boundary.mat, in_out.mat) is what Radia
should adopt for repeated Kelvin solves:
  - cache the inner-outer face pair list on disk (or in the .vol's
    Identifications section)
  - on subsequent loads, skip re-pairing
For meshes with 10000+ faces the pairing step would otherwise take
seconds; with the cache it is instant.
"""


ECT_EDDY_CURRENT_DISTRIBUTION = r"""
# Case Study #5 (2021-11) -- ECT Eddy-Current Distribution Post-Processing

Folder: S:/COMSOL/88_ケルビン変換/ (top-level) + ECT_01/ subfolder of #2

Lead files:
  2021_11_06_渦電流探傷_渦電流分布用.mph   13.3 MB   the eddy-current
                                                     distribution analysis run

  Inside ECT_01/ (case study #2 working folder):
  2021_11_06_渦電流探傷_信号分布用.mph    16.5 GB   the SIGNAL distribution
                                                     run (huge -- contains the
                                                     full transient sweep)
  渦電流分布.py                              1.8 kB   matplotlib post script
  contour表示.py                             2.3 kB   3D contour post script
  B..csv, B.txt                              ~30 kB   field-extraction data
  LR.txt                                     71 kB    inductance / resistance
                                                     sweep results
  case1.vtu                                  67 MB    VTK output for 3D viz
  3d-example_result1.pdf .. result3.pdf    23-302 kB  paper-grade figures
  figure.png, 磁場の変化.png                500 kB   summary plots

## Application

ECT = Eddy Current Testing.  An AC excitation coil is brought close to
a metallic workpiece.  The induced eddy currents in the workpiece
create a secondary magnetic field that perturbs the coil's
inductance.  By measuring R + jwL at the coil terminals as a function
of position, defects (cracks, voids) can be located.

The Kelvin transform here serves the same purpose as in study #2:
provides an exact open-boundary condition without truncating air
artificially.  Study #5 adds detailed POST-PROCESSING of the eddy
current density distribution INSIDE the workpiece, which is the
quantity actually scored against analytical / experimental references.

## Post-processing recipe (the .py scripts)

```python
# 渦電流分布.py -- 2D line plot of Im(Jx) vs y at fixed x, z
# (from a COMSOL "Cut Line" data export)

import matplotlib.pyplot as plt
from numpy import *
import os, sys, re

# Read COMSOL-exported text dump
y1, J1 = [], []
fid = open('case1', 'r')
strLines = fid.readlines()
for n in range(len(strLines)):
    sline = re.split(r"\s+|\n$", strLines[n])
    if sline[0][0] != '%':            # skip COMSOL header lines
        y1.append(float(sline[0]))
        J1.append(float(sline[1]))

# ... same for 'case2'

# IEEE-style 2-column figure (3 inch wide, 2 inch tall, 300 dpi)
plt.figure(figsize=(3, 2), dpi=300)
ax1 = plt.axes([0.16, 0.19, 0.80, 0.76])
ax1.tick_params(direction="in")
ax1.plot(y2, array(J2)/1000, 'k-', linewidth=0.7, label='a=8mm')
ax1.plot(y1, array(J1)/1000, 'r--', linewidth=0.7, label='a=6mm')
plt.savefig('3d-example_result2.pdf')
```

```python
# contour表示.py -- 3D surface plot of Bi(x, y) at fixed z
# (from a COMSOL "Cut Plane" 2D grid export)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

aa, xx, yy, Br, Bi = ...  # parsed from B.txt

fig = plt.figure(figsize=(4, 3.3), dpi=600)
ax1 = plt.axes([0.14, 0.03, 0.82, 0.98], projection="3d")
surf = ax1.plot_surface(xx, yy, Bi, cmap=cm.coolwarm,
                         linewidth=0, antialiased=False)
ax1.view_init(elev=40, azim=45)
ax1.set_xlabel('$\\itx$ (mm)', fontname='times new roman', fontsize=15)
ax1.set_ylabel('$\\ity$ (mm)', fontname='times new roman', fontsize=15)
ax1.set_zlabel('$\\itB_z (\\mu \\rm{T})$', fontname='times new roman', fontsize=15)
plt.savefig('3d-example_result3.pdf')
```

These scripts are TEMPLATES for paper-grade ECT post-processing:
  - IEEE-style 2-column figure dimensions (3 inch)
  - Times New Roman fonts (matches IEEE Author Center spec)
  - axis ticks pointing INWARD (lab-graph-style policy)
  - dpi >= 300 for raster / pdf for vector
  - all numeric labels with explicit units in italic

See ``radia_mcp.electromagnet.bibliography_index`` for the equivalent
Radia / NGSolve post recipe (which produces .msh v4.1 for GMSH instead
of pdf, but the same figure conventions apply for the final paper).

## Status

VALIDATED 2021-11 for ECT post-processing.  The two .py scripts are
the lab's canonical templates for ECT figure generation.  The
hundreds-of-MB .vtu and .pdf outputs document the validation runs.

## Cross-reference to Radia

Radia panel ``radia_ih.py`` (Layer 3) running ``calc_fem_kelvin.py``
(Layer 4) produces equivalent J distributions for the IH workpiece
problem.  The post-processing template in this study applies directly:
just substitute the COMSOL-exported text file with a GMSH .msh v4.1
exported by ``radia.gmsh_post_export.GmshPostExport``.
"""


LESSONS_LEARNED = r"""
# Cross-Cutting Lessons (from all 7 lab case studies)

These are the gotchas that EACH study hit, with the workaround that
worked.  They are arranged in priority order -- the most important
ones first.

## L1. Always use GeneralExtrusion with dstmap, NEVER LinearExtrusion

In 2020 (case #1) the lab used four separate LinearExtrusion couplings,
one per face quadrant, each pairing vertices manually.  By 2021-06
this had been replaced everywhere with two GeneralExtrusion + ONE
dstmap expression like ``{x, y, z + 3*a}``.

GeneralExtrusion + dstmap is:
  - 4x fewer cpl operations
  - portable across mesh refinements (vertex IDs don't matter)
  - readable (the dstmap is the actual map expression)

LESSON: never use LinearExtrusion for Kelvin BC.  Use GeneralExtrusion.

## L2. The Kelvin ball must STRICTLY contain all conductive geometry

From case #2 (ECT, 2021-06):
  "大地 Brickの、xy方向に110が最大で、それ以上がKelvin変換出来なかった"

When the earth slab extended beyond the Kelvin ball boundary, the
inversion mapped part of the earth INTO the outer ball geometry, which
broke the BC.  Solution: scale the Kelvin ball radius up so the earth
is strictly inside.

LESSON: the Kelvin ball is your "air bubble" -- nothing conductive /
magnetic can poke through.  If the geometry is large, the ball must be
larger.  Validate by listing every domain's bounding box and confirming
all are inside |x| < a.

## L3. dstmap z-sign is easy to get wrong (the "NG" vs "OK" mph)

In case #4 (WPT, 2021-10):
  NG_Kelvin_RF_3D_genext_PML_earth.mph (Aug 22 2021) -- coordinate
   sign issue in genext, field did not decay smoothly at inner-outer
   interface
  OK_Kelvin_RF_3D_genext_PML_earth.mph (Oct 1 2021) -- fixed by
   creating SEPARATE forward (z + up) AND inverse (z - up) genext
   operators

LESSON: each face pair needs TWO genext operators (forward and inverse).
Verify by symmetry: if A's image is B, then B's image (under the
inverse map) is A.  Check both directions are wired up.

## L4. method="closest" silently masks mesh-mismatch bugs

In case #4 the lab set ``method = "closest"`` on two of the eight
genext operators.  This works around topological mismatches between
inner and outer meshes, but masks them.

LESSON: use ``method = "closest"`` only as a localized last-resort
fix on the two faces that genuinely cannot be matched.  Document WHY
each "closest" is needed.  In Radia's NGSolve workflow, the equivalent
is the C++ Cubit plugin's tolerance-based vertex matching -- the
"all-or-nothing" mode is preferred; mesh-by-mesh tolerance overrides
should be flagged in the .jou comments.

## L5. Identify boundary face pairs automatically (case #3 pattern)

The forced-periodic-BC MATLAB script (case #3) demonstrates that
boundary pairing by CENTROID + COORDINATE OFFSET is robust enough for
production use.  Cache the result in a .mat (or in the .vol's
Identifications section) so subsequent solves skip the pairing step.

LESSON: never hand-pair faces in the GUI.  Always script the pairing.

## L6. Mesh the inner and outer side identically (CopyFace)

The forced-periodic-BC script uses CopyFace on every face pair to
enforce IDENTICAL inner / outer surface meshes.  Without CopyFace, the
default FreeTri produces non-conforming meshes on the two sides, which
forces the genext to fall back to interpolation -- inaccurate.

LESSON: always pair PeriodicCondition / genext with CopyFace.  In
NGSolve, the equivalent is the CopyFaceMesh option on Refine() or the
C++ Cubit plugin's vertex-matched mesh export.

## L7. Pick the offset = 3 * a (lab convention 2021+)

Case #1 used up = 2.1 * a (barely > 2a, balls touch with 0.1a gap).
Case #4 inherited this.  Case #6 (2023) firmly settled on up = 3 * a.

LESSON: use offset = 3a always.  It gives a 1a gap between balls,
which is visually unambiguous and leaves room for an air domain in
between.  2.1*a "works" but the gap is tiny and easy to misread.

## L8. Kelvin is exact ONLY in the Laplace / Coulomb-gauge limit

Case #6 (HF Hertz dipole) measured that Kelvin breaks down at high
frequency when k*a is no longer << 1 (k = wave number, a = ball
radius).  At k*a ~ 8 the Kelvin BC visibly distorts the field at the
inner-outer interface.

LESSON: for Helmholtz / full-wave problems, Kelvin works only when
k*a << 1 (quasi-static).  Beyond that, use PML or BEM.  For Radia
MQS / Darwin regime (CLAUDE.md: Laplace kernel only), this is never a
concern.

## L9. The Kelvin SAT trick (pre-bake the inversion)

Case #4 imported a pre-built Kelvin.sat (1.7 MB) that already contains
the Kelvin-inverted geometry of the wireless coil.  This was generated
once in Cubit by Kelvin-earth_300R.py (vertex-by-vertex inversion via
the patran-export-rewrite-reimport trick), then re-used in every
COMSOL rebuild.

LESSON: do the geometric inversion ONCE in Cubit, save as a SAT or
STEP, then import.  Re-doing the inversion per-COMSOL-rebuild is
slow and bug-prone.

## L10. Combine Kelvin + PML only when you must (case #4 unique)

The lab's WPT@12MHz study (case #4) is the unique case that combined
Kelvin AND PML in one mesh.  The reason: at 12 MHz the geometry is
borderline quasi-static, with some leakage propagating wave.  PML
catches the wave residue; Kelvin handles the quasi-static far-field.

LESSON: in plain MQS / Darwin work (Radia / IH), use Kelvin ONLY.
Adding PML offers no benefit and complicates the solver setup.  Only
combine the two when measured Kelvin leakage at the operating
frequency exceeds engineering tolerance.

## L11. Validate with an analytical reference (case #6 pattern)

Case #6 chose the Hertz dipole as the analytical reference because its
exterior field has a closed-form Hankel-function expression.  Comparing
the Kelvin-FEM field at multiple radii against this analytical formula
is THE definitive test of whether the Kelvin BC is correctly imposed.

LESSON: every new Kelvin setup should include an analytical-reference
comparison.  Cross-link: ``radia_mcp.radia_ngsolve.analytical_formulas``
has the Radia closed-form library; pick a relevant problem and validate
against it.
"""


ROBIN_LIKE_PATTERN = r"""
# The Robin-like Pattern of Kelvin BC in COMSOL

A subtle but important point that the case studies expose:

COMSOL's ``E0 = genext1(E)`` on the OUTER face is NOT a classical
Robin BC.  It is a STRONG (Dirichlet-like) CONTINUITY condition that
matches the OUTER-side field value to the INNER-side field's Kelvin
image.

This is implemented via "Electric Field" or "Magnetic Potential"
boundary conditions (efb1, efb2, ... in the .m files), which in COMSOL's
weak form translate to a tangential-field Dirichlet:
  n x E_outer  =  n x genext1(E_inner)

For the magnetic vector potential A:
  n x A_outer  =  n x genext1(A_inner)

This is exactly the Kelvin transformation written in WEAK FORM at the
inner-outer interface.

## Why this is "Robin-like" in spirit

If you think of the OUTER FEM domain as a "shell" that represents the
exterior of the inner ball via the Kelvin transform, then the inner
ball's solution is the SOURCE for the outer ball's BC.  The outer
ball's solution feeds back as the inner ball's far-field BC.  This is
EXACTLY the Robin coupling pattern between two domains -- with the
"reflection coefficient" being the Kelvin geometric weight (which
COMSOL absorbs into the genext map).

## Comparison: NGSolve weak-form Kelvin (Radia panels)

In NGSolve (radia ``calc_fem_kelvin.py``), the equivalent is:

```python
# Inner FEM space + outer FEM space (one ball each)
fes = HCurl(mesh, order=p)

# Identify inner Kelvin face with outer Kelvin face
# (via C++ Cubit plugin's auto-pairing during .vol export)
# The .vol file's Identifications section lists pair (face_inner, face_outer)

# Set up the bilinear form with the SAME weak Kelvin BC
a += curl(u) * curl(v) * weight_inner * dx("inner")
a += curl(u) * curl(v) * weight_outer * dx("outer")
# weight_outer = (a / |y - c|)^4 in 3D for HCurl (degree-1 form)
#              (verify via differential_forms_mathematica_recipes('kelvin'))
```

The Identifications section in the .vol file enforces n x A_outer =
n x A_inner at the paired faces, which IS the Kelvin BC.

## Practical implication for the Sugahara lab

A side-by-side Kelvin solve in COMSOL and in NGSolve / Radia on the
SAME geometry should give IDENTICAL inner-domain solutions, modulo
mesh resolution.  This is the SINGLE BEST cross-validation we have
between the two FEM tools, and the case studies above (1-7) provide
the COMSOL reference solutions.

To set up a comparison:
  1. Export the inner-ball mesh from COMSOL to Netgen .vol via the
     comsol2vol script (if available) or via STEP -> Netgen.
  2. Pair the inner-outer faces using the same centroid-offset logic
     as Set_periodic_COMSOL.m.
  3. Run NGSolve Kelvin solve with the matched Identifications.
  4. Extract the inner-domain field and compare with COMSOL's at
     the same observation points.
"""


CROSS_REFERENCE = r"""
# Cross-Reference: How Each Case Study Maps to the Radia Pipeline

## Case #1 (basic practice, 2020-06)
   -> Radia analog: any unit test that solves a sphere-in-uniform-field
                    with Kelvin (no equivalent exists yet; would be a
                    good addition to tests/panels/test_kelvin_*).

## Case #2 (ECT + earth, 2021-06)
   -> Radia analog: ``calc_fem_kelvin.py`` with earth-slab geometry +
                    induction current coil source.  The lab's COMSOL
                    .mph files (Kelvin_AC_coil.mph, COIL_case1_30.mph,
                    Case3_ゲージ固定無_無限遠.mph) are direct
                    validation targets.

## Case #3 (forced periodic BC script, 2021-09)
   -> Radia analog: the C++ Cubit plugin's auto-pairing of Kelvin faces
                    during ``radia_export netgen`` (see
                    ``src/cubit_plugin/ExportNetgenCommand.cpp`` and the
                    Identifications block written into the .vol).
                    The MATLAB script is the algorithmic blueprint.

## Case #4 (WPT@12MHz + earth + PML + Kelvin, 2021-10)
   -> Radia analog: ``radia.peec_topology.PEECCircuitSolver`` for the
                    coil + ``calc_fem_kelvin.py`` for the workpiece +
                    earth.  Combining PML and Kelvin is OUTSIDE Radia
                    scope (CLAUDE.md: Laplace kernel only).  At 12 MHz
                    with iron, use ``calc_inductance.py --coil-solver peec``
                    instead of trying to add PML.

## Case #5 (ECT eddy-current distribution post, 2021-11)
   -> Radia analog: post-process the FEM solve via
                    ``radia.gmsh_post_export.GmshPostExport`` to GMSH
                    .msh v4.1, then use the lab's matplotlib templates
                    (渦電流分布.py, contour表示.py) for paper figures.

## Case #6 (HF Hertz dipole validity check, 2023-03)
   -> Radia analog: NONE.  Radia is MQS / Darwin only (Laplace kernel,
                    no Helmholtz).  This case study is a one-time
                    validity check that confirms Kelvin is exact in the
                    Radia regime (k*a -> 0).  It does NOT need to be
                    re-run; the lesson is "Kelvin is always exact for
                    Radia."

## Case #7 (worked example case07_a_40_offset_20.mph, 2023)
   -> Radia analog: a parametric sweep tutorial.  Mirror this in
                    ``examples/induction_heating/`` with a similar
                    a=40 sweep of offset positions.

## Summary table

| Lab study | Radia component         | Status                               |
|-----------|-------------------------|--------------------------------------|
| #1 basic  | (none yet)              | Add as a test fixture                |
| #2 ECT    | calc_fem_kelvin.py      | Production validated 2021-06         |
| #3 script | Cubit plugin C++        | Production validated 2021-09         |
| #4 WPT    | calc_inductance.py PEEC | Production; skip PML side            |
| #5 ECT post | GmshPostExport        | Production; reuse matplotlib templates|
| #6 HF dipole | (none)                | One-time validity check; out of scope|
| #7 worked | example mph             | Reference; mirror in examples/       |
"""


LAB_FILES = r"""
# Lab File Catalog -- S:/COMSOL/88_ケルビン変換/

Verbatim NAS paths (Japanese folder names preserved).  All sizes / dates
captured 2026-05-25.

## Top-level files

```
S:/COMSOL/88_ケルビン変換/
    2021_11_06_渦電流探傷_渦電流分布用.mph            13.3 MB   case #5 short
    case07_a_40_offset_20.mph                          11.4 MB   case #7
```

## Case #1 (2020-06): basic Kelvin practice

```
S:/COMSOL/88_ケルビン変換/2020_06_12_Kelvin変換の練習/
    Kelvin_2D.mph                  2.0 MB     2D
    Kelvin_3D.mph                  7.0 MB     3D full
    Kelvin_3D_4Q.mph               5.9 MB     3D quarter-sym
    Kelvin_3D.m                    79 kB      MATLAB rebuild
    練習.mph                       3.7 MB     scratch
    周期境界練習.mph               1.2 MB     periodic BC practice
    case1.vtu                      419 kB     VTK output
    case2.vtu                      5.4 MB     VTK output
    isosurface.vtu                 34 MB      isosurface
    streamline.vtu                 16 MB      streamlines
    円周上の値の引継ぎ.png         1.1 MB     notes
    202142_0-7-13.glance           15 MB      glance recording
    202141_21-13-53.glance         1.4 kB     glance recording
```

## Case #2 (2021-06): ECT + earth

```
S:/COMSOL/88_ケルビン変換/2021_06_18_大地を模擬したケルビン変換@ECT/
    main.m                                         18.6 kB
    Kelvin_AC_coil.mph                             41.4 MB    production
    COIL_case1_30.mph                              54.7 MB
    case2_30.mph                                   22.9 MB
    Case3_ゲージ固定無_無限遠.mph                  13.0 MB    Kelvin version
    Kelvin_2D_mag_練習.mph                          2.0 MB
    Kelvin_2D_wav_練習.mph                          0.3 MB
    2021_06_16_ケルビン変換.xlsx                   16 kB      pair mapping

    ECT_01/                                                   subfolder
        2021_11_06_渦電流探傷_信号分布用.mph       16.5 GB    full transient
        2021_11_06_渦電流探傷_渦電流分布用.mph     93.2 MB
        渦電流分布.py                              1.8 kB     post script
        contour表示.py                             2.3 kB     post script
        B.txt, B..csv                              ~30 kB
        LR.txt                                     71 kB      L/R sweep
        case1.vtu                                  67 MB
        3d-example_result1.pdf                     302 kB
        3d-example_result2.pdf                      21 kB
        3d-example_result3.pdf                      24 kB
        figure.png, 磁場の変化.png                 ~500 kB

    コイルのみ/
        case1_matlab.m                             ~30 kB

    練習_Case1/                                     ゲージ固定有
        Case1_ゲージ固定有.mph                     16.4 MB
        case1_100.txt, case1_30.txt, ...           ~kB

    練習_Case2/
        Case2_ゲージ固定無.mph                     15.6 MB
        2021_07_10-4_case2_30.txt                  ~kB
        100/                                       sweep folder

    練習_Case3/
        Case3_ゲージ固定無_無限遠.mph              13.0 MB   ★ Kelvin
        coreform/
            2021_07_24_描写用.txt                  2 kB      view setup notes
            ※問題ありで休止中.txt                 ~bytes    KEY LESSON
                "大地 Brickの、xy方向に110が最大で、
                 それ以上がKelvin変換出来なかった"
```

## Case #3 (2021-09): forced-periodic-BC script

```
S:/COMSOL/88_ケルビン変換/2021_09_01_各面に対して強制的に周期境界条件を課すスクリプト/
    Set_periodic_COMSOL.m                          1.9 kB     ★ MATLAB driver
    boundary.java                                 383 kB     auto-generated
    Boundary.mat                                   68 kB     cached metadata
    in_out.mat                                      1.8 kB    pair cache
    case1.m                                        37 kB
    sphere.sat                                      2.8 MB    geometry input
    2021_09_01_粗.mph                             12.1 MB    coarse input
    2021_09_01_粗_peri.mph                        14.4 MB    coarse output
    2021_09_01_密.mph                             41.9 MB    fine input
    2021_09_01_密_peri.mph                        55.0 MB    fine output
    周期境界の練習.mph                            47.2 MB    practice
```

## Case #4 (2021-10): WPT + earth + PML + Kelvin

```
S:/COMSOL/88_ケルビン変換/2021_10_22_大地を模擬したケルビン変換@WPT/
    Kelvin.sat                                       1.7 MB    pre-built
    Kelvin_RF_3D_genext_PML_earth.java              71 kB     ★ COMSOL Java
    NG_Kelvin_RF_3D_genext_PML_earth.mph            18 MB     "not good"
    OK_Kelvin_RF_3D_genext_PML_earth.mph            33 MB     "OK"
    OK_Kelvin_RF_3D_genext_PML_earth_wireless.mph  2.4 GB     ★★ production
    Kelvin-earth_300R.py                            1.4 kB    Cubit Python
    Comsol用_Kelvin対面.xlsx                       10.6 kB   pair mapping
    airplane_antenna_crosstalk.pptx                 2.6 MB    application
    動画作成.py                                     394 B     GIF builder
    電界分布の式.pptx                              45 kB     theory

    00_coreform/                                              Cubit work
        ###########下の大地と空気.txt                          notes
        2021_08_25_Full-model.cub5                            ~11 MB
        2021_08_28_Kelvin-earth-air(tris).py                  2.3 kB
        2021_08_28_Kelvin-earth-air(quads).py                 2.3 kB
        Kelvin-earth_R300.cub5                                1.9 MB
        Kelvin_air-earth_300R.py                              2.5 kB
        sphere.py                                             1.8 kB
        資料表示用Kelvin-PML.py                              1.4 kB
        ... 50+ more .jou and .cub5 files

    NG_メッシュ読込_wireless_Kelvin_PML_銅損有_大地/          NG variant
    NG_メッシュ読込_wireless_Kelvin_no_hole_08_16/            NG variant
    SAT読込_wireless_Kelvin_PML_銅損有_大地_NG/               SAT NG
    SAT読込_wireless_Kelvin_PML_銅損有_大地_σ=1e7/            sigma sweep
    SAT読込_wireless_Kelvin_PML_銅損有_空気大地/              air+earth
    SAT読込_wireless_Kelvin_no_hole_銅損有_大地 _σ=1e7/       (no hole)
    SAT読込_wireless_Kelvin_no_hole_銅損有_大地 _σ=1e7
      _インピーダンス境界_NG/                                  impedance BC NG
    SAT読込_wireless_Kelvin_no_hole_銅損有_大地 _σ=1e7
      _インピーダンス境界_retry/                               impedance BC retry
    SAT読込_wireless_Kelvin_no_hole_銅損有_空気大地/
    COMSOL_CAD_wireless_Kelvin_no_hole_銅損有_空気大地/
    wireless_Kelvin_PML_銅損有/
    wireless_Kelvin_no_hole_完全導体/
    wireless_Kelvin_no_hole_銅損有/
    wireless_rad_remesh/
    png/                                                       all PNG outputs
```

## Case #5 (2021-11): ECT eddy current distribution

See ECT_01/ subfolder under case #2 above.

## Case #6 (2023-03): HF Hertz dipole

```
S:/COMSOL/88_ケルビン変換/2023_03_18_高周波のケルビン変換_ダイポール球/
    symbolic.m                                        842 B      Hertz dipole
    ケルビン変換の式.docx                            21 kB      Kelvin memo
    ヘルツダイポールの式.docx                        92 kB      dipole memo
    大地の電気定数.png                                2.0 MB     soil constants
    球のケルビン変換.cst                              94 kB      CST cross-ref
    球のケルビン変換.pptx                            249 kB      slides

    2D/
        case01.m, case01.py, case01.gif (45 MB)                2D case 1
        case02.m, case02.py, case02.gif (45 MB)                2D case 2
        case01/, case02/                                       PNG dirs

    axis/                                                       axisymmetric
        case01.gif (46 MB) ... case07.gif (92 MB)              animations
        case01.py ... case07.py                                GIF builders
        case01_10MHz_Kelvin.m                                  Kelvin 10MHz
        case01_10MHz_PML.m                                     PML 10MHz
        case01_100MHz.m                                        100 MHz
        case04_100MHz.m, case05_100MHz_offset.m
        case06_100MHz_GND.m, case06_a_3.m
        case07_100MHz_soil.m                                   ★ soil@100MHz
        contour_case01.m ... contour_case07.m                  contour scripts
        conference.png                                         summary fig

    3D/
        case01.m, case01.mph (200 MB)
        case01.py, case01/                                     PNGs
        case02.m, case03_a_30.m
        case03_a_40.m, case03_a_40/
        case04_a_40.m
        case05_a_40.m, case05_PML.m
        case06_a_30.mph (58 MB), case06_a_40.m
        case06_a_40_offset_00.mph, ..._10.mph, ..._20.mph
        case07_a_25.m .. case07_a_45.m                         a sweep
        case07_a_40.m, case07_a_40_offset_00.m ..._20.m        ★ production
        case07_a_40_offset_20.mph (89 MB)                     ★★ this is
                                                                also exposed
                                                                at top level
        case07_PML_10m.mph, case07_PML_12m.mph                 PML compare
        case07_PML_inner_25.mph, case07_PML_inner_35.mph       PML inner sweep
        contour_case01.m ... contour_case07.m                  contour scripts
        Figure_offset_*.m, Figure_size_*.m, Figure_scatter_*.m post scripts
        batch.m, batch.log, batch_2.log                        batch driver
```

## Case #7 (worked example)

See ``case07_a_40_offset_20.mph`` at the top level AND in 3D/.  Same
file, just two copies for convenience.
"""


TOPICS = {
    "overview": OVERVIEW,
    "basic_practice": KELVIN_BASIC_PRACTICE,
    "kelvin_basic_practice": KELVIN_BASIC_PRACTICE,
    "ect_earth": ECT_EARTH,
    "wpt_earth": WPT_EARTH,
    "hf_dipole_sphere": HF_DIPOLE_SPHERE,
    "forced_periodic_bc_script": FORCED_PERIODIC_BC_SCRIPT,
    "forced_periodic": FORCED_PERIODIC_BC_SCRIPT,
    "ect_eddy_current_distribution": ECT_EDDY_CURRENT_DISTRIBUTION,
    "ect_post": ECT_EDDY_CURRENT_DISTRIBUTION,
    "lessons_learned": LESSONS_LEARNED,
    "lessons": LESSONS_LEARNED,
    "robin_like_pattern": ROBIN_LIKE_PATTERN,
    "robin": ROBIN_LIKE_PATTERN,
    "cross_reference": CROSS_REFERENCE,
    "cross_ref": CROSS_REFERENCE,
    "lab_files": LAB_FILES,
    "files": LAB_FILES,
}


def get_kelvin_lab_studies(topic: str = "overview") -> str:
    """Dispatch Kelvin-lab-studies knowledge.

    Topics (aliases in parentheses):
        overview                  - catalog + lab convention (DEFAULT)
        basic_practice            - case #1: 2020-06 sphere practice
        ect_earth                 - case #2: 2021-06 ECT + earth
        forced_periodic_bc_script - case #3: 2021-09 MATLAB pairing script
        wpt_earth                 - case #4: 2021-10 WPT + earth + PML
        ect_post                  - case #5: 2021-11 ECT post-processing
        hf_dipole_sphere          - case #6: 2023-03 HF Kelvin validity
        lessons_learned (lessons) - cross-cutting pitfalls (L1-L11)
        robin_like_pattern        - genext-as-Robin physical picture
        cross_reference           - how each maps to Radia pipeline
        lab_files (files)         - NAS path catalog with sizes
        all                       - everything concatenated
    """
    topic = topic.lower().strip()
    if topic == "all":
        return "\n\n".join([
            OVERVIEW,
            KELVIN_BASIC_PRACTICE,
            ECT_EARTH,
            FORCED_PERIODIC_BC_SCRIPT,
            WPT_EARTH,
            ECT_EDDY_CURRENT_DISTRIBUTION,
            HF_DIPOLE_SPHERE,
            LESSONS_LEARNED,
            ROBIN_LIKE_PATTERN,
            CROSS_REFERENCE,
            LAB_FILES,
        ])
    if topic in TOPICS:
        return TOPICS[topic]
    return (
        f"Unknown topic '{topic}'.  Available: "
        + ", ".join(sorted(set(TOPICS.keys()))) + ", all."
    )
