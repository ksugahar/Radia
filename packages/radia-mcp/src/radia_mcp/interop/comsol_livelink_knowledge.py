"""COMSOL LiveLink knowledge for radia_mcp.interop.

Distilled from the Sugahara Lab COMSOL working tree at S:/COMSOL/.
COMSOL Multiphysics exposes the same model object through three
front-end APIs:

  * Java   -- the canonical "save as model.java" output (binary-compatible
              with the GUI; this is the lab's primary scripting path)
  * MATLAB -- "LiveLink for MATLAB" (mli), the same Java API wrapped in
              the MATLAB workspace via mphload/mphsave
  * Python -- MPh (third-party, MIT) over comsolmphserver, the most
              ergonomic option for headless / CI / parameter sweeps

This knowledge module documents all three, plus the lab's own working
examples (model1.java, model2.java, ForJavaScript.java) and the
PowerShell compile/run pipeline that ships with them.

Section map:
  overview            -- when to use Java vs MATLAB vs MPh Python
  java_api            -- com.comsol.model class hierarchy + idioms
  matlab_livelink     -- mli usage from a MATLAB session
  mph_python          -- MPh Python (third-option, lab default for sweeps)
  lab_examples        -- model1.java / model2.java / ForJavaScript.java
  compile_workflow    -- COMSOL_compile.ps1 recipe + comsolbatch
  when_to_use_which   -- decision table + perf notes
  pitfalls            -- license, encoding, escaping, version drift
  servers             -- comsolserver / comsolmphserver / batch ports
  references          -- doc paths + lab files

Use:
    from radia_mcp.interop.comsol_livelink_knowledge import get_comsol_livelink
    print(get_comsol_livelink())            # overview
    print(get_comsol_livelink("java_api"))
    print(get_comsol_livelink("all"))
"""

from __future__ import annotations


# -----------------------------------------------------------------
# Section content
# -----------------------------------------------------------------

OVERVIEW = """\
## COMSOL LiveLink -- overview (Sugahara Lab stance)

COMSOL Multiphysics is a commercial FEA package whose .mph project
file is, at its core, a serialized Java object graph. COMSOL exposes
that graph through three programmable front-ends -- Java, MATLAB
(mli), and Python (MPh) -- which all manipulate the same Model
object. Choose by ergonomics, not by capability:

| Front-end | Speed   | Ergonomics | Lab use case                          |
|-----------|---------|------------|----------------------------------------|
| Java      | Fastest | Verbose    | Reproducible exports from the GUI      |
| MATLAB    | Fast    | Familiar   | Lab members already in MATLAB workflow |
| MPh (Py)  | Fast    | Pythonic   | Parameter sweeps, CI, radia-mcp glue   |

**All three** call the same `com.comsol.model.Model` API under the
hood -- the Java SDK is the source of truth, and the lab's existing
.java files (model1.java, model2.java, ForJavaScript.java) can be
re-run from any of them.

### Lab pipeline (current)

```
COMSOL GUI                              MATLAB / Python / Java
  |                                       |
  v                                       v
  build model in GUI -- File > Save As ---+
  > model.java                            |
  |                                       |
  | compile (comsolcompile model.java)    |
  v                                       |
  model.class  <----------------- run via java -cp comsol.jar model
  |
  +-- load into GUI again (File > Open > .class) for visual check
  +-- or run headless via comsolbatch -class model
  +-- or load over comsolmphserver from MATLAB / MPh
```

The lab has standardized the **GUI -> Save As Java -> edit -> recompile
-> run** loop. Parameter sweeps and batch jobs use MPh Python on top of
comsolmphserver. MATLAB LiveLink is available but used less in 2026
because MPh covers the same ground with a smaller dependency surface
(no MATLAB license needed on the worker machine).

### Where this knowledge sits in radia-mcp

- `radia_mcp.interop` (this module) -- documentation only; does NOT
  invoke COMSOL. COMSOL is commercial and not auto-detectable across
  user installs the way OpenSCAD/FreeCAD are.
- For lab-internal solver pipelines (Radia + NGSolve + Cubit), there
  is no COMSOL dependency. COMSOL is treated as a **reference
  comparison code** (cross-check answers, replicate published
  benchmarks), not as a primary solver.

### Why document this at all

Lab members regularly need to:
1. Port a COMSOL .mph model to NGSolve/Radia (read the .java, mimic
   the geometry / BC / physics, run side-by-side comparison).
2. Drive COMSOL headlessly for parameter sweeps (sweep gap_length,
   compute Bz, regress against NGSolve).
3. Compile + open a GUI-edited model on a server with no display.

The three knowledge sections (`java_api`, `matlab_livelink`,
`mph_python`) plus the working `lab_examples` cover these.
"""


JAVA_API = """\
## com.comsol.model -- Java API class hierarchy

The COMSOL Java SDK exposes ~5000 classes; in practice >95 % of
scripting needs the dozen-class subgraph below. Every `.java` file
exported from the GUI is a `public class XXX extends nothing,
implements nothing` with two static methods (`run()` and `main`).

### Standard preamble

```java
import com.comsol.model.*;
import com.comsol.model.util.*;

/** Model exported on <DATE> by COMSOL 6.4.0.293. */
public class model1 {

  public static Model run() {
    Model model = ModelUtil.create("Model");
    model.modelPath("O:\\\\");
    model.label("model1.mph");
    // ... geometry / physics / mesh / study / sol / result ...
    return model;
  }

  public static void main(String[] args) {
    run();
  }
}
```

### Core classes (the dozen you actually use)

| Class / package                 | Purpose                              |
|---------------------------------|--------------------------------------|
| `ModelUtil`                     | Static factory: create / load / save |
| `Model`                         | Top-level model handle               |
| `ModelParam` (`model.param()`)  | Global parameters with units         |
| `ModelComponent` (`model.component(tag)`) | Component (geom + mesh + physics) |
| `GeomSequence` (`model.geom(tag)`)        | Geometry sequence (`comp.geom("geom1")`) |
| `GeomFeature` (`geom.feature(tag)`)       | One geometry primitive / op |
| `MeshSequence` (`model.mesh(tag)`)        | Mesh sequence (FreeQuad, FreeTet, Swept) |
| `PhysicsFeature`                | Physics interface (e.g. CoefficientFormPDE) |
| `StudyFeature` (`model.study(tag)`)       | Study (Stationary, Eigenfrequency, ...) |
| `SolverSequence` (`model.sol(tag)`)       | Solver dynamics (FullyCoupled, Direct, Iterative) |
| `ResultFeature` (`model.result(tag)`)     | Plot groups, evaluations |
| `NumericalFeature` (`result().numerical()`) | Numerical post-processing |

### The fluent navigation pattern

Every accessor returns the next-level object, so the API reads as a
chain:

```java
model.component("comp1")
     .geom("geom1")
     .create("blk1", "Block");
model.component("comp1")
     .geom("geom1")
     .feature("blk1")
     .set("size", new int[]{1, 1, 2});
model.component("comp1")
     .geom("geom1")
     .run("blk1");                       // run up to blk1
```

This is the form the GUI emits. Long chains can be cached as locals:

```java
GeomFeature blk = model.component("comp1").geom("geom1")
                       .create("blk1", "Block");
blk.set("size", new int[]{1, 1, 2});
```

### Setting values -- type discipline

| What            | How                                          |
|-----------------|----------------------------------------------|
| Scalar int      | `feat.set("h", 2)`                           |
| Scalar double   | `feat.set("r", 0.25)`                        |
| Vector int      | `feat.set("size", new int[]{1, 1, 2})`       |
| Vector double   | `feat.set("pos", new double[]{0.5, 0, 0})`   |
| String          | `feat.set("base", "center")`                 |
| String list     | `feat.set("input", new String[]{"blk1", "cyl1"})` |
| String matrix   | `feat.set("table", new String[][]{{"-pi", "r"}, ...})` |
| With unit       | `model.param().set("g", "0.3[mm]", "Gap")`  |

The matrix form `new String[][]{...}` is by far the most common
gotcha -- the GUI uses it for Polygon tables, coefficient PDE
material tensors (`c = [[1/y, 0, 0, y]]`), and parametric ranges.

### Physics / study / sol / result triplet

A normal solve has 4 objects in lockstep:

```java
// physics
model.component("comp1").physics().create("c", "CoefficientFormPDE", "geom1");

// study
model.study().create("std1");
model.study("std1").create("stat", "Stationary");

// solver sequence (auto-derived from the study)
model.sol().create("sol1");
model.sol("sol1").attach("std1");
model.sol("sol1").create("st1", "StudyStep");
model.sol("sol1").create("v1", "Variables");
model.sol("sol1").create("s1", "Stationary");
model.sol("sol1").feature("s1").create("fc1", "FullyCoupled");

// run
model.study("std1").runNoGen();    // build mesh + solve, skip GUI plots
// or
model.sol("sol1").runAll();        // run only the solver chain

// post-processing
model.result().create("pg1", "PlotGroup2D");
model.result("pg1").create("surf1", "Surface");
```

### Solver picks (lab examples)

Direct PARDISO (model1.java):
```java
model.sol("sol1").feature("s1").create("d1", "Direct");
model.sol("sol1").feature("s1").feature("d1").set("linsolver", "pardiso");
model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "d1");
```

Iterative GMRES + AMG (model2.java):
```java
model.sol("sol1").feature("s1").create("i1", "Iterative");
model.sol("sol1").feature("s1").feature("i1").set("linsolver", "gmres");
model.sol("sol1").feature("s1").feature("i1").set("prefuntype", "left");
model.sol("sol1").feature("s1").feature("i1").create("mg1", "Multigrid");
model.sol("sol1").feature("s1").feature("i1").feature("mg1").set("prefun", "saamg");
model.sol("sol1").feature("s1").feature("fc1").set("linsolver", "i1");
```

The `prefun = "saamg"` is COMSOL's smoothed-aggregation AMG. For
HCurl problems the Radia stack prefers Compact HX (see CLAUDE.md);
COMSOL's choice is internal and not exposed beyond the
`prefun` string.

### Numerical post-processing

```java
double[][] xy = model.result().numerical()
                     .create("pev1", "EvalPoint")
                     .setIndex("expr", "Az", 0)
                     .getReal();
```

In practice the lab evaluates fields via `numerical()` and writes
the result to a CSV for downstream comparison with Radia/NGSolve.

### Parameters with the lab's `setX` getter pattern

The lab's `model1.java` adds Java-level mutators that let an
outer Java/Python driver tweak parameters without re-editing the
`.mph` file:

```java
private static double gap_length = 0.3;  // default 0.3 mm

public static void setGapLength(double length) { gap_length = length; }
public static double getGapLength()           { return gap_length; }

public static Model run() {
    Model model = ModelUtil.create("Model");
    model.param().set("gap_length", gap_length + "[mm]", "Air gap length");
    // ...
}
```

This is the lab's idiomatic way to make a recompiled-once .class file
serve a parameter sweep -- the outer driver calls
`model1.setGapLength(0.4); model1.run();` per case.
"""


MATLAB_LIVELINK = """\
## LiveLink for MATLAB (mli)

LiveLink for MATLAB ("mli") is COMSOL's add-on that puts the same
`com.comsol.model.Model` object into a MATLAB workspace as a Java
handle. The license is sold separately from base COMSOL.

### Starting an mli session

Two routes:

A. From the COMSOL launcher: pick "COMSOL Multiphysics with MATLAB"
   -- launches MATLAB and pre-connects to a hidden comsolserver.
   The MATLAB command window already has `mphload`, `mphsave`,
   `mphmodel`, etc. on the path.

B. From a plain MATLAB: start `comsolserver` (or `comsolmphserver`)
   in a separate terminal, then in MATLAB:

```matlab
addpath('C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\mli');
mphstart(2036);   % connect to server on TCP port 2036
```

### Core mli commands

| Command                              | Java equivalent           |
|--------------------------------------|---------------------------|
| `model = mphload('model.mph')`       | `ModelUtil.load(...)`     |
| `model = mphload('Model.class')`     | `<Cls>.run()` (load .class) |
| `mphsave(model, 'out.mph')`          | `model.save("out.mph")`   |
| `mphmodel`                           | list loaded models        |
| `mphtags(model, 'geom')`             | list geom tags            |
| `mphmesh(model)`                     | visualize the mesh        |
| `mphgeom(model)`                     | visualize the geometry    |
| `mphplot(model, 'pg1')`              | plot a result group       |
| `mphint2(model, 'Az', 'surface', 2)` | integrate Az on domain 2  |
| `mpheval(model, 'Az', 'edim', 'point', 'selection', [0 0])` | evaluate at point |

### Mutating the model object

The MATLAB handle IS a Java object -- the API is identical to the
Java SDK, just with MATLAB syntax for method calls:

```matlab
% set a parameter
model.param.set('gap_length', '0.3[mm]', 'Air gap length');

% create a geometry feature
model.geom('geom1').create('blk1', 'Block');
model.geom('geom1').feature('blk1').set('size', [1 1 2]);
model.geom('geom1').run('blk1');

% run the study
model.study('std1').runNoGen();

% read a result
data = mphint2(model, 'Az', 'surface', 2);
```

Note `'std1'` (single quotes -- MATLAB strings) instead of `"std1"`
(Java strings); aside from that the chain looks the same.

### Why the lab uses MATLAB less than Java + MPh

1. **License surcharge**: LiveLink for MATLAB is a paid add-on on
   top of the COMSOL base license.
2. **Worker-machine licensing**: every machine running an mli
   parameter sweep needs both a COMSOL token and a MATLAB license.
3. **MPh Python covers it**: 95 % of mli features are available via
   MPh + Python with no extra MATLAB cost.

MATLAB LiveLink is still useful when:
- The user is already in a MATLAB workflow (signal processing,
  controls).
- The result post-processing leans on toolboxes (Statistics,
  Optimization).
- A figure must drop straight into a MATLAB live script.
"""


MPH_PYTHON = """\
## MPh (Python) -- third option, lab default for sweeps

`MPh` (https://mph.readthedocs.io/, https://github.com/MPh-py/MPh)
is a third-party MIT-licensed Python package that wraps the COMSOL
Java API by speaking the comsolmphserver protocol from JPype. It
is **not** the same as the official "LiveLink for Python" (LL4Py,
COMSOL 6.3+); MPh predates LL4Py, has wider version support
(5.4+), and is the lab's default.

### Install

```bash
pip install mph
# requires JDK + a working COMSOL install on the same machine
```

### Start the server, connect from Python

```powershell
# in a separate PowerShell window:
& 'C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolmphserver.exe' -port 2036
```

```python
import mph

client = mph.Client(port=2036)              # connect to running server
model = client.load("O:/COMSOL/model1.mph") # OR
model = client.load("O:/COMSOL/model1.class")
```

The lab's `COMSOL_compile.ps1` automates this connection -- pass
`-ServerPort 2036` and the script will load the freshly compiled
`.class` into the live server (see `compile_workflow`).

### Mutate, solve, post

```python
# parameters
model.parameter("gap_length", "0.4[mm]")    # mph helpfully exposes scalars

# study
model.solve("std1")                          # equivalent to study.runNoGen()

# evaluate
B_avg = model.evaluate("intop1(B_z)")
print(f"Average Bz = {B_avg} T")

# save
model.save("model_g0.4.mph")
client.remove(model)                         # release the model handle
```

### Parameter sweep recipe (lab pattern)

```python
import mph, csv

client = mph.Client(port=2036)
model = client.load("O:/COMSOL/model1.mph")

results = []
for g_mm in [0.2, 0.3, 0.4, 0.5, 0.6]:
    model.parameter("gap_length", f"{g_mm}[mm]")
    model.solve("std1")
    Bz = model.evaluate("intop1(B_z)")
    results.append((g_mm, Bz))

with open("sweep.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["gap_length_mm", "Bz_avg_T"])
    w.writerows(results)

client.remove(model)
client.disconnect()
```

### Why MPh is the lab default

| Feature                          | MPh | mli  | Pure Java |
|----------------------------------|-----|------|-----------|
| Free / OSS                       | yes | no   | yes       |
| Python ecosystem (numpy, scipy)  | yes | no   | partial   |
| Same machine as COMSOL only      | no (any TCP host) | yes | yes |
| Compile cost                     | 0   | 0    | per-file  |
| Worker license per sweep         | 1xCOMSOL | 1xCOMSOL+1xMATLAB | 1xCOMSOL |
| Direct Java jar access (escape)  | yes (via mph.session.java) | no | n/a |

### Versions

- MPh 1.2+ -- COMSOL 5.4 through 6.4 verified.
- COMSOL 6.3 ships official "LiveLink for Python" (LL4Py). MPh has
  not yet been retired; the lab continues with MPh because the
  install matrix is simpler and the script compatibility goes back
  to COMSOL 5.4 archived models.
"""


LAB_EXAMPLES = """\
## Lab examples on S:\\COMSOL\\

The Sugahara Lab maintains a small set of canonical .java files
exported from the COMSOL GUI. These are the "Rosetta Stone" --
read them whenever you need to remember the Java idiom for a
specific feature.

### S:/COMSOL/2026_01_14_Java経由の実行練習/

| File                           | What it shows                          |
|--------------------------------|----------------------------------------|
| `model1.java`  (14.3 KB)       | 2D rotor-slot CoefficientFormPDE + Direct (PARDISO) |
| `model1.class` (14.3 KB)       | compiled bytecode (output of comsolcompile) |
| `model1.mph`   (1.6 MB)        | saved .mph after one run               |
| `model2.java`  (13.9 KB)       | same geometry, Iterative (GMRES + SAAMG) |
| `model2.class` (13.8 KB)       | compiled bytecode                      |
| `COMSOL_compile.ps1` (2 KB)    | the lab's compile + load-via-MPh recipe |

Highlights from model1.java:
- Parameter-driven geometry (`r_rotor`, `r_stator`, `gap_length`).
- 11 Polygon features (`pol1` ... `pol11`) -- one per rotor slot.
- CoefficientFormPDE with `c = [[1/y, 0, 0, y]]` (axisymmetric
  curl-curl Laplacian for A_phi in (r,z); the `1/y` is the same
  Henrotte-style 1/r weighting documented in CLAUDE.md's
  "Axisymmetric FE" policy -- COMSOL handles it via plain integration
  with no basis-substitution because its Gauss order is high enough
  to absorb the singularity for typical mesh densities).
- Per-slot source term `f = +y` / `f = -y` (positive/negative coil).
- Periodic conditions `pc1/pc2/pc3` linking pairs of boundary edges.
- Mesh: `FreeQuad` (free quadrilateral), `hauto = 3` (predefined coarse).
- Java `setGapLength()` / `getGapLength()` mutators for outer driver.

Highlights from model2.java:
- Identical geometry, demonstrates the **solver-only swap**
  -- replace `Direct + pardiso` with `Iterative + gmres + Multigrid +
  saamg` while leaving the rest of the chain untouched.
- This is the lab's idiomatic A/B benchmark: build geometry once in
  the GUI, save-as-Java twice with different solver picks, time
  both via `comsolbatch`.

### S:/COMSOL/03_マニュアル/2022_12_13-14_スタートアップセミナ/

| File                                       | What it shows                |
|--------------------------------------------|------------------------------|
| `data-startup1/ForJavaScript.java` (1.7 KB) | Smallest "minimal Block + Cylinder + Mirror + Array + Difference" example -- the COMSOL training course's "hello world" for Java scripting |
| `data-startup1/ForJavaScript.mph`           | matching .mph                |
| `data-startup1/geometry_programming.mph`    | parameterized geometry demo  |
| `data-startup2/ForJavaScript_saved_Model.mph` | model saved after running the .class |

ForJavaScript.java is the cleanest place to start reading: 60 lines
that build a hex-array-of-cylinders subtracted from a slab, demonstrating
`Block`, `Cylinder`, `Copy`, `Mirror`, `Union`, `Array`, `Difference`
-- the whole 3D constructive-solid-geometry vocabulary in one file.

### S:/COMSOL/99_Tips/Comsol_Tips.docx

Symmetry-plane cheat sheet (bilingual EN+JP):
- Magnetic Insulation: applied where current flows **normal** to
  the plane. When using voltage excitation, divide V by 2 per
  symmetry plane. When using current excitation, leave I, but
  post-process voltage *= 2 per plane.
- Perfect Magnetic Conductor: applied where current flows
  **tangentially** to the plane. When using current excitation,
  divide I by 2 per plane.

Apply the same algebra when porting to a Radia/NGSolve symmetry
quarter or eighth model.

### Lab knowledge propagation rule

The "Promotion-After-Verify" policy (see CLAUDE.md) applies here:
when a COMSOL-vs-Radia comparison settles on a number, write the
verified result into `examples/<topic>/README.md` AND into the
relevant `radia_mcp.<topic>` knowledge file. Do NOT leave the
verification only in a .mph file -- the next contributor has no
COMSOL license.
"""


COMPILE_WORKFLOW = """\
## Compile + run workflow

### The full pipeline

```
.java   --[comsolcompile]-->  .class   --[load]-->  Model object
                                          |
                                          +-> GUI (File > Open .class)
                                          +-> comsolbatch -class <Cls>
                                          +-> comsolmphserver + MPh client
```

### Lab's COMSOL_compile.ps1 (verbatim)

```powershell
# COMSOL_compile.ps1  --  compile a .java + open it
# Usage: .\\COMSOL_compile.ps1 <ModelName> [-ServerPort <port>]
# Example: .\\COMSOL_compile.ps1 model1
# Example: .\\COMSOL_compile.ps1 model1 -ServerPort 2036

param(
    [Parameter(Mandatory=$true, Position=0)]
    [string]$ModelName,

    [Parameter(Mandatory=$false)]
    [int]$ServerPort = 0
)

$ComsolCompile = "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolcompile.exe"
$WorkDir = "O:\\COMSOL"
$JavaFile  = Join-Path $WorkDir "$ModelName.java"
$ClassFile = Join-Path $WorkDir "$ModelName.class"

Set-Location $WorkDir
if (-not (Test-Path $JavaFile)) {
    Write-Error "Error: $JavaFile not found"
    exit 1
}

Write-Host "Compiling $JavaFile ..." -ForegroundColor Cyan
& $ComsolCompile $JavaFile

if ($LASTEXITCODE -eq 0 -and (Test-Path $ClassFile)) {
    Write-Host "Compilation successful: $ClassFile" -ForegroundColor Green

    if ($ServerPort -gt 0) {
        # Load via comsolmphserver + MPh
        $pythonScript = @"
import mph
client = mph.Client(port=$ServerPort)
print('Connected to COMSOL server')
model = client.load('$ClassFile'.replace('\\\\', '/'))
print(f'Model loaded: {model}')
"@
        python -c $pythonScript
    } else {
        Write-Host "=== Next Steps ===" -ForegroundColor Yellow
        Write-Host "Option 1: Open in running COMSOL GUI manually"
        Write-Host "  File > Open > $ClassFile"
        Write-Host ""
        Write-Host "Option 2: Start comsolmphserver and use -ServerPort"
        Write-Host "  1. Start: comsolmphserver -port 2036"
        Write-Host "  2. Run:   .\\COMSOL_compile.ps1 $ModelName -ServerPort 2036"
    }
} else {
    Write-Error "Compilation failed"
    exit 1
}
```

Source: `S:/COMSOL/2026_01_14_Java経由の実行練習/COMSOL_compile.ps1`

### Three ways to actually RUN a compiled .class

#### A. Open in the GUI (interactive)

```
File > Open > model1.class
```

The GUI calls `model1.run()` internally, then deserializes the
returned `Model` and lets you inspect / edit it. This is the
standard developer loop.

#### B. comsolbatch (headless, single run)

```powershell
& "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolbatch.exe" `
    -class model1 -classpath O:\\COMSOL `
    -outputfile O:\\COMSOL\\model1_out.mph
```

Flags:
- `-class <Cls>` : the FQCN of the class to run (must extend the
                   pattern `public static Model run()`).
- `-classpath`   : directory containing the .class
- `-inputfile`   : (alt) load a .mph instead of a .class
- `-outputfile`  : write the post-solve .mph here
- `-batchlog`    : redirect log to file
- `-np <N>`      : number of compute threads

#### C. comsolmphserver + Python/MATLAB client (long-lived)

```powershell
# In window A: start the server
& "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolmphserver.exe" `
    -port 2036
```

```python
# In window B (Python): connect and drive
import mph
client = mph.Client(port=2036)
model = client.load("O:/COMSOL/model1.class")
model.solve("std1")
client.disconnect()
```

This is the model the lab uses for parameter sweeps -- one server,
many `client.load()` calls with different parameters, no per-case
COMSOL startup cost.

### Directly running the .class with plain java (works but rarely used)

```powershell
$cp = "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\plugins\\com.comsol.model.jar;"
$cp += "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\plugins\\com.comsol.model.util.jar;"
$cp += "O:\\COMSOL"
& java -cp $cp model1
```

This works because `model1` has a `main(String[] args)` that calls
`run()`. However: the actual jar list COMSOL needs is ~30 entries;
in practice `comsolbatch` is easier than building the classpath
manually.

### Sweep loop with the lab's mutator pattern

```powershell
# Compile once
.\\COMSOL_compile.ps1 model1

# Sweep via comsolbatch's -pname / -plist
$gaps = 0.2, 0.3, 0.4, 0.5, 0.6
foreach ($g in $gaps) {
    & "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolbatch.exe" `
        -class model1 -classpath O:\\COMSOL `
        -pname gap_length -plist "$g[mm]" `
        -outputfile "O:\\COMSOL\\model1_g$g.mph" `
        -batchlog "O:\\COMSOL\\model1_g$g.log"
}
```

The `-pname` / `-plist` flags inject parameter overrides without
needing the Java-level `setGapLength()` mutator -- the lab uses
both depending on whether the sweep needs to touch parameters that
exist as COMSOL `model.param()` entries (use `-pname`/`-plist`) or
needs to vary geometry that is NOT a COMSOL parameter (use the Java
`setX()` mutator).
"""


WHEN_TO_USE_WHICH = """\
## When to use which front-end

### Decision table

| Situation                                     | Recommendation               |
|-----------------------------------------------|------------------------------|
| GUI just exported a .java, want to verify it  | Open .class in GUI (option A) |
| One-shot batch on cluster                     | `comsolbatch -class` (option B) |
| Parameter sweep, 10-10000 cases               | MPh Python + comsolmphserver  |
| Same machine has MATLAB, want figures fast    | LiveLink for MATLAB (mli)    |
| Need to wrap into radia-mcp tool              | MPh Python                   |
| Live coding in Jupyter                        | MPh Python                   |
| Read-only inspection of an .mph from CI       | MPh Python                   |
| Need to mutate geometry that is NOT a COMSOL parameter | Java-level setter + recompile |
| Cross-machine driver (LAB drives 100号機 COMSOL) | MPh over TCP                 |

### Performance notes

- **Compile cost**: `comsolcompile` is ~5-15 s for a medium .java.
  Recompile only when the .java changes; otherwise the .class is
  cheap to load.
- **Server startup**: `comsolmphserver -port 2036` takes ~10-30 s
  the first time (Java VM warm-up + license check). Subsequent
  client connections are sub-second.
- **License**: every concurrent `comsolbatch` and every active
  `comsolmphserver` consumes one COMSOL Multiphysics token. The
  GUI consumes one too. The lab has a finite token pool -- the
  `client.disconnect()` discipline matters for sweeps.
- **MATLAB**: mli ADDITIONALLY needs a MATLAB Computational Toolbox
  token; budget accordingly.

### Re-running the same model across versions

A .java exported by COMSOL 6.4 will usually compile and run on 6.3
and 6.2 with no edits -- the API is largely additive. A .mph saved
on 6.4 will silently upgrade when opened on 6.4+; older COMSOL
versions refuse to load it.

For lab archival: keep both the .java (text, version-control-able)
AND the .mph (binary, faster to open in GUI). The .java is the
canonical artifact; the .mph is the convenience cache.

### MPh vs LL4Py (COMSOL 6.3+)

LL4Py is COMSOL's official Python LiveLink, shipped from 6.3 onward.
The lab continues to use MPh because:
- MPh works on 5.4-6.4 (LL4Py is 6.3+ only).
- MPh's API is dictionary-of-everything (read `model.parameters` as
  a dict); LL4Py mirrors the Java fluent style (more verbose).
- Switching costs are real -- existing lab scripts use MPh idioms.

When LL4Py matures and the install base passes 6.3, the lab will
re-evaluate. For now: **MPh is the recommended Python entry point**.
"""


PITFALLS = """\
## Pitfalls (license, encoding, escaping, version drift)

### 1. License token exhaustion

Every of these consumes ONE Multiphysics token:
- An open GUI window
- A running `comsolbatch` process
- A running `comsolmphserver`
- A running `comsol server -multi` worker

The lab pool is finite. In MPh sweeps, ALWAYS:
```python
client.remove(model)
client.disconnect()
```
in a `finally:` block. Use `try / finally` around every `client.load()`.

### 2. Windows path escaping in .java

A .java exported by the GUI on Windows contains:
```java
model.modelPath("O:\\\\");
```
The double-backslash is correct Java string escaping for a single
`\\`. The COMSOL_compile.ps1 helper passes single backslashes via
PowerShell quoting, but if you HAND-EDIT a .java in Python:
```python
java_src = java_src.replace("O:\\", "O:\\\\\\\\")   # 1 backslash -> 2 in Java
```

### 3. cp932 console encoding (Sugahara Lab specific)

The lab's Windows machines default to cp932 (Shift-JIS). COMSOL log
output through `comsolbatch -batchlog` is fine, but if you tee
through PowerShell on the console you may see mojibake for any
Japanese parameter names. Either:
- Set `[Console]::OutputEncoding = [System.Text.Encoding]::UTF8` in the
  PowerShell session before launching, OR
- Use English-only parameter names in COMSOL (preferred).

### 4. Apostrophes in unit strings

```java
model.param().set("g", "0.3[mm]", "Gap");        // OK
model.param().set("g", "0.3 [mm]", "Gap");       // OK
model.param().set("g", "0.3 mm", "Gap");         // WRONG -- not a unit
```

Square brackets are required.

### 5. Tag uniqueness within a sequence

Tags (`"blk1"`, `"std1"`, `"pg1"`) must be unique within their
parent sequence. If you `create("blk1", ...)` twice, COMSOL throws
`com.comsol.util.exceptions.FlException: Tag is already used`. The
GUI auto-renames; programmatic `create()` does not.

### 6. selection().set(...) is 1-indexed

```java
model.component("comp1").physics("c").feature("cfeq2").selection().set(2);
```
The `2` is a 1-indexed domain number (or boundary / edge / vertex
number for lower-dimensional selections). The numbering is assigned
by the geometry kernel after `geom.run()` and CAN CHANGE if you
re-order features. The lab's recommended discipline (matches
`CLAUDE.md`'s Journal File Portability Policy for Cubit):
- **Do not** hardcode domain numbers in long-lived .java; instead
  use Selection nodes (`model.selection().create(...)`) with
  geometric predicates (`selectBox`, `selectAdjacent`, ...).

In practice, GUI exports DO hardcode domain numbers -- the lab
accepts this for one-shot benchmarks but rewrites with Selection
nodes for any production sweep.

### 7. Version drift (.mph forward-only)

A .mph saved on COMSOL 6.4 cannot be opened on 6.3 or earlier.
This is an asymmetric drift -- newer COMSOL opens older .mph
automatically, but never the reverse. The .java text export is
the portable artifact.

### 8. modelPath sets the WORKING directory, not a file

```java
model.modelPath("O:\\\\");   // not "O:/COMSOL/model1.mph"
```
`modelPath` is the directory where COMSOL resolves relative paths
(e.g. external CSV imports). For the .mph file itself, use
`model.save("model1.mph")`.

### 9. study.runNoGen() vs sol.runAll()

- `study("std1").runNoGen()` -- runs the study including
  auto-generating any missing solver / mesh sequences.
- `sol("sol1").runAll()` -- runs ONLY the named solver sequence;
  fails silently if the mesh or variables are not already built.

For a fresh load, prefer `runNoGen()`. For a recompute after
mutating only a parameter, `sol("sol1").runAll()` is faster.

### 10. No transactions

There is no `try { ... } commit` / `rollback` for the Model object.
A failed `create()` may leave the model in a half-built state.
The lab convention: serialize the model BEFORE mutating, and
re-load on failure.

```python
import mph
client = mph.Client(port=2036)
backup = client.load("model.mph")
try:
    backup.parameter("g", "0.4[mm]")
    backup.solve("std1")
    backup.save("model_committed.mph")
except Exception:
    backup = client.load("model.mph")   # reload the on-disk version
```
"""


SERVERS = """\
## comsolserver / comsolmphserver / comsolbatch

Three server-mode binaries ship with COMSOL Multiphysics. They are
NOT interchangeable; pick by use case.

### comsolserver (license server)

Manages license tokens. Started automatically by the install on
license-server machines. Lab users normally never touch this.

### comsolmphserver (mph protocol -- MPh / mli)

Long-lived headless COMSOL kernel that accepts connections from
MPh Python clients and MATLAB mli clients.

```powershell
& "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolmphserver.exe" `
    -port 2036 -np 8 -multi on
```

Flags:
- `-port`    : TCP port (default 2036)
- `-np`      : threads per solve
- `-multi`   : `on` lets multiple clients share one server (each
               gets its own Model handle, but they share the JVM
               and the COMSOL license token)
- `-silent`  : suppress banner

The lab convention: one comsolmphserver per worker machine, all
sweeps share it.

### comsolbatch (one-shot headless run)

Compiles + runs a single .java / .class / .mph and exits.

```powershell
& "C:\\Program Files\\COMSOL\\COMSOL64\\Multiphysics\\bin\\win64\\comsolbatch.exe" `
    -class model1 -classpath O:\\COMSOL `
    -outputfile O:\\COMSOL\\model1_out.mph `
    -batchlog  O:\\COMSOL\\model1.log
```

Useful flags (full set: `comsolbatch -help`):
- `-inputfile`  : path to .mph (alternative to `-class`)
- `-class`      : FQCN of a compiled Java class
- `-classpath`  : directory containing the .class
- `-job <tag>`  : the StudyJob tag to run (default: run all studies)
- `-pname`/`-plist` : parameter override (single name + comma list)
- `-paramfile`  : CSV of parameter sweep cases
- `-np`         : compute threads
- `-batchlog`   : log file
- `-outputfile` : save the post-solve .mph here

### When to use which

| Use case                              | Binary             |
|---------------------------------------|--------------------|
| One off solve from CI                 | comsolbatch        |
| Interactive sweep from notebook       | comsolmphserver + MPh client |
| Multi-tenant lab sweep machine        | comsolmphserver -multi on    |
| License server (admin)                | comsolserver       |
"""


REFERENCES = """\
## References

### Lab files (Sugahara Lab, S: drive)

| Path | What |
|------|------|
| `S:/COMSOL/2026_01_14_Java経由の実行練習/model1.java` | Rotor-slot 2D, Direct PARDISO solver |
| `S:/COMSOL/2026_01_14_Java経由の実行練習/model2.java` | Same geometry, Iterative GMRES+SAAMG |
| `S:/COMSOL/2026_01_14_Java経由の実行練習/COMSOL_compile.ps1` | Compile + MPh-load PowerShell helper |
| `S:/COMSOL/2026_01_14_Java経由の実行練習/Batch/` | (reserved for batch sweep scripts; currently empty) |
| `S:/COMSOL/COMSOL_Launcher.ps1` | Remote-PC COMSOL GUI launcher (WMI + Task Scheduler) |
| `S:/COMSOL/03_マニュアル/COMSOL_Multiphysics.pdf` | Main COMSOL reference (16.8 MB) |
| `S:/COMSOL/03_マニュアル/2022_12_13-14_スタートアップセミナ/.../ForJavaScript.java` | Minimal 3D CSG hello world |
| `S:/COMSOL/99_Tips/Comsol_Tips.docx` | Symmetry-plane current/voltage scaling cheat sheet |

### COMSOL official docs (paths on a default Windows install)

| Path                                                         | What                 |
|--------------------------------------------------------------|----------------------|
| `C:/Program Files/COMSOL/COMSOL64/Multiphysics/doc/pdf/...`  | All product manuals  |
| `COMSOL_ProgrammingReferenceManual.pdf`                      | Java API reference (full class tree) |
| `COMSOL_LiveLinkForMATLABUsersGuide.pdf`                     | mli reference         |
| `COMSOL_ApplicationProgrammingGuide.pdf`                     | Building apps         |
| `bin/win64/comsolcompile.exe`                                | Java compiler         |
| `bin/win64/comsolbatch.exe`                                  | Headless runner       |
| `bin/win64/comsolmphserver.exe`                              | MPh/mli server        |
| `mli/`                                                       | MATLAB add-path target |

### Third-party

| Item                                  | URL                                   |
|---------------------------------------|----------------------------------------|
| MPh -- COMSOL Python wrapper          | https://mph.readthedocs.io/            |
| MPh source                            | https://github.com/MPh-py/MPh          |

### Related radia-mcp knowledge

- `radia_mcp.cubit` -- the lab's hex-mesh authoring (the Radia side
  of any COMSOL/Radia cross-check).
- `radia_mcp.radia_ngsolve` -- the FEM solver the lab uses for
  most NGSolve-side validations.
- `radia_mcp.motor.onelab_knowledge` -- ONELAB/GetDP, the
  open-source motor FEA reference; complements COMSOL for
  reproducibility-mandated comparisons.
- `radia_mcp.interop` (this subpackage) -- routes other CAD MCPs
  (build123d, OpenSCAD, FreeCAD) through Cubit. COMSOL is NOT
  routed through this dispatcher; COMSOL stays at the documentation
  layer because it is closed-source and cannot be auto-detected
  reliably.
"""


# -----------------------------------------------------------------
# Dispatch
# -----------------------------------------------------------------

TOPICS = {
    "overview":            OVERVIEW,
    "java_api":            JAVA_API,
    "matlab_livelink":     MATLAB_LIVELINK,
    "mph_python":          MPH_PYTHON,
    "lab_examples":        LAB_EXAMPLES,
    "compile_workflow":    COMPILE_WORKFLOW,
    "when_to_use_which":   WHEN_TO_USE_WHICH,
    "pitfalls":            PITFALLS,
    "servers":             SERVERS,
    "references":          REFERENCES,
}


def get_comsol_livelink(topic: str = "overview") -> str:
    """Return COMSOL LiveLink (Java + MATLAB + MPh Python) knowledge.

    Args:
        topic: section name (see TOPICS keys) or "all" for everything.
    """
    t = topic.strip().lower()
    if t == "all":
        return "\n\n---\n\n".join(TOPICS[k] for k in TOPICS)
    if t not in TOPICS:
        valid = ", ".join(sorted(TOPICS))
        return (
            f"Unknown topic {t!r}. Valid topics: {valid}, or 'all'\n"
        )
    return TOPICS[t]
