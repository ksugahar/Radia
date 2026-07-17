"""Regenerate the three .vol meshes the ESIM cross-formulation lane
consumes.  Requires Coreform Cubit (license) -- run once before
`run_cross_formulation.py`.

All three are regenerated from TRACKED generators (the .vol files
themselves are gitignored):

  1. samples/ih_fem_kelvin_demo.vol   (dia 20 x 20 small cylinder,
     + ih_fem_kelvin_demo_coil.step)  -- from samples/ih_fem_kelvin_demo.py
  2. samples/ih_bem_sample_p1.vol     (dia 50 x 25 BIE surface mesh)
     -- from samples/ih_bem_sample.jou
  3. <lane>/ih_fem_kelvin_50mm.vol    (dia 50 x 25 FEM-Kelvin mesh)
     -- from this lane's make_ih_fem_kelvin_50mm.py

Mesh (1) is shared by both the BIE and FEM small-cylinder runs (same
surface).  Meshes (2) and (3) are INDEPENDENT surface discretisations
of the same 50 mm cylinder -- comparing P_wp across them isolates the
outer-field formulation from the mesh.

NOTE (LAB / 100号機): `coreform_cubit -batch <script>` fails here with
an RLM license error (-102) in headless batch, but the in-process
Cubit Python API (`cubit.init([... , '-commandplugindir', plugins])`)
authenticates correctly.  This script therefore drives everything in a
single cubit.init() session rather than spawning batch subprocesses.

Usage:
    python make_meshes.py                # regenerate all that are missing
    python make_meshes.py --force        # regenerate all unconditionally
    python make_meshes.py --only fem50   # just the 50 mm FEM mesh
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
SAMPLES = os.path.join(REPO, "src", "radia", "panels", "samples")

DEMO_PY = os.path.join(SAMPLES, "ih_fem_kelvin_demo.py")
DEMO_VOL = os.path.join(SAMPLES, "ih_fem_kelvin_demo.vol")
BEM_JOU = os.path.join(SAMPLES, "ih_bem_sample.jou")
BEM_VOL = os.path.join(SAMPLES, "ih_bem_sample_p1.vol")
FEM50_PY = os.path.join(HERE, "make_ih_fem_kelvin_50mm.py")
FEM50_VOL = os.path.join(HERE, "ih_fem_kelvin_50mm.vol")

CUBIT_BIN = os.environ.get(
    "CUBIT_BIN", r"C:\Program Files\Coreform Cubit 2025.12\bin")
CUBIT_PLUGINS = os.path.join(CUBIT_BIN, "plugins")


def _init_cubit():
    if CUBIT_BIN not in sys.path:
        sys.path.insert(0, CUBIT_BIN)
    import cubit
    cubit.init(["cubit", "-nojournal", "-batch", "-nographics",
                "-commandplugindir", CUBIT_PLUGINS])
    return cubit


def _exec_generator(path, extra_globals=None):
    """exec() a Cubit-embedded generator script in a namespace with no
    __file__ (so its own NameError fallback picks up cwd / env override)."""
    g = {"__name__": "__cubit_gen__"}
    if extra_globals:
        g.update(extra_globals)
    with open(path, "r", encoding="utf-8") as f:
        code = f.read()
    exec(compile(code, path, "exec"), g)


def gen_demo(cubit):
    """samples/ih_fem_kelvin_demo.py writes its outputs next to getcwd()
    when __file__ is absent; chdir to samples so the .vol + coil .step
    land in the samples dir."""
    cubit.cmd("reset")
    prev = os.getcwd()
    os.chdir(SAMPLES)
    try:
        _exec_generator(DEMO_PY)
    finally:
        os.chdir(prev)
    _require(DEMO_VOL, "demo (dia 20 small cylinder)")


def gen_bem(cubit):
    cubit.cmd("reset")
    cubit.cmd('play "%s"' % BEM_JOU.replace("\\", "/"))
    cubit.cmd('export netgen "%s" order 1 overwrite'
              % BEM_VOL.replace("\\", "/"))
    _require(BEM_VOL, "BEM (dia 50 BIE surface mesh)")


def gen_fem50(cubit):
    cubit.cmd("reset")
    os.environ["RADIA_50MM_VOL_OUT"] = FEM50_VOL.replace("\\", "/")
    _exec_generator(FEM50_PY)
    _require(FEM50_VOL, "FEM-Kelvin (dia 50 volume mesh)")


def _require(path, label):
    if not os.path.isfile(path):
        raise RuntimeError("mesh generation FAILED for %s: %s not written"
                           % (label, path))
    print("  OK %-40s %s" % (label, path))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--force", action="store_true",
                    help="regenerate even if the .vol already exists")
    ap.add_argument("--only", choices=["demo", "bem", "fem50"],
                    help="regenerate only one mesh")
    args = ap.parse_args()

    targets = {
        "demo": (DEMO_VOL, gen_demo),
        "bem": (BEM_VOL, gen_bem),
        "fem50": (FEM50_VOL, gen_fem50),
    }
    if args.only:
        targets = {args.only: targets[args.only]}

    todo = [(k, fn) for k, (vol, fn) in targets.items()
            if args.force or not os.path.isfile(vol)]
    if not todo:
        print("All meshes present (use --force to regenerate).")
        for k, (vol, _) in targets.items():
            print("  have %-6s %s" % (k, vol))
        return 0

    cubit = _init_cubit()
    print("Regenerating:", ", ".join(k for k, _ in todo))
    for _, fn in todo:
        fn(cubit)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
