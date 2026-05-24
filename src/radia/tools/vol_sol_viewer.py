"""
vol_sol_viewer.py -- double-click viewer for .vol and .sol files.

User-side counterpart to upstream PR
https://github.com/NGSolve/netgen/pull/242 -- closed by author
(head repo deleted), NOT rejected by maintainers (zero review
comments).  This file replicates PR #242's `_vol_handler` verbatim
for .vol and adds an NGSolve-binary loader for .sol.

Architecture:

  .vol  -> Ng_LoadMesh "<path>"             (Tcl, native loader)
        -> set selectvisual mesh             (*** required ***)
        -> Ng_SetVisParameters / redraw / Ng_ReadStatus
        -> g.win.mainloop()

  .sol  -> [find companion .vol]             (Radia naming convention)
        -> ngsolve.Mesh(vol_path)            (so NGSolve sees the mesh)
        -> probe FES by ndof = sol_size/8    (HDiv > H1 scalar > H1 dim=3 > HCurl)
        -> GridFunction(fes); gfu.Load(sol)
        -> Draw(gfu, mesh, name)             (NGSolve registers viz)
        -> set selectvisual solution         (switch viewer to solution mode)
        -> redraw / Ng_ReadStatus
        -> g.win.mainloop()

Companion .vol discovery (.sol -> .vol):

  Radia panel output convention (calc_fem_kelvin / calc_peec etc.):
    <basename>_<solver>_fem.vol          <- the FEM mesh
    <basename>_<solver>_B.sol            <- B field      (HDiv o=1)
    <basename>_<solver>_Jsurf.sol        <- |J| at surf  (H1 scalar)
    <basename>_<solver>_qsurf.sol        <- charge       (H1 scalar)
    <basename>_<solver>_Jvec.sol         <- J vector     (HDivSurf or H1 dim=3)

  Heuristic: peel trailing _<word> suffixes from the stem, try
  `<stem>.vol` and `<stem>_fem.vol` at each level until a match is
  found.

Usage:
    python vol_sol_viewer.py path/to/file.vol
    python vol_sol_viewer.py path/to/field.sol

Register as Windows file handler:
    python tools/vol_sol_viewer.py --register

Restore default Netgen association:
    python tools/vol_sol_viewer.py --unregister
"""

import sys
import os
import argparse
from pathlib import Path


# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

def _tcl_path(p):
    """Tcl prefers forward slashes."""
    return str(p).replace("\\", "/")


def _find_companion_vol(sol_path):
    """Locate the .vol file that goes with sol_path.

    Tries in order:
      <stem>.vol
      <stem>_fem.vol                    (Radia panel-output convention)
      <stem-stripped>.vol               (peel trailing _word)
      <stem-stripped>_fem.vol
      ... repeat until no underscore remains

    Returns None if nothing found.
    """
    p = Path(sol_path)
    parent = p.parent
    stem = p.stem
    while True:
        for suffix in ("", "_fem"):
            cand = parent / f"{stem}{suffix}.vol"
            if cand.exists():
                return cand
        if "_" not in stem:
            return None
        stem = stem.rsplit("_", 1)[0]


def _infer_fes(mesh, sol_path):
    """Return (fes, label) by matching FES ndof to sol size.

    .sol files are raw float64 (or float64 pairs for complex) coefficient
    dumps from NGSolve's GridFunction.Save().  Match ndof to file size.
    """
    from ngsolve import H1, HCurl, HDiv

    sol_size = os.path.getsize(sol_path)
    n_real = sol_size // 8       # float64
    n_complex = sol_size // 16   # complex float64

    candidates = [
        ("HDiv order=1",        lambda: HDiv (mesh, order=1)),
        ("HDiv order=2",        lambda: HDiv (mesh, order=2)),
        ("HCurl order=1",       lambda: HCurl(mesh, order=1)),
        ("HCurl order=2",       lambda: HCurl(mesh, order=2)),
        ("H1 order=1 scalar",   lambda: H1   (mesh, order=1)),
        ("H1 order=2 scalar",   lambda: H1   (mesh, order=2)),
        ("H1 order=1 dim=3",    lambda: H1   (mesh, order=1, dim=3)),
        ("H1 order=2 dim=3",    lambda: H1   (mesh, order=2, dim=3)),
    ]
    for label, mk in candidates:
        try:
            fes = mk()
        except Exception as e:
            print(f"  {label}: construction failed ({type(e).__name__})", flush=True)
            continue
        if fes.ndof == n_real:
            print(f"  {label}: ndof={fes.ndof} MATCHES real (sol_size/8)",
                  flush=True)
            return fes, label
        if fes.ndof == n_complex:
            print(f"  {label}: ndof={fes.ndof} MATCHES complex (sol_size/16) — trying complex FES",
                  flush=True)
            # Re-build complex variant
            try:
                # Strip "(...)" if any; just try common complex flag
                if "HCurl" in label:
                    fes_c = HCurl(mesh, order=fes.order, complex=True)
                elif "HDiv" in label:
                    fes_c = HDiv(mesh, order=fes.order, complex=True)
                elif "H1" in label and "dim" in label:
                    fes_c = H1(mesh, order=fes.order, dim=3, complex=True)
                else:
                    fes_c = H1(mesh, order=fes.order, complex=True)
                return fes_c, label + " complex"
            except Exception:
                return fes, label

    # Nothing matched
    print(f"  WARNING: no FES candidate matched (sol={sol_size} B = "
          f"{n_real} doubles or {n_complex} complex)", flush=True)
    return None, None


# -------------------------------------------------------------------------
# View .vol
# -------------------------------------------------------------------------

def view_vol(vol_path):
    """Display .vol mesh -- PR #242 verbatim sequence."""
    import netgen.gui as g

    fp = _tcl_path(os.path.abspath(vol_path))
    print(f"load {fp}", flush=True)

    def _load():
        g.win.tk.eval(f'Ng_LoadMesh "{fp}"')
        g.win.tk.eval("set selectvisual mesh")
        g.win.tk.eval("Ng_SetVisParameters")
        g.win.tk.eval("redraw")
        g.win.tk.eval("Ng_ReadStatus")

    g.win.after(100, _load)
    g.win.mainloop()


# -------------------------------------------------------------------------
# View .sol
# -------------------------------------------------------------------------

def view_sol(sol_path, vol_override=None):
    """Display .vol + .sol via NGSolve gfu.Load() + Draw().

    .sol files are NGSolve GridFunction binary (gfu.Save() raw coefficient
    dumps), NOT Netgen's text "solution" format -- so Ng_ImportSolution
    cannot read them.  We load via NGSolve and let Draw() register the
    visualization, then switch the Tk viewer mode to "solution".
    """
    import netgen.gui as g
    from ngsolve import Mesh, GridFunction, Draw

    sol_abs = os.path.abspath(sol_path)
    if vol_override:
        vol = Path(os.path.abspath(vol_override))
    else:
        vol = _find_companion_vol(sol_abs)
    if vol is None or not vol.exists():
        print(f"ERROR: cannot find companion .vol for {sol_abs}", flush=True)
        print(f"  Pass --vol explicitly, or place .vol next to .sol", flush=True)
        sys.exit(1)

    print(f"mesh: {vol}", flush=True)
    print(f"sol:  {sol_abs}", flush=True)

    # NGSolve must own the mesh for the GridFunction it loads.  Loading via
    # ngsolve.Mesh() also populates Netgen's global mesh, so the Tk viewer
    # sees it too.
    mesh = Mesh(str(vol))
    print(f"  mesh ne={mesh.ne}, nv={mesh.nv}, "
          f"materials={mesh.GetMaterials()[:3]}{'...' if len(mesh.GetMaterials())>3 else ''}",
          flush=True)

    fes, label = _infer_fes(mesh, sol_abs)
    field_name = Path(sol_abs).stem.rsplit("_", 1)[-1]   # _B -> "B", _Jsurf -> "Jsurf"

    if fes is None:
        print("Could not determine FES — falling back to mesh-only display",
              flush=True)
        gfu = None
    else:
        gfu = GridFunction(fes)
        gfu.Load(sol_abs)
        print(f"  loaded {label}: ndof={fes.ndof}", flush=True)
        Draw(gfu, mesh, field_name)
        print(f"  Draw(gfu, mesh, '{field_name}') registered", flush=True)

    def _show():
        if gfu is not None:
            g.win.tk.eval("set selectvisual solution")
        else:
            g.win.tk.eval("set selectvisual mesh")
        g.win.tk.eval("Ng_SetVisParameters")
        g.win.tk.eval("redraw")
        g.win.tk.eval("Ng_ReadStatus")

    g.win.after(200, _show)
    print(f"Field '{field_name}' displayed. Close the Netgen window to exit.",
          flush=True)
    g.win.mainloop()


# -------------------------------------------------------------------------
# Windows file association
# -------------------------------------------------------------------------

def _find_gui_script():
    """Return path to radia-vol-viewer-gui.exe if shipped by this install.

    The gui-script .exe is built by setuptools from
    `[project.gui-scripts]` in pyproject.toml; it launches via
    pythonw.exe (no console window).  Falling back to
    `pythonw.exe + this_script.py` works too but is brittle (path
    contains spaces, two arguments) -- prefer the .exe when present.
    """
    import shutil
    exe = shutil.which("radia-vol-viewer-gui")
    if exe and os.path.exists(exe):
        return exe
    scripts_dir = os.path.join(os.path.dirname(sys.executable), "Scripts")
    candidate = os.path.join(scripts_dir, "radia-vol-viewer-gui.exe")
    if os.path.exists(candidate):
        return candidate
    return None


def register_associations():
    """Register .vol / .sol -> radia-vol-viewer-gui.exe (preferred) or
    pythonw + this_script.py (fallback for standalone tools/ checkout)."""
    import subprocess

    gui_exe = _find_gui_script()
    if gui_exe:
        cmd_vol = f'Radia.VolViewer="{gui_exe}" "%1"'
        cmd_sol = f'Radia.SolViewer="{gui_exe}" "%1"'
        print("Registering file associations (console-less gui-script):")
        print(f"  Launcher: {gui_exe}")
    else:
        py_exe = sys.executable
        pyw = py_exe.replace("python.exe", "pythonw.exe")
        if not os.path.exists(pyw):
            pyw = py_exe
        this_script = os.path.abspath(__file__)
        cmd_vol = f'Radia.VolViewer="{pyw}" "{this_script}" "%1"'
        cmd_sol = f'Radia.SolViewer="{pyw}" "{this_script}" "%1"'
        print("Registering file associations (pythonw fallback -- "
              "radia-vol-viewer-gui.exe not found):")
        print(f"  Python (windowless): {pyw}")
        print(f"  Script:              {this_script}")

    subprocess.run(["cmd", "/c", "assoc", ".vol=Radia.VolViewer"], check=False)
    subprocess.run(["cmd", "/c", "ftype", cmd_vol], check=False)
    subprocess.run(["cmd", "/c", "assoc", ".sol=Radia.SolViewer"], check=False)
    subprocess.run(["cmd", "/c", "ftype", cmd_sol], check=False)
    print()
    print("Done.  Double-click any .vol or .sol to view.")
    print()
    print("Verification:")
    print('  cmd /c "assoc .vol"            -> .vol=Radia.VolViewer')
    print('  cmd /c "ftype Radia.VolViewer" -> launcher path')


def unregister_associations():
    """Restore default Netgen.VolFile association."""
    import subprocess

    print("Restoring default Netgen.VolFile association...")
    py_exe = sys.executable
    subprocess.run(["cmd", "/c", "assoc", ".vol=Netgen.VolFile"], check=False)
    subprocess.run(
        ["cmd", "/c", "ftype",
         f'Netgen.VolFile="{py_exe}" -m netgen "%1"'],
        check=False,
    )
    subprocess.run(["cmd", "/c", "assoc", ".sol="], check=False)
    print("Done.")


# -------------------------------------------------------------------------
# CLI
# -------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("path", nargs="?",
                        help="Path to .vol or .sol file to view")
    parser.add_argument("--vol", default=None,
                        help="Override companion .vol path (only with .sol)")
    parser.add_argument("--register", action="store_true",
                        help="Register this script as the .vol/.sol Windows handler")
    parser.add_argument("--unregister", action="store_true",
                        help="Restore default Netgen .vol association")
    args = parser.parse_args(argv)

    if args.register:
        register_associations()
        return 0
    if args.unregister:
        unregister_associations()
        return 0

    if not args.path:
        parser.print_help()
        return 1

    path = args.path
    if path.lower().endswith(".vol"):
        view_vol(path)
    elif path.lower().endswith(".sol"):
        view_sol(path, vol_override=args.vol)
    else:
        print(f"Unknown file extension: {path}")
        print("Supported: .vol, .sol")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
