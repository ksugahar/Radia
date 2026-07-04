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

Legacy/helper registration as Windows file handler:
    python tools/vol_sol_viewer.py --register

Restore default Netgen association:
    python tools/vol_sol_viewer.py --unregister

Policy note:
    A raw Windows association to `netgen.exe "%1"` can open a blank GUI
    because it may not call Netgen's native mesh-loader path.  For robust
    double-click use either this helper (`radia-vol-viewer --register`) or
    an equivalent Netgen startup hook that runs `Ng_LoadMesh`, selects
    mesh/solution visual mode, redraws, and reads status.  `.sol` always
    needs the companion `.vol` because NGSolve GridFunction `.sol` files do
    not contain a mesh.
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

    # Each tuple: (label, builder, is_vector).  is_vector tells the viewer
    # whether to use scalfunction="<name>:0" + vecfunction="<name>" (vector
    # field) or scalfunction="<name>" alone (scalar field).
    candidates = [
        ("HDiv order=1",        lambda: HDiv (mesh, order=1),         True),
        ("HDiv order=2",        lambda: HDiv (mesh, order=2),         True),
        ("HCurl order=1",       lambda: HCurl(mesh, order=1),         True),
        ("HCurl order=2",       lambda: HCurl(mesh, order=2),         True),
        ("H1 order=1 scalar",   lambda: H1   (mesh, order=1),         False),
        ("H1 order=2 scalar",   lambda: H1   (mesh, order=2),         False),
        ("H1 order=1 dim=3",    lambda: H1   (mesh, order=1, dim=3),  True),
        ("H1 order=2 dim=3",    lambda: H1   (mesh, order=2, dim=3),  True),
    ]
    for label, mk, is_vector in candidates:
        try:
            fes = mk()
        except Exception as e:
            print(f"  {label}: construction failed ({type(e).__name__})", flush=True)
            continue
        if fes.ndof == n_real:
            print(f"  {label}: ndof={fes.ndof} MATCHES real (sol_size/8)",
                  flush=True)
            return fes, label, is_vector
        if fes.ndof == n_complex:
            print(f"  {label}: ndof={fes.ndof} MATCHES complex (sol_size/16) — trying complex FES",
                  flush=True)
            try:
                if "HCurl" in label:
                    fes_c = HCurl(mesh, order=fes.order, complex=True)
                elif "HDiv" in label:
                    fes_c = HDiv(mesh, order=fes.order, complex=True)
                elif "H1" in label and "dim" in label:
                    fes_c = H1(mesh, order=fes.order, dim=3, complex=True)
                else:
                    fes_c = H1(mesh, order=fes.order, complex=True)
                return fes_c, label + " complex", is_vector
            except Exception:
                return fes, label, is_vector

    print(f"  WARNING: no FES candidate matched (sol={sol_size} B = "
          f"{n_real} doubles or {n_complex} complex)", flush=True)
    return None, None, False


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

    fes, label, is_vector = _infer_fes(mesh, sol_abs)
    field_name = Path(sol_abs).stem.rsplit("_", 1)[-1]   # _B -> "B", _Jsurf -> "Jsurf"

    if fes is None:
        print("Could not determine FES — falling back to mesh-only display",
              flush=True)
        gfu = None
    else:
        gfu = GridFunction(fes)
        gfu.Load(sol_abs)
        print(f"  loaded {label}: ndof={fes.ndof} {'(vector)' if is_vector else '(scalar)'}",
              flush=True)
        Draw(gfu, mesh, field_name)
        print(f"  Draw(gfu, mesh, '{field_name}') registered", flush=True)

    def _show():
        if gfu is None:
            # Mesh-only fallback
            g.win.tk.eval("set selectvisual mesh")
            g.win.tk.eval("Ng_SetVisParameters")
            g.win.tk.eval("redraw")
            g.win.tk.eval("Ng_ReadStatus")
            return

        # Switch the viewer to "Solution" mode AND bind which functions
        # to use for the colour-mapped scalar and the vector overlay.
        # Without scalfunction the togl viewer renders uniform-coloured
        # surfaces (the 2026-05-24 "色がついていない" symptom).
        g.win.tk.eval("set selectvisual solution")
        if is_vector:
            # Vector field (HDiv / HCurl / H1 dim=3): colour by magnitude
            # (NGSolve registers scalar component slots :0 :1 :2 + a
            # synthetic "abs" component reached via visoptions.evaluate).
            g.win.tk.eval(f'set visoptions.scalfunction "{field_name}:0"')
            g.win.tk.eval(f'set visoptions.vecfunction "{field_name}"')
            g.win.tk.eval('set visoptions.evaluate "abs"')
        else:
            # Scalar field (H1 order=p): colour by value directly.
            g.win.tk.eval(f'set visoptions.scalfunction "{field_name}"')
            g.win.tk.eval('set visoptions.vecfunction "none"')
        # Surface projection is required on volume meshes -- otherwise
        # the togl viewer has no surface to paint the field on (it does
        # not slice into the volume by default).
        g.win.tk.eval("set visoptions.showsurfacesolution 1")
        g.win.tk.eval("set visoptions.subdivision 1")
        g.win.tk.eval("Ng_Vis_Set parameters")
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
    r"""Register .vol / .sol -> radia-vol-viewer-gui.exe (preferred) or
    pythonw + this_script.py (fallback for standalone tools/ checkout).

    Writes the association entries directly via the `winreg` stdlib
    module instead of cmd.exe `assoc` + `ftype`.  The latter mangles
    embedded quotes in the argv -> command-line round-trip: argv
    `['cmd','/c','ftype', 'Radia.VolViewer="C:\\path\\app.exe" "%1"']`
    reaches cmd.exe as `ftype Radia.VolViewer=\C:\path\app.exe\ \%1\`,
    which silently fails to register the ftype handler (assoc is set
    but ftype is empty).  Symptom seen 2026-05-24 deploying to 100号機
    + mdx: `radia-vol-viewer --register` prints "Done" but
    `ftype Radia.VolViewer` then returns "file type not found".

    Writes to HKLM\SOFTWARE\Classes (machine-wide), which is what the
    original `assoc`/`ftype` already targeted.  Requires admin.
    """
    import winreg

    gui_exe = _find_gui_script()
    if gui_exe:
        command_str = f'"{gui_exe}" "%1"'
        print("Registering file associations (console-less gui-script):")
        print(f"  Launcher: {gui_exe}")
    else:
        py_exe = sys.executable
        pyw = py_exe.replace("python.exe", "pythonw.exe")
        if not os.path.exists(pyw):
            pyw = py_exe
        this_script = os.path.abspath(__file__)
        command_str = f'"{pyw}" "{this_script}" "%1"'
        print("Registering file associations (pythonw fallback -- "
              "radia-vol-viewer-gui.exe not found):")
        print(f"  Python (windowless): {pyw}")
        print(f"  Script:              {this_script}")

    try:
        # .vol -> Radia.VolViewer
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\.vol",
            winreg.REG_SZ, "Radia.VolViewer",
        )
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\Radia.VolViewer\shell\open\command",
            winreg.REG_SZ, command_str,
        )
        # .sol -> Radia.SolViewer
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\.sol",
            winreg.REG_SZ, "Radia.SolViewer",
        )
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\Radia.SolViewer\shell\open\command",
            winreg.REG_SZ, command_str,
        )
    except PermissionError as e:
        print(f"ERROR: cannot write to HKLM\\SOFTWARE\\Classes ({e})",
              flush=True)
        print("  Re-run from an elevated (admin) shell.", flush=True)
        sys.exit(1)

    print()
    print("Done.  Double-click any .vol or .sol to view.")
    print()
    print("Verification:")
    print('  cmd /c "assoc .vol"            -> .vol=Radia.VolViewer')
    print('  cmd /c "ftype Radia.VolViewer" -> "<launcher>" "%1"')


def unregister_associations():
    """Restore default Netgen.VolFile association and drop .sol mapping.

    Uses winreg for the same quote-mangling reason as register_associations.
    """
    import winreg

    print("Restoring default Netgen.VolFile association...")
    py_exe = sys.executable
    netgen_command = f'"{py_exe}" -m netgen "%1"'
    try:
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\.vol",
            winreg.REG_SZ, "Netgen.VolFile",
        )
        winreg.SetValue(
            winreg.HKEY_LOCAL_MACHINE,
            r"SOFTWARE\Classes\Netgen.VolFile\shell\open\command",
            winreg.REG_SZ, netgen_command,
        )
        try:
            with winreg.OpenKey(
                winreg.HKEY_LOCAL_MACHINE,
                r"SOFTWARE\Classes",
                0, winreg.KEY_WRITE,
            ) as classes:
                winreg.DeleteKey(classes, ".sol")
        except FileNotFoundError:
            pass
    except PermissionError as e:
        print(f"ERROR: cannot write to HKLM\\SOFTWARE\\Classes ({e})",
              flush=True)
        print("  Re-run from an elevated (admin) shell.", flush=True)
        sys.exit(1)
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
