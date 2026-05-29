"""Per-DOF Z_s + |H_t| map (3-panel hotspot visualisation).

Renders the per-DOF complex surface impedance Z_s and the driver
|H_t| field on the workpiece surface, unrolled to (theta, z) at
the cylindrical side wall.  Three panels:

    (a) Re Z_s [mOhm]        loss component (drives P_wp_local)
    (b) Im Z_s [mOhm]        reactive (tracks mu(|H|) more strongly)
    (c) |H_t| [A/m]          driver field (66x max/min contrast)

For the IH headline benchmark (steel cylinder, 50 kHz,
I_port = 100 A) this figure visually explains the per-element
ESIM contribution: |H_t| spans 16-1084 A/m (66x max/min, 2.86x
max/mean), so the surface is strongly non-uniform; where |H_t| is
largest the BH curve drives the steel deeper into saturation,
lowering the local Re Z_s (13.4-19.1 mOhm, min/mean 0.75x).  The
local loss density ~ Re(Z_s) * |H_t|^2 is dominated by the |H_t|^2
spread.  Scalar mesh-RMS Karl collapses this map to ONE Z_s; per-
element Karl preserves it, giving P_wp = 18.75 W vs scalar 30.51 W
(-38.5%).  (The earlier "94x / 2.1x" figures were a Galerkin-
localisation artifact, corrected to the triangle-gradient |H_t|.)

Inputs:
    JSON with esim_per_panel_Z_s_real/imag and esim_per_panel_H_t
    arrays (produced by calc_inductance.py --esim-per-panel; needs
    radia >= 4.67 for the corrected triangle-gradient |H_t|).
    Default path points at the committed
    sweep_data/I100_f50k_per_panel.json (the 100 A, 50 kHz headline
    case) -- regenerable from the repo, no C:/temp dependency.

    Workpiece .vol mesh used to recover the per-DOF (x, y, z)
    coordinates.  Default points at samples/ih_bem_sample_p1.vol.

Outputs:
    Zs_per_dof_map.png       (3-panel figure)

Usage:
    python examples/ih_esim_benchmark/plot_zs_per_dof_map.py \\
        [JSON_PATH] [VOL_PATH]

The default paths are configured for the IGTE 2026 paper
reproducer; pass alternative paths to visualise a different run.
"""
import sys
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from mcp_server_document.graph.tools import apply_lab_style, lab_savefig

# Data Persistence Policy: default to the committed sweep_data case,
# not a transient C:/temp snapshot.
DEFAULT_JSON = Path(__file__).parent / "sweep_data" / "I100_f50k_per_panel.json"
DEFAULT_VOL = Path(__file__).resolve().parent.parent.parent / \
              "src" / "radia" / "panels" / "samples" / "ih_bem_sample_p1.vol"
OUT_PNG = Path(__file__).parent / "Zs_per_dof_map.png"


def main():
    json_path = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT_JSON
    vol_path = Path(sys.argv[2]) if len(sys.argv) > 2 else DEFAULT_VOL

    if not json_path.exists():
        print(f"ERROR: {json_path} not found.")
        print("Generate it by running:")
        print("  python src/radia/panels/calc_inductance.py \\")
        print("    --coil-step samples/ih_fem_kelvin_demo_coil.step \\")
        print("    --coil-solver peec --vol samples/ih_bem_sample_p1.vol \\")
        print("    --wp-label sibc --sigma 2e6 --mu-r 100 --half-thickness 0.005 \\")
        print("    --frequency 50000 --current 100.0 --coil-sigma 5.8e7 \\")
        print("    --impedance-model esim --bh-file samples/em_sample_bh.txt \\")
        print("    --esim-per-panel --h1-order 1 --wp-bem-backend intree-dense \\")
        print(f"    --output {json_path}")
        sys.exit(1)

    from ngsolve import Mesh, BND
    from ngsolve import TaskManager

    with open(json_path) as f:
        d = json.load(f)
    Zs_re = np.asarray(d["esim_per_panel_Z_s_real"])
    Zs_im = np.asarray(d["esim_per_panel_Z_s_imag"])
    if "esim_per_panel_H_t" not in d:
        raise RuntimeError(
            "JSON missing esim_per_panel_H_t array.  Re-run "
            "calc_inductance.py with radia >= 4.55.x.")
    H_t = np.asarray(d["esim_per_panel_H_t"])
    print(f"JSON: {json_path.name}, ndof = {len(Zs_re)}")
    print(f"Re(Z_s): {Zs_re.min()*1e3:.2f} -- {Zs_re.max()*1e3:.2f} mOhm "
          f"(ratio {Zs_re.max()/Zs_re.min():.2f}x)")
    print(f"Im(Z_s): {Zs_im.min()*1e3:.2f} -- {Zs_im.max()*1e3:.2f} mOhm "
          f"(ratio {Zs_im.max()/Zs_im.min():.2f}x)")
    print(f"|H_t|:   {H_t.min():.1f} -- {H_t.max():.1f} A/m "
          f"(ratio {H_t.max()/H_t.min():.2f}x)")

    with TaskManager():
        mesh = Mesh(str(vol_path))
        bnd_labels = list(mesh.GetBoundaries())
        sibc_idx = {i for i, n in enumerate(bnd_labels) if n.lower() == "sibc"}
        bnd_verts = set()
        for el in mesh.Elements(BND):
            if int(el.index) in sibc_idx:
                for v in el.vertices:
                    bnd_verts.add(int(v.nr))
        bnd_verts = sorted(bnd_verts)
        pts = np.array([list(mesh.ngmesh.Points()[i+1]) for i in bnd_verts])
        x, y, z = pts.T

        # Identify side wall (cylindrical, r ~ R_max) and unroll to (theta, z)
        r_cyl = np.sqrt(x**2 + y**2)
        R_max = float(r_cyl.max())
        side = r_cyl > 0.95 * R_max
        theta = np.arctan2(y[side], x[side])
        z_side = z[side] * 1e3
        theta_deg = np.degrees(theta)

        figsize = apply_lab_style(target="paper_double_column", aspect=0.30)
        fig, axes = plt.subplots(1, 3, figsize=figsize,
                                  constrained_layout=True)
        panels = [
            (r"$\mathrm{Re}\,Z_s$ (m$\Omega$)", Zs_re[side]*1e3, "inferno"),
            (r"$\mathrm{Im}\,Z_s$ (m$\Omega$)", Zs_im[side]*1e3, "viridis"),
            (r"$|H_t|$ (A/m)",                  H_t[side],       "magma"),
        ]
        panel_tags = ["(a)", "(b)", "(c)"]
        for ax, tag, (clabel, data, cmap) in zip(axes, panel_tags, panels):
            sc = ax.scatter(theta_deg, z_side, c=data, cmap=cmap,
                            s=6, edgecolors="none")
            cb = plt.colorbar(sc, ax=ax, pad=0.02, fraction=0.05)
            cb.set_label(clabel)
            ax.set_xlabel(r"$\theta$ (deg)")
            ax.set_xticks([-180, -90, 0, 90, 180])
            # Subfigure label below each panel (IEEE convention).
            ax.text(0.5, -0.32, tag, transform=ax.transAxes,
                    ha="center", va="top")
        axes[0].set_ylabel(r"$z$ (mm)")
        lab_savefig(fig, str(OUT_PNG.with_suffix("")))
        print(f"Saved {OUT_PNG}")


if __name__ == "__main__":
    main()
