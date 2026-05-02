"""Export an interactive 3D HTML view of the combined BEM-A surface mesh
+ PEEC filaments via Plotly.  Open the resulting .html in any browser
to rotate / zoom / toggle layers without needing GMSH GUI.

Robust to environments where ``gmsh.fltk.run()`` doesn't show a window
(headless servers, SSH without X / RDP, pip-gmsh missing FLTK on
Windows).
"""
from __future__ import annotations

import os
import sys

import numpy as np
import plotly.graph_objects as go


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
SRC_RADIA = os.path.join(ROOT, "src", "radia")
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)


def main():
    step = os.path.join(ROOT, "tests", "panels", "golden",
                         "rect_torus_lofted_united.step")

    # --- BEM-A surface mesh ---
    from netgen.occ import OCCGeometry, Glue
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh, BND
    geo = OCCGeometry(step)
    shape = geo.shape
    fs = sorted(enumerate(shape.faces), key=lambda x: x[1].mass)
    ca, fa = fs[0]
    cb, fb = fs[1]
    if abs(fa.center[1]) < abs(fb.center[1]):
        fa.name, fb.name = "source", "sink"
    else:
        fa.name, fb.name = "sink", "source"
    for i, f in enumerate(shape.faces):
        if i not in (ca, cb):
            f.name = "body"
    ngmesh = OCCGeometry(Glue(list(shape.faces))).GenerateMesh(
        mp=MeshingParameters(maxh=0.012, curvaturesafety=1.0,
                              segmentsperedge=1))
    mesh = Mesh(ngmesh)
    bnd_names = list(mesh.GetBoundaries())

    verts = np.array(
        [[mesh.vertices[i].point[k] for k in range(3)]
         for i in range(mesh.nv)], dtype=np.float64)
    tris_body = []
    tris_src = []
    tris_snk = []
    for el in mesh.Elements(BND):
        v3 = [int(v.nr) for v in el.vertices]
        lbl = bnd_names[el.index]
        if lbl == "source":
            tris_src.append(v3)
        elif lbl == "sink":
            tris_snk.append(v3)
        else:
            tris_body.append(v3)
    tris_body = np.asarray(tris_body, dtype=np.int64)
    tris_src = np.asarray(tris_src, dtype=np.int64)
    tris_snk = np.asarray(tris_snk, dtype=np.int64)

    # --- PEEC filaments ---
    from radia.coil_from_cad import filaments_from_step
    topo = filaments_from_step(step, sigma=5.8e7, n_peri=16,
                                cad_units_per_meter=1.0)
    fps = topo["filament_paths"]   # list of (n_seg, 2, 3)

    # --- Plotly traces ---
    fig = go.Figure()

    def mesh3d(tris, name, color, opacity):
        if len(tris) == 0:
            return
        fig.add_trace(go.Mesh3d(
            x=verts[:, 0], y=verts[:, 1], z=verts[:, 2],
            i=tris[:, 0], j=tris[:, 1], k=tris[:, 2],
            color=color, opacity=opacity, name=name,
            flatshading=True, showscale=False,
            hoverinfo="name"))

    mesh3d(tris_body, "body (BEM-A 298 tris)", "lightblue", 0.55)
    mesh3d(tris_src, "source cap (BEM-A 2 tris)", "red", 0.95)
    mesh3d(tris_snk, "sink cap (BEM-A 2 tris)", "darkorange", 0.95)

    # PEEC filaments: each filament is a polyline; we can collapse the
    # (n_seg, 2, 3) into (n_seg+1, 3) by taking unique points along the
    # spine of one filament.
    def filament_polyline(fp_arr):
        # fp_arr: (n_seg, 2, 3); first point of each segment + last point
        pts = list(fp_arr[:, 0])
        pts.append(fp_arr[-1, 1])
        return np.asarray(pts)

    for k, fp in enumerate(fps):
        line = filament_polyline(np.asarray(fp))
        fig.add_trace(go.Scatter3d(
            x=line[:, 0], y=line[:, 1], z=line[:, 2],
            mode="lines",
            line=dict(color="black", width=4),
            name=("PEEC filaments (n_peri=16)" if k == 0 else None),
            legendgroup="peec_fil",
            showlegend=(k == 0),
            hoverinfo="name"))

    fig.update_layout(
        title=("rect_torus_lofted_united -- BEM-A surface mesh + "
               "PEEC perimeter filaments"),
        scene=dict(
            aspectmode="data",
            xaxis_title="x [m]", yaxis_title="y [m]", zaxis_title="z [m]"),
        legend=dict(itemclick="toggle", itemdoubleclick="toggleothers"))
    out_html = os.path.join(HERE, "rect_torus_loft_combined.html")
    fig.write_html(out_html, include_plotlyjs="cdn", full_html=True)
    print(f"saved: {out_html}  ({os.path.getsize(out_html)} bytes)")
    print()
    print("Open in any browser to interactively rotate / zoom / toggle layers.")
    print("Click legend entries to hide/show each physical group.")


if __name__ == "__main__":
    main()
