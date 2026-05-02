"""Export a combined .msh v4.1 showing BOTH the BEM-A surface mesh
AND the PEEC filaments on the same coil, so they can be inspected
side-by-side in GMSH (open with ``gmsh combined.msh``).

Physical groups in the output file:
  body       (2D tris) -- BEM-A lateral surface mesh
  source     (2D tris) -- BEM-A source cap
  sink       (2D tris) -- BEM-A sink cap
  peec_fil   (1D lines) -- PEEC perimeter filament polylines (n_peri x N_seg)

Coil: ``tests/panels/golden/rect_torus_lofted_united.step`` (build123d
355deg multi-loft, 8 x 6 mm rect cross-section, R = 50 mm).
"""
from __future__ import annotations

import os
import sys

import gmsh
import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))

# Make peec_bundle and friends importable.
SRC_RADIA = os.path.join(ROOT, "src", "radia")
if SRC_RADIA not in sys.path:
    sys.path.insert(0, SRC_RADIA)


def main():
    step = os.path.join(ROOT, "tests", "panels", "golden",
                         "rect_torus_lofted_united.step")
    out_msh = os.path.join(HERE, "rect_torus_loft_combined.msh")

    # ------------------------------------------------------------
    # Get the BEM-A surface mesh (verts, tris, BND labels)
    # ------------------------------------------------------------
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
    # Mesh density tuned for VISUALISATION: caps and lateral panels
    # both subdivide into ~10 tris each so the rectangular cross-section
    # geometry is visually clear in GMSH.  (BEM-A solve at this density
    # would be slow in pure Python; for L extraction we re-run on a
    # coarser mesh in the test.)
    ngmesh = OCCGeometry(Glue(list(shape.faces))).GenerateMesh(
        mp=MeshingParameters(maxh=0.003, curvaturesafety=1.0,
                              segmentsperedge=1))
    mesh = Mesh(ngmesh)
    bnd_names = list(mesh.GetBoundaries())

    verts = np.array(
        [[mesh.vertices[i].point[k] for k in range(3)]
         for i in range(mesh.nv)], dtype=np.float64)
    tris_body = []
    tris_source = []
    tris_sink = []
    for el in mesh.Elements(BND):
        v3 = [int(v.nr) for v in el.vertices]
        lbl = bnd_names[el.index]
        if lbl == "source":
            tris_source.append(v3)
        elif lbl == "sink":
            tris_sink.append(v3)
        else:
            tris_body.append(v3)
    tris_body = np.array(tris_body, dtype=np.int64)
    tris_source = np.array(tris_source, dtype=np.int64)
    tris_sink = np.array(tris_sink, dtype=np.int64)
    print(f"BEM-A mesh: nv={mesh.nv}, "
          f"body={len(tris_body)}, source={len(tris_source)}, "
          f"sink={len(tris_sink)}")

    # ------------------------------------------------------------
    # Get PEEC filaments on the same STEP
    # ------------------------------------------------------------
    from radia.coil_from_cad import filaments_from_step
    topo = filaments_from_step(step, sigma=5.8e7, n_peri=16,
                                cad_units_per_meter=1.0)
    fps = topo["filament_paths"]
    n_fil = len(fps)
    n_seg_per = np.asarray(fps[0]).shape[0]
    print(f"PEEC filaments: n_peri={n_fil}, n_seg/fil={n_seg_per}, "
          f"total segments={n_fil * n_seg_per}")

    # ------------------------------------------------------------
    # Build combined .msh via gmsh API
    # ------------------------------------------------------------
    gmsh.initialize()
    gmsh.model.add("rect_torus_combined")

    # Create one entity per dim+role so we can group nicely.
    # Surfaces: 3 entities (body=1, source=2, sink=3) all "discrete".
    # Curves: 1 entity (peec_fil=10).
    SURF_BODY = gmsh.model.addDiscreteEntity(2, 1)
    SURF_SRC = gmsh.model.addDiscreteEntity(2, 2)
    SURF_SNK = gmsh.model.addDiscreteEntity(2, 3)
    CURVE_FIL = gmsh.model.addDiscreteEntity(1, 10)

    # Physical groups.
    gmsh.model.addPhysicalGroup(2, [SURF_BODY], tag=1, name="body")
    gmsh.model.addPhysicalGroup(2, [SURF_SRC], tag=2, name="source")
    gmsh.model.addPhysicalGroup(2, [SURF_SNK], tag=3, name="sink")
    gmsh.model.addPhysicalGroup(1, [CURVE_FIL], tag=10, name="peec_fil")

    # Nodes for the BEM-A mesh: 1-indexed.  We use the same global node
    # IDs across all surface entities since the verts are shared.
    bem_node_tags = list(range(1, mesh.nv + 1))
    bem_node_xyz = verts.flatten().tolist()
    # Add nodes attached to one of the surfaces (any will do; we just
    # want them in the global pool).  GMSH .msh v4 wants nodes per
    # entity, so we attach them all to BODY.
    gmsh.model.mesh.addNodes(2, SURF_BODY, bem_node_tags, bem_node_xyz)

    # Surface elements.  GMSH element type 2 = Tri3.
    def tri_elements(tris, start_tag):
        # gmsh.model.mesh.addElementsByType wants flat connectivity
        # array of node tags (1-indexed).
        n = len(tris)
        if n == 0:
            return start_tag
        elem_tags = list(range(start_tag, start_tag + n))
        node_tags = (tris + 1).flatten().tolist()
        return elem_tags, node_tags

    elem_id = 1
    if len(tris_body) > 0:
        et, nt = tri_elements(tris_body, elem_id)
        gmsh.model.mesh.addElementsByType(SURF_BODY, 2, et, nt)
        elem_id += len(et)
    if len(tris_source) > 0:
        et, nt = tri_elements(tris_source, elem_id)
        gmsh.model.mesh.addElementsByType(SURF_SRC, 2, et, nt)
        elem_id += len(et)
    if len(tris_sink) > 0:
        et, nt = tri_elements(tris_sink, elem_id)
        gmsh.model.mesh.addElementsByType(SURF_SNK, 2, et, nt)
        elem_id += len(et)

    # Filament line segments.  GMSH element type 1 = Line2 (2-node line).
    # Each segment has 2 endpoints; multiple segments share endpoints
    # within one filament but we keep them flat for simplicity.
    fil_node_tags = []
    fil_node_xyz = []
    fil_elem_tags = []
    fil_node_pairs = []  # for the 2-node-line connectivity
    next_node = mesh.nv + 1   # continue numbering after BEM nodes
    for f_idx, fp in enumerate(fps):
        segs = np.asarray(fp)   # (N_seg, 2, 3)
        for s_idx, (p0, p1) in enumerate(segs):
            n0_tag = next_node
            n1_tag = next_node + 1
            next_node += 2
            fil_node_tags.extend([n0_tag, n1_tag])
            fil_node_xyz.extend([float(p0[0]), float(p0[1]),
                                  float(p0[2]),
                                  float(p1[0]), float(p1[1]),
                                  float(p1[2])])
            fil_elem_tags.append(elem_id)
            elem_id += 1
            fil_node_pairs.extend([n0_tag, n1_tag])
    # Attach filament nodes to the curve entity.
    gmsh.model.mesh.addNodes(1, CURVE_FIL, fil_node_tags, fil_node_xyz)
    gmsh.model.mesh.addElementsByType(CURVE_FIL, 1, fil_elem_tags,
                                       fil_node_pairs)

    gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
    gmsh.write(out_msh)
    gmsh.finalize()
    print(f"\nsaved: {out_msh}  ({os.path.getsize(out_msh)} bytes)")
    print(f"\nOpen in GMSH:")
    print(f"  gmsh {out_msh}")
    print()
    print("In GMSH: Tools -> Visibility shows 4 physical groups:")
    print("  body, source, sink (surfaces), peec_fil (lines)")
    print("Toggle each group to compare BEM-A surface mesh vs PEEC filaments.")


if __name__ == "__main__":
    main()
