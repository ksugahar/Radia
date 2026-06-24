"""
Shared pipeline helper: build123d (OCC) -> Netgen (tet) -> Gmsh (post).

Lab policy (2026-04-19):
    - CAD:  build123d (Python / OCCT) — main
    - Mesh: Netgen (tet) via netgen.occ — tet-only main path
    - Post: gmsh (Python API) — visualization-only, no mesh generation
    - FreeCAD / ParaView / VTK: not used

The pipeline is designed to be called repeatedly from a sweep script so
the lab mcp-servers (build123d / gmsh) get exercised on many inputs.

Stages
------
1. build123d Part   -> .brep                 (export_brep)
2. .brep            -> netgen.occ.OCCGeometry -> generate tet mesh
3. netgen Mesh      -> .msh (Gmsh2 Format)    (ngsolve Mesh.Export)
4. gmsh             -> opens .msh, attaches a dummy scalar field
                       (proxy for post-processing from a solver),
                       writes the augmented .msh back out.

Each stage returns a small dict so the sweep script can aggregate
pass/fail + geometry/mesh stats per run.
"""

from __future__ import annotations

import os
import json
import time
import traceback
from pathlib import Path


def _stage_cad(part, out_dir: Path, label: str) -> dict:
    """Stage 1: export build123d Part to .brep, capture CAD stats."""
    from build123d import export_brep

    brep = out_dir / f"{label}.brep"
    export_brep(part, str(brep))

    edges = part.edges()
    faces = part.faces()
    bb = part.bounding_box()
    stats = {
        "stage": "cad",
        "ok": part.is_valid,
        "brep": str(brep),
        "volume": float(part.volume),
        "area": float(part.area),
        "faces": len(faces),
        "edges": len(edges),
        "min_edge": float(min(e.length for e in edges)),
        "bbox_size": [round(s, 6) for s in bb.size],
    }
    return stats


def _stage_mesh(brep_path: Path, out_dir: Path, label: str,
                maxh: float = 5.0) -> dict:
    """Stage 2+3: Netgen mesh the BREP, export as .msh (Gmsh2 Format)."""
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh as NGMesh
    from ngsolve import TaskManager

    t0 = time.perf_counter()
    geo = OCCGeometry(str(brep_path))
    with TaskManager():
        ng_mesh = geo.GenerateMesh(maxh=maxh)
        dt_gen = time.perf_counter() - t0

        msh_path = out_dir / f"{label}.msh"
        # Netgen Mesh exposes .Export(filename, format_name)
        # Format names include "Gmsh Format" (v4) and "Gmsh2 Format" (v2.2)
        ng_mesh.Export(str(msh_path), "Gmsh2 Format")

        # Wrap in NGSolve Mesh for quick element counts (cheaper than iterating).
        ngs = NGMesh(ng_mesh)
        stats = {
            "stage": "mesh",
            "ok": True,
            "msh": str(msh_path),
            "gen_seconds": round(dt_gen, 3),
            "maxh": maxh,
            "nv": ngs.nv,
            "ne": ngs.ne,
            "nface": ngs.nface,
            "nedge": ngs.nedge,
        }
        return stats


def _stage_post(msh_path: Path, out_dir: Path, label: str) -> dict:
    """Stage 4: gmsh reads the .msh, attaches a dummy scalar field,
    writes an augmented .msh for visualization/inspection.

    The dummy field is f(x,y,z) = x + 2y + 3z — a trivial stand-in for
    solver output (Radia B-field, NGSolve potential, etc.).  Real
    workflows replace this with actual data via gmsh.view.addModelData.
    """
    import gmsh

    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.open(str(msh_path))
        model_name = gmsh.model.getCurrent()

        # Collect nodes of the loaded mesh.
        node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
        # node_coords is a flat [x0,y0,z0, x1,y1,z1, ...] array.
        n_nodes = len(node_tags)
        xyz = [
            (node_coords[3 * i], node_coords[3 * i + 1], node_coords[3 * i + 2])
            for i in range(n_nodes)
        ]
        values = [[x + 2.0 * y + 3.0 * z] for (x, y, z) in xyz]

        view = gmsh.view.add(f"{label}_dummy_scalar")
        gmsh.view.addModelData(
            tag=view,
            step=0,
            modelName=model_name,
            dataType="NodeData",
            tags=list(node_tags),
            data=values,
            numComponents=1,
        )

        out_msh = out_dir / f"{label}_post.msh"
        gmsh.write(str(out_msh))

        # Sanity: probe the dummy field at the origin (well-defined for
        # any geometry centered near 0). Expect ~0 (probe_distance == 0
        # means the sample point is inside a mesh element).
        probe_val, probe_dist = gmsh.view.probe(view, 0.0, 0.0, 0.0)
        probe_val = list(probe_val) if len(probe_val) else None

        stats = {
            "stage": "post",
            "ok": True,
            "msh_post": str(out_msh),
            "view_tag": view,
            "n_nodes": n_nodes,
            "probe_at_origin": probe_val,
            "probe_distance": round(probe_dist, 6),
        }
    finally:
        gmsh.finalize()
    return stats


def run_pipeline(part, out_dir, label: str, maxh: float = 5.0) -> dict:
    """Run the full build123d -> Netgen -> Gmsh pipeline.

    Parameters
    ----------
    part : build123d.Part (or any build123d Shape with volume)
    out_dir : str or Path — directory to write .brep, .msh, *_post.msh
    label : str — filename stem (one label = one run)
    maxh : float — Netgen global mesh size

    Returns
    -------
    dict with keys: label, status ("ok"/"error"), stages (list of stage
    dicts), and (on error) error (traceback).
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    record = {"label": label, "status": "ok", "stages": []}
    try:
        s1 = _stage_cad(part, out_dir, label)
        record["stages"].append(s1)
        s2 = _stage_mesh(Path(s1["brep"]), out_dir, label, maxh=maxh)
        record["stages"].append(s2)
        s3 = _stage_post(Path(s2["msh"]), out_dir, label)
        record["stages"].append(s3)
    except Exception:
        record["status"] = "error"
        record["error"] = traceback.format_exc()
    return record


def save_record(record: dict, out_dir) -> Path:
    """Persist the run record as JSON next to the outputs."""
    out_dir = Path(out_dir)
    p = out_dir / f"{record['label']}.json"
    p.write_text(json.dumps(record, indent=2, default=str))
    return p


# ---------------------------------------------------------------------------
# Multi-region flow
# ---------------------------------------------------------------------------
#
# build123d doesn't emit material names into STEP/BREP that Netgen will pick
# up, so we:
#   1. Export each region as its own Part in a build123d Compound to STEP
#      (carries solid geometry + order, but region names are dropped).
#   2. Reload with netgen.occ, iterate `.solids` in order, call `.mat(name)`,
#      then `Glue(...)` them into a single multi-domain shape. This is the
#      key: Glue preserves material names at the mesh-generation level.
#   3. Export to Gmsh v4 ("Gmsh Format"). Per-domain tags are preserved but
#      names do not survive the exporter.
#   4. Open in gmsh and call `setPhysicalName(3, tag, name)` against the
#      sorted tag list to re-attach region names, then write back.
# Region order is the contract: regions[i] owns physical-group sorted-index i.


def _stage_cad_multi(regions, out_dir: Path, label: str) -> dict:
    """Stage 1 (multi): write one STEP carrying all regions as a Compound."""
    from build123d import Compound, export_step

    parts = [p for (p, _n) in regions]
    names = [n for (_p, n) in regions]

    # Compound preserves the ordering of children, which is our region-id contract.
    compound = Compound(children=parts)
    compound.label = label

    step = out_dir / f"{label}.step"
    export_step(compound, str(step))

    stats = {
        "stage": "cad",
        "ok": all(p.is_valid for p in parts),
        "step": str(step),
        "regions": [
            {"index": i, "name": names[i],
             "volume": float(parts[i].volume),
             "faces": len(parts[i].faces()),
             "is_valid": bool(parts[i].is_valid)}
            for i in range(len(parts))
        ],
    }
    return stats


def _stage_mesh_multi(step_path: Path, region_names, out_dir: Path,
                      label: str, maxh: float = 5.0) -> dict:
    """Stage 2+3 (multi): load STEP, Glue named solids, mesh, export Gmsh v4."""
    from netgen.occ import OCCGeometry, Glue
    from ngsolve import Mesh as NGMesh
    from ngsolve import TaskManager

    t0 = time.perf_counter()
    geo = OCCGeometry(str(step_path))
    sols = list(geo.shape.solids)
    if len(sols) != len(region_names):
        raise RuntimeError(
            f"STEP has {len(sols)} solids but {len(region_names)} region "
            f"names were given. Order is the contract — check the "
            f"build123d Compound children list."
        )
    for s, name in zip(sols, region_names):
        s.mat(name)
    geo2 = OCCGeometry(Glue(sols))
    ng_mesh = geo2.GenerateMesh(maxh=maxh)
    dt_gen = time.perf_counter() - t0

    msh_path = out_dir / f"{label}.msh"
    # Gmsh v4 format preserves per-domain physical tags (needed for
    # multi-region); v2.2 flattens them.
    ng_mesh.Export(str(msh_path), "Gmsh Format")

    ngs = NGMesh(ng_mesh)
    materials = list(ngs.GetMaterials())
    stats = {
        "stage": "mesh",
        "ok": True,
        "msh": str(msh_path),
        "gen_seconds": round(dt_gen, 3),
        "maxh": maxh,
        "materials": materials,
        "nv": ngs.nv,
        "ne": ngs.ne,
        "nface": ngs.nface,
        "nedge": ngs.nedge,
    }
    return stats


def _stage_post_multi(msh_path: Path, region_names, out_dir: Path,
                      label: str) -> dict:
    """Stage 4 (multi): open in gmsh, re-tag physical group names (by
    sorted-tag order matching `region_names`), attach a region-id scalar
    field so each region renders in a different color, write back."""
    import gmsh

    gmsh.initialize(["gmsh", "-nopopup"])
    try:
        gmsh.open(str(msh_path))
        model_name = gmsh.model.getCurrent()

        # Physical groups for volumes (dim=3). Sort by tag so the order
        # matches the build123d Compound children order (= region_names order).
        vol_groups = sorted(
            gmsh.model.getPhysicalGroups(dim=3),
            key=lambda pg: pg[1],
        )
        if len(vol_groups) != len(region_names):
            raise RuntimeError(
                f"gmsh reports {len(vol_groups)} volume groups but "
                f"{len(region_names)} region names were given."
            )
        for (dim, tag), name in zip(vol_groups, region_names):
            gmsh.model.setPhysicalName(dim, tag, name)

        # Build a per-element region-id field: each element gets the
        # 1-based index of its owning region. Renders as flat-shaded
        # regions in gmsh — handy for visual region QA.
        elem_tags = []
        elem_values = []
        for region_idx, (dim, tag) in enumerate(vol_groups, start=1):
            entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            for ent_tag in entities:
                e_types, e_tags, _ = gmsh.model.mesh.getElements(dim, ent_tag)
                for t_arr in e_tags:
                    for t in t_arr:
                        elem_tags.append(int(t))
                        elem_values.append([float(region_idx)])

        view = gmsh.view.add(f"{label}_region_id")
        gmsh.view.addModelData(
            tag=view,
            step=0,
            modelName=model_name,
            dataType="ElementData",
            tags=elem_tags,
            data=elem_values,
            numComponents=1,
        )

        out_msh = out_dir / f"{label}_post.msh"
        gmsh.write(str(out_msh))

        # Per-region node counts — useful for mesh-balance sanity checks.
        per_region = []
        for region_idx, (dim, tag) in enumerate(vol_groups, start=1):
            name = gmsh.model.getPhysicalName(dim, tag)
            ents = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
            n_elem = sum(
                sum(len(t) for t in gmsh.model.mesh.getElements(dim, e)[1])
                for e in ents
            )
            per_region.append({
                "index": region_idx,
                "physical_tag": tag,
                "name": name,
                "entities": list(map(int, ents)),
                "n_elem": n_elem,
            })

        stats = {
            "stage": "post",
            "ok": True,
            "msh_post": str(out_msh),
            "view_tag": view,
            "n_regions": len(vol_groups),
            "regions": per_region,
        }
    finally:
        gmsh.finalize()
    return stats


def run_pipeline_multi(regions, out_dir, label: str,
                       maxh: float = 5.0) -> dict:
    """Multi-region variant of run_pipeline.

    Parameters
    ----------
    regions : list of (build123d.Part, str)
        Ordered list of (part, region_name) tuples. Order is the contract:
        the i-th region corresponds to the i-th volume group in the output.
    out_dir : str or Path
    label : str
    maxh : float — Netgen global mesh size

    Output files (in out_dir):
        <label>.step                  CAD (Compound)
        <label>.msh                   Netgen output, Gmsh v4 format
        <label>_post.msh              same mesh + named physical groups
                                       + per-region 'region_id' ElementData view
        <label>.json                  run record
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = [n for (_p, n) in regions]

    record = {"label": label, "status": "ok", "mode": "multi_region",
              "region_names": names, "stages": []}
    try:
        s1 = _stage_cad_multi(regions, out_dir, label)
        record["stages"].append(s1)
        s2 = _stage_mesh_multi(Path(s1["step"]), names, out_dir, label,
                               maxh=maxh)
        record["stages"].append(s2)
        s3 = _stage_post_multi(Path(s2["msh"]), names, out_dir, label)
        record["stages"].append(s3)
    except Exception:
        record["status"] = "error"
        record["error"] = traceback.format_exc()
    return record
