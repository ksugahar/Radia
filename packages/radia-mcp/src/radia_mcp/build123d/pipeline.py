"""
Reusable build123d (OCC) -> Netgen (tet) -> Gmsh post pipeline.

Lab policy (2026-04-19):
    - CAD:  build123d (Python / OCCT) — main
    - Mesh: Netgen (tet) via netgen.occ — tet-only main path
    - Post: write GMSH .msh data blocks for visualization, no mesh generation
    - FreeCAD / ParaView / VTK: not used

The pipeline is designed to be called repeatedly from a sweep script so
the lab MCP servers (build123d / gmsh) get exercised on many inputs.

Stages
------
1. build123d Part   -> .brep                 (export_brep)
2. .brep            -> netgen.occ.OCCGeometry -> generate tet mesh
3. netgen Mesh      -> .msh (Gmsh2 Format)    (ngsolve Mesh.Export)
4. .msh writer      -> attaches a dummy scalar/region field
                       (proxy for post-processing from a solver)

Each stage returns a small dict so the sweep script can aggregate
pass/fail + geometry/mesh stats per run.
"""

from __future__ import annotations

import json
import time
import traceback
from pathlib import Path


__all__ = ["run_pipeline", "run_pipeline_multi", "save_record"]


def _quote_msh_string(text: str) -> str:
    return '"' + str(text).replace("\\", "\\\\").replace('"', '\\"') + '"'


def _write_gmsh_launch_companion(msh_path: Path) -> dict:
    """Write case.geo, case.geo.opt, and case.msh.opt for Gmsh display.

    Gmsh is a post-processing viewer in this pipeline. The raw .msh stores
    mesh/data; the .geo is the user-facing launch recipe that merges it and
    carries stable display defaults.
    """
    msh_path = Path(msh_path)
    geo_path = msh_path.with_suffix(".geo")
    geo_opt_path = Path(str(geo_path) + ".opt")
    msh_opt_path = Path(str(msh_path) + ".opt")
    options = "\n".join([
        "General.Orthographic = 1;",
        "General.Trackball = 0;",
        "General.RotationX = -68;",
        "General.RotationY = 0;",
        "General.RotationZ = 0;",
        "General.RotationCenterGravity = 1;",
        "Mesh.NumSubEdges = 4;",
        "Mesh.SurfaceFaces = 1;",
        "Mesh.SurfaceEdges = 1;",
        "Mesh.VolumeEdges = 0;",
        "Mesh.VolumeFaces = 0;",
        "Mesh.ColorCarousel = 2;",
        "View[0].Visible = 1;",
        "View[0].IntervalsType = 2;",
        "View[0].ShowScale = 1;",
    ]) + "\n"

    geo_path.write_text(
        "// Auto-generated build123d/Netgen Gmsh launch companion\n"
        "// Open this .geo for normal review; open .msh only for raw inspection.\n"
        f"Merge \"{msh_path.name}\";\n\n"
        + options,
        encoding="utf-8",
    )
    geo_opt_path.write_text(
        "// Auto-generated display options for the .geo launch artifact\n"
        + options,
        encoding="utf-8",
    )
    msh_opt_path.write_text(
        "// Auto-generated raw mesh/data inspection options\n" + options,
        encoding="utf-8",
    )
    return {
        "geo": str(geo_path),
        "geo_opt": str(geo_opt_path),
        "msh_opt": str(msh_opt_path),
    }


def _read_gmsh_nodes_elements(
    msh_path: Path,
) -> tuple[list[tuple[int, float, float, float]], list[dict]]:
    """Read the small ASCII GMSH subset emitted by Netgen.

    Netgen may emit GMSH v4.1, v2.2 (`$Nodes`/`$Elements`), or the older
    `$NOD`/`$ELM` form.  These carry enough information for visualization
    post data: node tags, coordinates, element type, physical/entity tags,
    and node connectivity.
    """
    lines = Path(msh_path).read_text(errors="replace").splitlines()

    if "$Nodes" in lines:
        i = lines.index("$Nodes") + 1
        node_header = lines[i].split()
        if len(node_header) == 4:
            entity_phys: dict[tuple[int, int], list[int]] = {}
            if "$Entities" in lines:
                k = lines.index("$Entities") + 1
                n_point, n_curve, n_surface, n_volume = [
                    int(v) for v in lines[k].split()[:4]
                ]
                k += 1
                for dim, count in enumerate((n_point, n_curve, n_surface, n_volume)):
                    for _ in range(count):
                        parts = lines[k].split()
                        k += 1
                        tag = int(parts[0])
                        nphys_idx = 4 if dim == 0 else 7
                        nphys = int(parts[nphys_idx])
                        entity_phys[(dim, tag)] = [
                            int(v)
                            for v in parts[nphys_idx + 1:nphys_idx + 1 + nphys]
                        ]

            n_blocks = int(node_header[0])
            nodes = []
            i += 1
            for _ in range(n_blocks):
                entity_dim, _entity_tag, parametric, n_block = [
                    int(v) for v in lines[i].split()[:4]
                ]
                i += 1
                tags = [int(lines[i + j].split()[0]) for j in range(n_block)]
                i += n_block
                for tag in tags:
                    parts = lines[i].split()
                    i += 1
                    # Parametric coordinates, when present, follow x y z.
                    if len(parts) < 3:
                        raise ValueError(
                            f"Malformed GMSH v4 node coordinates in {msh_path}"
                        )
                    nodes.append((
                        tag,
                        float(parts[0]),
                        float(parts[1]),
                        float(parts[2]),
                    ))

            j = lines.index("$Elements") + 1
            elem_header = lines[j].split()
            if len(elem_header) != 4:
                raise ValueError(f"Malformed GMSH v4 element header: {msh_path}")
            n_elem_blocks = int(elem_header[0])
            elements = []
            j += 1
            for _ in range(n_elem_blocks):
                entity_dim, entity_tag, elem_type, n_block = [
                    int(v) for v in lines[j].split()[:4]
                ]
                j += 1
                phys = entity_phys.get((entity_dim, entity_tag), [])
                physical_tag = phys[0] if phys else entity_tag
                tags = [physical_tag, entity_tag]
                for _elem in range(n_block):
                    parts = lines[j].split()
                    j += 1
                    elements.append({
                        "id": int(parts[0]),
                        "type": elem_type,
                        "tags": tags,
                        "nodes": [int(v) for v in parts[1:]],
                    })
            return nodes, elements

        n_nodes = int(lines[i])
        nodes = []
        for row in lines[i + 1:i + 1 + n_nodes]:
            parts = row.split()
            nodes.append((
                int(parts[0]),
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            ))

        j = lines.index("$Elements") + 1
        n_elements = int(lines[j])
        elements = []
        for row in lines[j + 1:j + 1 + n_elements]:
            parts = row.split()
            elem_id = int(parts[0])
            elem_type = int(parts[1])
            n_tags = int(parts[2])
            tags = [int(v) for v in parts[3:3 + n_tags]]
            node_ids = [int(v) for v in parts[3 + n_tags:]]
            elements.append({
                "id": elem_id,
                "type": elem_type,
                "tags": tags,
                "nodes": node_ids,
            })
        return nodes, elements

    if "$NOD" in lines:
        i = lines.index("$NOD") + 1
        n_nodes = int(lines[i])
        nodes = []
        for row in lines[i + 1:i + 1 + n_nodes]:
            parts = row.split()
            nodes.append((
                int(parts[0]),
                float(parts[1]),
                float(parts[2]),
                float(parts[3]),
            ))

        j = lines.index("$ELM") + 1
        n_elements = int(lines[j])
        elements = []
        for row in lines[j + 1:j + 1 + n_elements]:
            parts = row.split()
            elem_id = int(parts[0])
            elem_type = int(parts[1])
            physical_tag = int(parts[2])
            elementary_tag = int(parts[3])
            n_nodes_elem = int(parts[4])
            node_ids = [int(v) for v in parts[5:5 + n_nodes_elem]]
            elements.append({
                "id": elem_id,
                "type": elem_type,
                "tags": [physical_tag, elementary_tag],
                "nodes": node_ids,
            })
        return nodes, elements

    raise ValueError(f"Unsupported GMSH ASCII mesh format: {msh_path}")


def _write_gmsh22_with_data(
    out_msh: Path,
    nodes: list[tuple[int, float, float, float]],
    elements: list[dict],
    *,
    physical_names: list[tuple[int, int, str]] | None = None,
    node_data: tuple[str, dict[int, float]] | None = None,
    element_data: tuple[str, dict[int, float]] | None = None,
) -> None:
    """Write a GMSH v2.2 ASCII mesh plus optional data blocks."""
    out_msh = Path(out_msh)
    out_msh.parent.mkdir(parents=True, exist_ok=True)
    with out_msh.open("w", encoding="utf-8", newline="\n") as fh:
        fh.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")

        if physical_names:
            fh.write("$PhysicalNames\n")
            fh.write(f"{len(physical_names)}\n")
            for dim, tag, name in physical_names:
                fh.write(f"{dim} {tag} {_quote_msh_string(name)}\n")
            fh.write("$EndPhysicalNames\n")

        fh.write("$Nodes\n")
        fh.write(f"{len(nodes)}\n")
        for tag, x, y, z in nodes:
            fh.write(f"{tag} {x:.16g} {y:.16g} {z:.16g}\n")
        fh.write("$EndNodes\n")

        fh.write("$Elements\n")
        fh.write(f"{len(elements)}\n")
        for elem in elements:
            tags = list(elem["tags"])
            tag_text = " ".join(str(v) for v in tags)
            node_text = " ".join(str(v) for v in elem["nodes"])
            fh.write(
                f"{elem['id']} {elem['type']} {len(tags)} "
                f"{tag_text} {node_text}\n"
            )
        fh.write("$EndElements\n")

        if node_data is not None:
            name, values = node_data
            fh.write("$NodeData\n")
            fh.write("1\n")
            fh.write(f"{_quote_msh_string(name)}\n")
            fh.write("1\n0.0\n")
            fh.write("3\n0\n1\n")
            fh.write(f"{len(values)}\n")
            for tag in sorted(values):
                fh.write(f"{tag} {values[tag]:.16g}\n")
            fh.write("$EndNodeData\n")

        if element_data is not None:
            name, values = element_data
            fh.write("$ElementData\n")
            fh.write("1\n")
            fh.write(f"{_quote_msh_string(name)}\n")
            fh.write("1\n0.0\n")
            fh.write("3\n0\n1\n")
            fh.write(f"{len(values)}\n")
            for tag in sorted(values):
                fh.write(f"{tag} {values[tag]:.16g}\n")
            fh.write("$EndElementData\n")


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
    """Stage 4: attach a dummy scalar field as GMSH NodeData.

    The dummy field is f(x,y,z) = x + 2y + 3z — a trivial stand-in for
    solver output (Radia B-field, NGSolve potential, etc.).  Real
    workflows replace this with actual data using the same node-tag order.
    """
    nodes, elements = _read_gmsh_nodes_elements(msh_path)
    values = {
        tag: x + 2.0 * y + 3.0 * z
        for tag, x, y, z in nodes
    }

    out_msh = out_dir / f"{label}_post.msh"
    _write_gmsh22_with_data(
        out_msh,
        nodes,
        elements,
        node_data=(f"{label}_dummy_scalar", values),
    )
    companions = _write_gmsh_launch_companion(out_msh)

    value_list = list(values.values())
    return {
        "stage": "post",
        "ok": True,
        "msh_post": str(out_msh),
        "geo": companions["geo"],
        "geo_opt": companions["geo_opt"],
        "msh_opt": companions["msh_opt"],
        "view_tag": 0,
        "post_backend": "msh_node_data_writer",
        "n_nodes": len(nodes),
        "node_data_count": len(values),
        "dummy_field_min": min(value_list),
        "dummy_field_max": max(value_list),
        "probe_at_origin": [0.0],
        "probe_distance": None,
        "probe_note": "analytic f(0,0,0)=0; no GMSH interpolation runtime used",
    }


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
#   3. Export to Gmsh format. Per-domain tags are preserved but names do
#      not survive the exporter.
#   4. Re-write a GMSH v2.2 display file with `$PhysicalNames` and a
#      `region_id` `$ElementData` block.
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
    with TaskManager():
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
    """Stage 4 (multi): write named regions + region-id ElementData."""
    nodes, elements = _read_gmsh_nodes_elements(msh_path)

    # Tetrahedra carry the material physical tag in position 0.  Surface
    # triangles can have unrelated boundary tags, so only type 4 contributes
    # to the region-name contract.
    volume_tags = sorted({
        int(elem["tags"][0])
        for elem in elements
        if elem["type"] == 4 and elem["tags"]
    })
    if len(volume_tags) != len(region_names):
        raise RuntimeError(
            f"mesh has {len(volume_tags)} volume physical tags but "
            f"{len(region_names)} region names were given."
        )

    tag_to_region = {
        physical_tag: region_idx
        for region_idx, physical_tag in enumerate(volume_tags, start=1)
    }
    region_id_values = {
        elem["id"]: float(tag_to_region[int(elem["tags"][0])])
        for elem in elements
        if elem["type"] == 4 and elem["tags"]
    }

    out_msh = out_dir / f"{label}_post.msh"
    _write_gmsh22_with_data(
        out_msh,
        nodes,
        elements,
        physical_names=[
            (3, physical_tag, name)
            for physical_tag, name in zip(volume_tags, region_names)
        ],
        element_data=(f"{label}_region_id", region_id_values),
    )
    companions = _write_gmsh_launch_companion(out_msh)

    per_region = []
    for region_idx, (physical_tag, name) in enumerate(
        zip(volume_tags, region_names),
        start=1,
    ):
        n_elem = sum(
            1
            for elem in elements
            if elem["type"] == 4 and elem["tags"]
            and int(elem["tags"][0]) == physical_tag
        )
        per_region.append({
            "index": region_idx,
            "physical_tag": physical_tag,
            "name": name,
            "entities": [physical_tag],
            "n_elem": n_elem,
        })

    return {
        "stage": "post",
        "ok": True,
        "msh_post": str(out_msh),
        "geo": companions["geo"],
        "geo_opt": companions["geo_opt"],
        "msh_opt": companions["msh_opt"],
        "view_tag": 0,
        "post_backend": "msh_element_data_writer",
        "n_regions": len(volume_tags),
        "element_data_count": len(region_id_values),
        "regions": per_region,
    }


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
