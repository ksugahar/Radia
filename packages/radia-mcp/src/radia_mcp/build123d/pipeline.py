"""build123d -> Netgen -> checked .vol, with derived MSH 4.2 post output.

The analysis mesh is always the native Netgen .vol plus its versioned label
contract and check-vol report. Gmsh files are output only. No MSH serializer
or mesh-to-solver import is implemented here: the existing Radia exporter
owns result serialization. The standalone check-vol API needs no Cubit
installation or licence. Optional dependencies are imported only on use.

Coordinates and maxh must be in meters for analysis. No unit conversion or
physical boundary condition is inferred; the neutral mesh_boundary group
does not identify electrical ports or a physical excitation.
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

__all__ = ["run_pipeline", "run_pipeline_multi", "save_record"]

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

def _load_checked_vol(vol_path: Path):
    """Load only a checked SI solver mesh; never infer a mesh format."""
    vol_path = Path(vol_path)
    if vol_path.suffix.lower() != ".vol":
        raise ValueError("build123d analysis requires Netgen .vol; MSH is post-output only")
    contract_path = vol_path.with_suffix(".labels.json")
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("units") != "m":
        raise ValueError("Netgen solver label contract must declare SI units 'm'")
    from cubit_mesh_export.check import check_consistency
    from ngsolve import Mesh

    report = check_consistency(vol_path, contract=contract)
    check_path = vol_path.with_suffix(".vol-check.json")
    check_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not report["passed"]:
        raise ValueError(f"check-vol rejected solver input: {report['warnings']}")
    return Mesh(str(vol_path)), check_path


def _save_solver_mesh(ng_mesh, out_dir: Path, label: str, region_names=None) -> dict:
    """Save Netgen's authoritative mesh and run the canonical standalone gate."""
    from ngsolve import Mesh

    mesh = Mesh(ng_mesh)
    if mesh.dim != 3:
        raise ValueError("The build123d tet pipeline requires a three-dimensional mesh")
    if region_names is None:
        if len(mesh.GetMaterials()) != 1:
            raise ValueError("Multiple materials require the explicitly named multi-region pipeline")
        ng_mesh.SetMaterial(1, "domain")
    elif (len(set(region_names)) != len(region_names)
          or sorted(mesh.GetMaterials()) != sorted(region_names)):
        raise ValueError("Netgen materials do not match the declared distinct region names")
    # A neutral group for unnamed mesh faces; never infer a physical condition.
    for index, name in enumerate(mesh.GetBoundaries()):
        if name in ("", "default", "boundary"):
            ng_mesh.SetBCName(index, "mesh_boundary")
    vol_path = Path(out_dir) / f"{label}.vol"
    ng_mesh.Save(str(vol_path))
    contract_path = vol_path.with_suffix(".labels.json")
    contract = {
        "schema": "radia.vol-label-contract.v1", "application": "build123d_netgen",
        "strict_labels": True, "units": "m",
        "required": {"materials": sorted(set(mesh.GetMaterials())),
                     "boundaries": sorted(set(mesh.GetBoundaries()))},
    }
    contract_path.write_text(json.dumps(contract, indent=2) + "\n", encoding="utf-8")
    checked, check_path = _load_checked_vol(vol_path)
    return {"vol": str(vol_path), "label_contract": str(contract_path),
            "vol_check": str(check_path), "materials": list(checked.GetMaterials()),
            "nv": checked.nv, "ne": checked.ne, "nface": checked.nface,
            "nedge": checked.nedge}


def _validate_post_output(path: Path) -> dict:
    """Inspect derived results without passing them into a solver."""
    from radia_mcp.gmsh.msh_inspect import inspect_msh, validate_msh

    report = validate_msh(path)
    if not report["ok"]:
        raise ValueError(f"Post exporter must produce valid ASCII MSH 4.2: {report}")
    return inspect_msh(path)


def _post_exporter(mesh):
    """Require the current typed result writer before creating any MSH output."""
    try:
        from radia.gmsh_post_export import GmshPostExport, PayloadV42
    except ImportError as exc:
        raise RuntimeError("build123d post output requires the current Radia MSH 4.2 exporter") from exc
    # The typed payload is the writer's format capability, independent of
    # the native solver ABI or an unchanged development package version.
    if not isinstance(PayloadV42, type):
        raise RuntimeError("Radia MSH 4.2 payload capability is unavailable")
    return GmshPostExport(mesh)


def _stage_mesh(brep_path: Path, out_dir: Path, label: str, maxh: float = 5.0) -> dict:
    """Generate with Netgen and save checked native .vol, never an MSH intermediate."""
    from netgen.occ import OCCGeometry
    from ngsolve import TaskManager

    start = time.perf_counter()
    geometry = OCCGeometry(str(brep_path))
    with TaskManager():
        ng_mesh = geometry.GenerateMesh(maxh=maxh)
    seconds = time.perf_counter() - start
    return {"stage": "mesh", "ok": True, "gen_seconds": round(seconds, 3),
            "maxh": maxh, **_save_solver_mesh(ng_mesh, out_dir, label)}


def _stage_mesh_multi(step_path: Path, region_names, out_dir: Path,
                      label: str, maxh: float = 5.0) -> dict:
    """Keep named Netgen solids and use .vol for the analysis hand-off."""
    from netgen.occ import OCCGeometry, Glue
    from ngsolve import TaskManager

    start = time.perf_counter()
    geometry = OCCGeometry(str(step_path))
    solids = list(geometry.shape.solids)
    if len(solids) != len(region_names):
        raise ValueError(f"STEP has {len(solids)} solids but {len(region_names)} region names")
    for solid, name in zip(solids, region_names):
        solid.mat(name)
    geometry = OCCGeometry(Glue(solids))
    with TaskManager():
        ng_mesh = geometry.GenerateMesh(maxh=maxh)
    seconds = time.perf_counter() - start
    return {"stage": "mesh", "ok": True, "gen_seconds": round(seconds, 3),
            "maxh": maxh, **_save_solver_mesh(ng_mesh, out_dir, label, region_names)}


def _stage_post(vol_path: Path, out_dir: Path, label: str) -> dict:
    """Attach the existing analytic demo field to a derived MSH 4.2 result."""
    mesh, check_path = _load_checked_vol(vol_path)
    from ngsolve import x, y, z

    output = Path(out_dir) / f"{label}_post.msh"
    post = _post_exporter(mesh)
    post.add_scalar_field(f"{label}_dummy_scalar", x + 2.0*y + 3.0*z)
    post.write(str(output))
    summary = _validate_post_output(output)
    companions = _write_gmsh_launch_companion(output)
    values = [sum(a*b for a,b in zip(vertex.point, (1.,2.,3.)))
              for vertex in mesh.vertices]
    return {"stage": "post", "ok": True, "msh_post": str(output),
            **companions, "view_tag": 0, "post_backend": "radia_msh42_exporter",
            "n_nodes": summary["nodes"]["count"],
            "node_data_count": summary["nodes"]["count"],
            "dummy_field_min": min(values), "dummy_field_max": max(values),
            "probe_at_origin": [0.0], "probe_distance": None,
            "probe_note": "analytic f(0,0,0)=0; no native interpolation used",
            "solver_vol": str(vol_path), "vol_check": str(check_path)}


def _stage_post_multi(vol_path: Path, region_names, out_dir: Path, label: str) -> dict:
    """Export region display data from checked .vol; never read MSH as analysis input."""
    mesh, check_path = _load_checked_vol(vol_path)
    from ngsolve import VOL

    materials = list(mesh.GetMaterials())
    if len(set(region_names)) != len(region_names) or sorted(materials) != sorted(region_names):
        raise ValueError("Checked .vol materials do not match the declared distinct regions")
    output = Path(out_dir) / f"{label}_post.msh"
    post = _post_exporter(mesh)
    post.add_scalar_field(f"{label}_region_id",
                          mesh.MaterialCF({name: float(i) for i,name in enumerate(region_names,1)}),
                          cell_data=True)
    post.write(str(output))
    summary = _validate_post_output(output)
    physical = {row["name"]: row["tag"] for row in summary["physical_names"] if row["dim"] == 3}
    if set(physical) != set(region_names):
        raise ValueError("Post output physical names do not match the checked .vol regions")
    companions = _write_gmsh_launch_companion(output)
    regions = [{"index": i, "physical_tag": physical[name], "name": name,
                "entities": [physical[name]],
                "n_elem": sum(element.mat == name for element in mesh.Elements(VOL))}
               for i,name in enumerate(region_names,1)]
    return {"stage": "post", "ok": True, "msh_post": str(output), **companions,
            "view_tag": 0, "post_backend": "radia_msh42_exporter",
            "n_regions": len(regions), "element_data_count": mesh.ne,
            "regions": regions, "solver_vol": str(vol_path), "vol_check": str(check_path)}


def run_pipeline(part, out_dir, label: str, maxh: float = 5.0) -> dict:
    """Build CAD, checked Netgen .vol, and a derived MSH 4.2 scalar display.

    part coordinates and maxh are caller-owned SI lengths. No Cubit runtime
    is used. check-vol, NGSolve and Radia's exporter are required on this path.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record = {"label": label, "status": "ok", "stages": []}
    try:
        cad = _stage_cad(part, out_dir, label)
        record["stages"].append(cad)
        mesh = _stage_mesh(Path(cad["brep"]), out_dir, label, maxh=maxh)
        record["stages"].append(mesh)
        record["stages"].append(_stage_post(Path(mesh["vol"]), out_dir, label))
    except Exception:
        record["status"] = "error"
        record["error"] = traceback.format_exc()
    return record


def run_pipeline_multi(regions, out_dir, label: str, maxh: float = 5.0) -> dict:
    """Build named Netgen .vol regions and derived region-id post output.

    Outputs: CAD .step, authoritative .vol, .labels.json, .vol-check.json,
    and derived _post.msh with exact display launch companions. The i-th
    declared region owns region_id i regardless of the mesh's material order.
    Coordinates/maxh must be meters; no physical condition is inferred.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = [name for _part,name in regions]
    record = {"label": label, "status": "ok", "mode": "multi_region",
              "region_names": names, "stages": []}
    try:
        cad = _stage_cad_multi(regions, out_dir, label)
        record["stages"].append(cad)
        mesh = _stage_mesh_multi(Path(cad["step"]), names, out_dir, label, maxh=maxh)
        record["stages"].append(mesh)
        record["stages"].append(_stage_post_multi(Path(mesh["vol"]), names, out_dir, label))
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
