"""build123d -> Netgen -> checked .vol, with derived MSH 4.2 post output.

The analysis mesh is always the native Netgen .vol plus its versioned label
contract and check-vol report. Gmsh files are output only. No MSH serializer
or mesh-to-solver import is implemented here: the existing Radia exporter
owns result serialization. The standalone check-vol API needs no Cubit
installation or licence. Optional dependencies are imported only on use.

Callers must explicitly declare whether CAD coordinates and maxh are in
meters or millimeters. Geometry and maxh are converted together before meshing
so the authoritative .vol and embedded CAD are SI. No physical boundary
condition is inferred; the neutral mesh_boundary group
does not identify electrical ports or a physical excitation.
"""
from __future__ import annotations

import json
import time
import traceback
from pathlib import Path

__all__ = ["run_pipeline", "run_pipeline_multi", "save_record"]


_DEPENDENCY_HELP = (
    "Install the radia-mcp build123d extra. The canonical "
    "cubit-mesh-export.check checker is available only on Windows/Python 3.12; "
    "the analysis path requires that platform and NGSolve. Post output also needs a Radia "
    "release newer than 5.3.3 exposing GmshPostExport and PayloadV42. "
    "The published 5.3.3 writer is unsupported; a current source serializer "
    "can be validated separately until the post-migration Radia release is published."
)


def _length_scale_to_m(units: str) -> float:
    """Interpret an explicit unit declaration; never infer it from dimensions."""
    if not isinstance(units, str) or units not in ("m", "mm"):
        raise ValueError("units must be explicitly declared as 'm' or 'mm'")
    return 1.0 if units == "m" else 1.0e-3

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

def _stage_cad_multi(regions, out_dir: Path, label: str, *, units: str) -> dict:
    """Stage 1 (multi): write one STEP carrying explicitly unit-labelled solids."""
    _length_scale_to_m(units)
    from build123d import Compound, Unit, export_step

    parts = [p for (p, _n) in regions]
    names = [n for (_p, n) in regions]

    # Compound preserves the ordering of children, which is our region-id contract.
    compound = Compound(children=parts)
    compound.label = label

    step = out_dir / f"{label}.step"
    # build123d declares the shape units in its STEP document. OCCT STEP
    # import then normalizes coordinates to millimeters, independent of input units.
    export_step(compound, str(step), unit=Unit.M if units == "m" else Unit.MM)

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
    try:
        from cubit_mesh_export.check import check_consistency
        from ngsolve import Mesh
    except ImportError as exc:
        raise RuntimeError(_DEPENDENCY_HELP) from exc

    report = check_consistency(vol_path, contract=contract)
    check_path = vol_path.with_suffix(".vol-check.json")
    check_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if not report["passed"]:
        raise ValueError(f"check-vol rejected solver input: {report['warnings']}")
    return Mesh(str(vol_path)), check_path


def _save_solver_mesh(ng_mesh, out_dir: Path, label: str, region_names=None, *, units: str) -> dict:
    """Save an explicitly SI Netgen mesh and run the canonical standalone gate."""
    if _length_scale_to_m(units) != 1.0:
        raise ValueError("Solver mesh must already use SI units 'm'; convert its CAD before meshing")
    try:
        from ngsolve import Mesh
    except ImportError as exc:
        raise RuntimeError(_DEPENDENCY_HELP) from exc

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
    # Netgen boundary names retain CR on Linux; keep this native artifact portable.
    vol_path.write_bytes(vol_path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b""))
    contract_path = vol_path.with_suffix(".labels.json")
    contract = {
        "schema": "radia.vol-label-contract.v1", "application": "build123d_netgen",
        "strict_labels": True, "units": units,
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
        raise RuntimeError("build123d post output requires the current Radia MSH 4.2 exporter (PayloadV42). " + _DEPENDENCY_HELP) from exc
    # The typed payload is the writer's format capability, independent of
    # the native solver ABI or an unchanged development package version.
    if not isinstance(PayloadV42, type):
        raise RuntimeError("Radia PayloadV42 capability is unavailable. " + _DEPENDENCY_HELP)
    return GmshPostExport(mesh)


def _stage_mesh(brep_path: Path, out_dir: Path, label: str, maxh: float = 5.0,
                *, units: str) -> dict:
    """Interpret unitless BREP coordinates explicitly and mesh SI geometry."""
    scale = _length_scale_to_m(units)
    try:
        from netgen.occ import OCCGeometry, Pnt
        from ngsolve import TaskManager
    except ImportError as exc:
        raise RuntimeError(_DEPENDENCY_HELP) from exc

    start = time.perf_counter()
    geometry = OCCGeometry(str(brep_path))
    if scale != 1.0:
        geometry = OCCGeometry(geometry.shape.Scale(Pnt(0, 0, 0), scale))
    # Scale CAD before generation, rather than scaling mesh points afterwards:
    # the latter would leave the .vol's embedded CAD in the old units.
    maxh_m = maxh * scale
    with TaskManager():
        ng_mesh = geometry.GenerateMesh(maxh=maxh_m)
    seconds = time.perf_counter() - start
    return {"stage": "mesh", "ok": True, "gen_seconds": round(seconds, 3),
            "maxh": maxh, "input_units": units, "length_scale_to_m": scale,
            "maxh_m": maxh_m,
            **_save_solver_mesh(ng_mesh, out_dir, label, units="m")}


def _stage_mesh_multi(step_path: Path, region_names, out_dir: Path,
                      label: str, maxh: float = 5.0, *, units: str) -> dict:
    """Read STEP's declared units, then preserve named SI Netgen solids."""
    input_scale = _length_scale_to_m(units)
    try:
        from netgen.occ import OCCGeometry, Glue, Pnt
        from ngsolve import TaskManager
    except ImportError as exc:
        raise RuntimeError(_DEPENDENCY_HELP) from exc

    start = time.perf_counter()
    # OCCT STEP import normalizes document units to millimeters. The CAD stage
    # declares its source shape units via export_step(unit=...), so this reader
    # conversion is fixed while maxh still follows the caller's declared units.
    geometry = OCCGeometry(str(step_path))
    geometry = OCCGeometry(geometry.shape.Scale(Pnt(0, 0, 0), 1.0e-3))
    solids = list(geometry.shape.solids)
    if len(solids) != len(region_names):
        raise ValueError(f"STEP has {len(solids)} solids but {len(region_names)} region names")
    for solid, name in zip(solids, region_names):
        solid.mat(name)
    geometry = OCCGeometry(Glue(solids))
    maxh_m = maxh * input_scale
    with TaskManager():
        ng_mesh = geometry.GenerateMesh(maxh=maxh_m)
    seconds = time.perf_counter() - start
    return {"stage": "mesh", "ok": True, "gen_seconds": round(seconds, 3),
            "maxh": maxh, "input_units": units, "length_scale_to_m": input_scale,
            "cad_reader_units": "mm", "cad_reader_scale_to_m": 1.0e-3,
            "maxh_m": maxh_m,
            **_save_solver_mesh(ng_mesh, out_dir, label, region_names, units="m")}


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


def run_pipeline(part, out_dir, label: str, maxh: float = 5.0, *, units: str) -> dict:
    """Build CAD, checked Netgen .vol, and a derived MSH 4.2 scalar display.

    ``units`` is required: CAD coordinates and maxh must share 'm' or 'mm'.
    CAD statistics keep those input units; solver/post artifacts are SI.
    check-vol, NGSolve and the current Radia exporter are required. No Cubit
    runtime or unit inference is used.
    """
    scale = _length_scale_to_m(units)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    record = {"label": label, "status": "ok", "input_units": units,
              "length_scale_to_m": scale, "stages": []}
    try:
        cad = _stage_cad(part, out_dir, label)
        cad.update(input_units=units, volume_units=f"{units}^3", area_units=f"{units}^2")
        record["stages"].append(cad)
        mesh = _stage_mesh(Path(cad["brep"]), out_dir, label, maxh=maxh, units=units)
        record["stages"].append(mesh)
        record["stages"].append(_stage_post(Path(mesh["vol"]), out_dir, label))
    except Exception:
        record["status"] = "error"
        record["error"] = traceback.format_exc()
    return record


def run_pipeline_multi(regions, out_dir, label: str, maxh: float = 5.0, *, units: str) -> dict:
    """Build named Netgen .vol regions and derived region-id post output.

    Outputs: CAD .step, authoritative .vol, .labels.json, .vol-check.json,
    and derived _post.msh with exact display launch companions. The i-th
    declared region owns region_id i regardless of the mesh's material order.
    Required ``units`` declares CAD coordinates/maxh as 'm' or 'mm'.
    CAD statistics retain input units; checked .vol and post output are SI.
    No physical condition or unit is inferred.
    """
    scale = _length_scale_to_m(units)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    names = [name for _part,name in regions]
    record = {"label": label, "status": "ok", "mode": "multi_region",
              "region_names": names, "input_units": units,
              "length_scale_to_m": scale, "stages": []}
    try:
        cad = _stage_cad_multi(regions, out_dir, label, units=units)
        cad.update(input_units=units, volume_units=f"{units}^3", area_units=f"{units}^2")
        record["stages"].append(cad)
        mesh = _stage_mesh_multi(Path(cad["step"]), names, out_dir, label, maxh=maxh, units=units)
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
