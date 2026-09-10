"""Build the canonical Cubit 2025.12 meshes for the C-type comparison.

Two meshes are deliberate:

* ``iron.vol`` contains only the exact ACIS iron and is the HDiv-MMM open-
  boundary input. HDiv-MMM must not acquire an artificial air domain.
* ``kelvin_domain.vol`` contains the same iron, a spherical physical-air
  domain, and the translated Kelvin exterior sphere. Reduced-A and
  Omega-reduced-Omega share this open-boundary mesh.

Both artifacts originate from ``cad/c_type_iron.jou``. Python never authors a
replacement C-yoke with ``netgen.occ``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import time
from pathlib import Path


HERE = Path(__file__).resolve().parent
CAD_JOURNAL = HERE / "cad" / "c_type_iron.jou"
CUBIT_DEFAULT = Path(
    r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.com"
)
EXACT_IRON_VOLUME_M3 = 1_446_095.333333e-9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run(command: list[str], cwd: Path, log: Path) -> dict[str, object]:
    started = time.perf_counter()
    completed = subprocess.run(
        command,
        cwd=cwd,
        text=True,
        encoding="utf-8",
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log.write_text(completed.stdout, encoding="utf-8")
    return {
        "command": command,
        "returncode": int(completed.returncode),
        "runtime_s": float(time.perf_counter() - started),
        "log": str(log),
    }


def _write_journal(
    path: Path,
    output: Path,
    *,
    iron_size_m: float,
    air_size_m: float | None,
    gap_size_m: float,
    kelvin_radius_m: float,
    kelvin_mesh_size_m: float,
    curve_order: int,
) -> None:
    cad_lines = CAD_JOURNAL.read_text(encoding="utf-8").rstrip().splitlines()
    reflect_command = "volume all copy reflect z"
    if not cad_lines or cad_lines[-1].strip().lower() != reflect_command:
        raise RuntimeError(
            f"{CAD_JOURNAL} must end with {reflect_command!r}; the canonical "
            "builder replaces that geometry-only reflection with a reflected "
            "mesh copy"
        )
    positive_cad = "\n".join(cad_lines[:-1]).rstrip()
    physical_helper = HERE / "cubit_reflection_mesh.py"
    physical_script = path.with_suffix(".physical.py")
    physical_script.write_text(
        "import importlib.util\n"
        f'_spec = importlib.util.spec_from_file_location("_radia_c_type_reflection", r"{physical_helper}")\n'
        "_module = importlib.util.module_from_spec(_spec)\n"
        "_spec.loader.exec_module(_module)\n"
        "_physical = _module.build_reflection_invariant_physical_mesh("
        f"iron_size={float(iron_size_m):.17g}, "
        f"air_size={None if air_size_m is None else float(air_size_m)!r}, "
        f"gap_size={float(gap_size_m):.17g}, "
        f"kelvin_radius={float(kelvin_radius_m):.17g})\n"
        "assert _physical is not None\n",
        encoding="utf-8",
    )
    lines = [positive_cad, "", f'play "{physical_script.as_posix()}"']
    if air_size_m is None:
        export = f'export netgen "{output.as_posix()}" order 1 overwrite'
    else:
        helper = (
            HERE.parents[1]
            / "packages"
            / "cubit-mesh-export"
            / "src"
            / "cubit_mesh_export"
            / "cubit_helpers"
            / "add_kelvin.py"
        )
        kelvin_script = path.with_suffix(".kelvin.py")
        kelvin_script.write_text(
            "import importlib.util\n"
            f'_spec = importlib.util.spec_from_file_location("_radia_c_type_kelvin", r"{helper}")\n'
            "_module = importlib.util.module_from_spec(_spec)\n"
            "_spec.loader.exec_module(_module)\n"
            f'_info = _module.add_kelvin_cubit(R={float(kelvin_radius_m):.17g}, air_block="air", symmetry=["z"], mesh_size={float(kelvin_mesh_size_m):.17g}, kelvin_block="kelvin")\n'
            "assert _info is not None\n",
            encoding="utf-8",
        )
        lines.append(f'play "{kelvin_script.as_posix()}"')
        export = f'export netgen "{output.as_posix()}" order {int(curve_order)} overwrite'
    lines.extend([export, "exit"])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _check_vol(vol: Path, output_dir: Path, *, kelvin_domain: bool) -> dict[str, object]:
    report = output_dir / f"{vol.stem}.vol-check.json"
    command = [
        "check-vol",
        str(vol),
        "--strict-labels",
        "--required-materials",
        "iron,air,kelvin" if kelvin_domain else "iron",
        "--required-boundaries",
        "iron_air_interface,kelvin_int,kelvin_ext" if kelvin_domain else "iron_boundary",
    ]
    if kelvin_domain:
        command.extend(
            [
                "--air-materials",
                "air,kelvin",
            ]
        )
    command.extend(["--report-json", str(report), "--format", "json"])
    run = _run(command, output_dir, output_dir / f"{vol.stem}.check-vol.log")
    if not report.is_file():
        raise RuntimeError(f"check-vol did not create {report}")
    payload = json.loads(report.read_text(encoding="utf-8"))
    if not payload.get("passed", False):
        raise RuntimeError(f"check-vol rejected {vol}; see {report}")
    return {"run": run, "report": str(report), "payload": payload}


def _ngsolve_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng

    mesh = ng.Mesh(str(path))
    vertex_counts = sorted({len(element.vertices) for element in mesh.Elements(ng.VOL)})
    return {
        "elements": int(mesh.ne),
        "vertices": int(mesh.nv),
        "materials": list(mesh.GetMaterials()),
        "boundaries": sorted(set(mesh.GetBoundaries())),
        "volume_element_vertex_counts": vertex_counts,
    }


def _reflection_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    coordinates = np.asarray([vertex.point for vertex in mesh.vertices], dtype=float)
    vertex_keys = {tuple(np.round(point, 14)) for point in coordinates}
    missing_vertices = 0
    maximum_vertex_error = 0.0
    for point in coordinates:
        reflected = np.asarray((point[0], point[1], -point[2]), dtype=float)
        key = tuple(np.round(reflected, 14))
        if key not in vertex_keys:
            missing_vertices += 1
            continue
        candidates = coordinates[
            np.all(np.isclose(coordinates, reflected, rtol=0.0, atol=1e-13), axis=1)
        ]
        if candidates.size:
            maximum_vertex_error = max(
                maximum_vertex_error,
                float(np.min(np.linalg.norm(candidates - reflected, axis=1))),
            )

    def element_signature(element, *, reflect: bool) -> tuple[object, ...]:
        points = []
        for vertex in element.vertices:
            point = np.asarray(mesh.vertices[vertex.nr].point, dtype=float)
            if reflect:
                point[2] = -point[2]
            points.append(tuple(np.round(point, 14)))
        return str(element.mat), tuple(sorted(points))

    elements = tuple(mesh.Elements(ng.VOL))
    signatures = {element_signature(element, reflect=False) for element in elements}
    missing_elements = sum(
        element_signature(element, reflect=True) not in signatures
        for element in elements
    )
    return {
        "vertex_count": int(len(coordinates)),
        "missing_reflected_vertices": int(missing_vertices),
        "maximum_reflected_vertex_error_m": maximum_vertex_error,
        "volume_element_count": int(len(elements)),
        "missing_reflected_volume_elements": int(missing_elements),
    }


def _kelvin_identification_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    pairs = mesh.ngmesh.GetIdentifications()
    points = mesh.ngmesh.Points()
    displacements = np.asarray(
        [
            np.asarray(points[slave].p, dtype=float)
            - np.asarray(points[master].p, dtype=float)
            for master, slave in pairs
        ],
        dtype=float,
    )
    if displacements.size == 0:
        return {
            "pair_count": 0,
            "translation_m": [0.0, 0.0, 0.0],
            "maximum_pair_translation_error_m": None,
        }
    translation = np.mean(displacements, axis=0)
    return {
        "pair_count": int(len(pairs)),
        "translation_m": translation.tolist(),
        "maximum_pair_translation_error_m": float(
            np.max(np.linalg.norm(displacements - translation, axis=1))
        ),
    }


def _kelvin_fes_inventory(path: Path) -> dict[str, object]:
    """Verify that the persisted identification actually couples H1 traces."""
    import ngsolve as ng

    mesh = ng.Mesh(str(path))
    base = ng.H1(mesh, order=1, dirichlet_bbbnd="GND")
    periodic = ng.Periodic(base)
    slaved_dofs = int(sum(base.FreeDofs()) - sum(periodic.FreeDofs()))

    trace = ng.GridFunction(periodic)
    trace.vec[:] = 0.0
    trace.Set(1.0, definedon=mesh.Boundaries("kelvin_int"))
    inner_norm = float(
        ng.Integrate(
            trace * trace,
            mesh,
            definedon=mesh.Boundaries("kelvin_int"),
        )
    )
    outer_norm = float(
        ng.Integrate(
            trace * trace,
            mesh,
            definedon=mesh.Boundaries("kelvin_ext"),
        )
    )
    ratio = outer_norm / inner_norm if inner_norm > 0.0 else None
    return {
        "base_h1_ndof": int(base.ndof),
        "periodic_h1_ndof": int(periodic.ndof),
        "slaved_free_dofs": slaved_dofs,
        "kelvin_int_trace_norm": inner_norm,
        "kelvin_ext_trace_norm": outer_norm,
        "trace_norm_ratio": ratio,
    }


def _gap_inventory(path: Path) -> dict[str, object]:
    import ngsolve as ng
    import numpy as np

    mesh = ng.Mesh(str(path))
    elements = 0
    maximum_edge = 0.0
    maximum_z_span = 0.0
    for element in mesh.Elements(ng.VOL):
        material = str(element.mat) if hasattr(element, "mat") else ""
        if material != "air":
            continue
        points = np.asarray(
            [mesh.vertices[vertex.nr].point for vertex in element.vertices],
            dtype=float,
        )
        centroid = np.mean(points, axis=0)
        if not (
            abs(float(centroid[0])) <= 0.0171
            and abs(float(centroid[1])) <= 0.0121
            and abs(float(centroid[2])) <= 0.0051
        ):
            continue
        elements += 1
        maximum_z_span = max(
            maximum_z_span, float(np.max(points[:, 2]) - np.min(points[:, 2]))
        )
        for left in range(len(points)):
            for right in range(left + 1, len(points)):
                maximum_edge = max(
                    maximum_edge,
                    float(np.linalg.norm(points[left] - points[right])),
                )
    gap_height = 0.010
    profile = _gap_line_profile(mesh, gap_height)
    return {
        "elements": elements,
        "maximum_edge_m": maximum_edge,
        "maximum_z_span_m": maximum_z_span,
        "gap_height_m": gap_height,
        # A RESOLUTION INDICATOR derived from the widest element in z, not a
        # layer count: gap_height / max_z_span.  It is reported because
        # `gap_size` does not control it -- that value is applied as a surface
        # size on the pole faces, so it refines the mesh IN PLANE while the
        # tets stay free to span the whole half-gap in z.  A family can shrink
        # `gap_size` from 2.5 mm to 1.28 mm and move this from 2.00 to 2.55.
        "z_resolution_indicator": (gap_height / maximum_z_span
                                   if maximum_z_span > 0.0 else 0.0),
        # The count that the indicator only approximates: how many distinct
        # elements a line actually crosses on its way through the gap.
        "line_profile": profile,
    }


def _curved_line_segments(mesh, x0: float, y0: float, half: float,
                          coarse_samples: int = 2001, tolerance: float = 1.0e-9):
    """Segments along a gap line on the CURVED geometry.

    The vertex-tet clipping works on the straightened cell, so it can only
    report straight-geometry segment lengths.  Agreeing element COUNTS do not
    show that the endpoints or the longest segment agree.  Here the switch
    positions themselves are located: a coarse walk finds each change of
    element id, then bisection drives the change to ``tolerance`` so the
    segment boundaries -- and therefore the lengths -- are measured on the
    geometry the solver actually integrates over.

    LIMIT: bisection only refines transitions the COARSE WALK already found.
    A segment shorter than the coarse spacing can be stepped over entirely,
    so a disagreement in the segment COUNT against the straight-geometry
    clipping is not attributable to curvature without further work -- it can
    equally be a missed short segment here.  Measured on the C-type family
    2026-09-10: the two agree on the maximum segment length to about 1 um on
    the five representative lines, while one line disagreed 14 against 13 in
    the count, cause undetermined.
    """
    import numpy as np
    import ngsolve as ng

    materials = mesh.GetMaterials()

    def identify(z):
        try:
            point = mesh(float(x0), float(y0), float(z))
        except Exception:                                       # noqa: BLE001
            return None
        number = int(point.nr)
        if str(materials[mesh[ng.ElementId(ng.VOL, number)].index]) != "air":
            return None
        return number

    grid = np.linspace(-half + 1e-9, half - 1e-9, int(coarse_samples))
    ids = [identify(z) for z in grid]
    boundaries = []
    for index in range(1, len(grid)):
        if ids[index] == ids[index - 1]:
            continue
        low, high = float(grid[index - 1]), float(grid[index])
        low_id = ids[index - 1]
        while high - low > tolerance:
            middle = 0.5 * (low + high)
            if identify(middle) == low_id:
                low = middle
            else:
                high = middle
        boundaries.append(0.5 * (low + high))
    edges = [-half] + boundaries + [half]
    segments = []
    for index in range(len(edges) - 1):
        centre = 0.5 * (edges[index] + edges[index + 1])
        if identify(centre) is None:
            continue
        segments.append(edges[index + 1] - edges[index])
    return {
        "segments": len(segments),
        "maximum_segment_m": float(max(segments)) if segments else 0.0,
        "minimum_segment_m": float(min(segments)) if segments else 0.0,
        "covered_m": float(sum(segments)),
        "switch_tolerance_m": float(tolerance),
    }


def _gap_line_profile(mesh, gap_height: float,
                      offsets=((0.0, 0.0), (0.008, 0.0), (0.0, 0.005),
                               (0.008, 0.005), (-0.008, -0.005))) -> dict:
    """Measure what a line actually meets on its way across the gap.

    A crossing COUNT is not a division count: a line can clip several tets
    over a short span and still run through one coarse element elsewhere, so
    the count alone cannot say whether a coarse stretch remains.  The segment
    the line spends inside each element is therefore computed EXACTLY, by
    clipping the line against each tetrahedron's four face half-spaces, which
    also removes the dependence on a sampling density.  Coverage is checked
    too: the union of the segments must fill the gap without a hole and
    without overlap, or the measurement is not describing one clean traversal.

    The clipping uses the straight tetrahedron through the element's vertices.
    On a curved mesh that is an approximation of the curved cell, so the
    sampled counts -- which go through NGSolve's own point location and see
    the curved geometry -- are reported next to it as a cross-check.
    """
    import numpy as np
    import ngsolve as ng

    half = 0.5 * float(gap_height)
    materials = mesh.GetMaterials()
    air_boxes = []
    for element in mesh.Elements(ng.VOL):
        if str(materials[element.index]) != "air":
            continue
        points = np.asarray(
            [mesh.vertices[v.nr].point for v in element.vertices], dtype=float)
        if points.shape[0] != 4:
            continue
        lo, hi = points.min(axis=0), points.max(axis=0)
        if hi[2] < -half or lo[2] > half:
            continue
        air_boxes.append((lo, hi, points))

    def _segments(x0, y0):
        found = []
        for lo, hi, points in air_boxes:
            if not (lo[0] <= x0 <= hi[0] and lo[1] <= y0 <= hi[1]):
                continue
            low, high = -half, half
            ok = True
            for drop in range(4):
                face = [points[i] for i in range(4) if i != drop]
                normal = np.cross(face[1] - face[0], face[2] - face[0])
                inward = points[drop] - face[0]
                if float(normal @ inward) < 0.0:
                    normal = -normal
                offset = float(normal @ face[0])
                a = float(normal[2])
                b = float(normal[0] * x0 + normal[1] * y0) - offset
                if abs(a) < 1.0e-300:
                    if b < 0.0:
                        ok = False
                        break
                    continue
                bound = -b / a
                if a > 0.0:
                    low = max(low, bound)
                else:
                    high = min(high, bound)
                if low >= high:
                    ok = False
                    break
            if ok and high - low > 1.0e-12:
                found.append((low, high))
        found.sort()
        return found

    def _sampled_count(x0, y0, samples):
        previous, count = None, 0
        for value in np.linspace(-half + 1e-9, half - 1e-9, samples):
            try:
                point = mesh(float(x0), float(y0), float(value))
            except Exception:                                   # noqa: BLE001
                previous = None
                continue
            number = int(point.nr)
            if str(materials[mesh[ng.ElementId(ng.VOL, number)].index]) != "air":
                previous = None
                continue
            if number != previous:
                count += 1
                previous = number
        return count

    lines = []
    for x0, y0 in offsets:
        found = _segments(float(x0), float(y0))
        lengths = [high - low for low, high in found]
        covered = float(sum(lengths))
        holes, overlaps = 0.0, 0.0
        cursor = -half
        for low, high in found:
            if low > cursor + 1.0e-10:
                holes += low - cursor
            elif low < cursor - 1.0e-10:
                overlaps += cursor - low
            cursor = max(cursor, high)
        if cursor < half - 1.0e-10:
            holes += half - cursor
        lines.append({
            "x_m": float(x0), "y_m": float(y0),
            "segments": len(found),
            "maximum_segment_m": float(max(lengths)) if lengths else 0.0,
            "minimum_segment_m": float(min(lengths)) if lengths else 0.0,
            "covered_m": covered,
            "uncovered_m": float(holes),
            "overlap_m": float(overlaps),
            "sampled_4001": _sampled_count(x0, y0, 4001),
            "sampled_16001": _sampled_count(x0, y0, 16001),
            "curved": _curved_line_segments(mesh, float(x0), float(y0), half),
        })
    counts = [row["segments"] for row in lines if row["segments"]]
    return {
        "lines": lines,
        "minimum_segments": int(min(counts)) if counts else 0,
        "maximum_segment_m": float(max(row["maximum_segment_m"] for row in lines))
        if lines else 0.0,
        "clean_traversal": all(row["uncovered_m"] < 1.0e-9
                               and row["overlap_m"] < 1.0e-9 for row in lines),
        "sampling_is_stable": all(row["sampled_4001"] == row["sampled_16001"]
                                  for row in lines),
        # Curved geometry, measured by bisecting the element-id switch rather
        # than inferred from an agreeing count.
        "curved_minimum_segments": int(min(row["curved"]["segments"]
                                           for row in lines)) if lines else 0,
        "curved_maximum_segment_m": float(max(row["curved"]["maximum_segment_m"]
                                              for row in lines)) if lines else 0.0,
        "geometry_of_lengths": "straight vertex tetrahedra; the 'curved' "
                               "entries repeat the measurement on the curved "
                               "cells",
        "note": "measured on the listed representative lines only; not a "
                "statement about the whole gap",
    }

def build(options: argparse.Namespace) -> dict[str, object]:
    import ngsolve as ng

    output_dir = options.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    cubit = options.cubit.resolve()
    if not cubit.is_file():
        raise FileNotFoundError(cubit)
    if not CAD_JOURNAL.is_file():
        raise FileNotFoundError(CAD_JOURNAL)
    command_plugin_dir = (
        options.command_plugin_dir.resolve()
        if options.command_plugin_dir is not None
        else None
    )
    if command_plugin_dir is not None:
        command_plugin = command_plugin_dir / "cubit_mesh_export.ccm"
        if not command_plugin.is_file():
            raise FileNotFoundError(command_plugin)

    iron_vol = output_dir / "iron.vol"
    kelvin_vol = output_dir / "kelvin_domain.vol"
    iron_journal = output_dir / "build_iron.jou"
    kelvin_journal = output_dir / "build_kelvin_domain.jou"

    _write_journal(
        iron_journal,
        iron_vol,
        iron_size_m=options.iron_size,
        air_size_m=None,
        gap_size_m=options.gap_size,
        kelvin_radius_m=options.kelvin_radius,
        kelvin_mesh_size_m=options.kelvin_mesh_size,
        curve_order=options.curve_order,
    )
    _write_journal(
        kelvin_journal,
        kelvin_vol,
        iron_size_m=options.iron_size,
        air_size_m=options.air_size,
        gap_size_m=options.gap_size,
        kelvin_radius_m=options.kelvin_radius,
        kelvin_mesh_size_m=options.kelvin_mesh_size,
        curve_order=options.curve_order,
    )

    runs = {}
    for name, journal, vol in (
        ("iron", iron_journal, iron_vol),
        ("kelvin_domain", kelvin_journal, kelvin_vol),
    ):
        command = [str(cubit), "-nographics", "-batch", "-nojournal"]
        if command_plugin_dir is not None:
            command.extend(["-commandplugindir", str(command_plugin_dir)])
        command.extend(["-input", str(journal)])
        run = _run(
            command,
            output_dir,
            output_dir / f"{name}.cubit.log",
        )
        runs[name] = run
        if not vol.is_file():
            raise RuntimeError(f"Cubit did not create {vol}; see {run['log']}")

    checks = {
        "iron": _check_vol(iron_vol, output_dir, kelvin_domain=False),
        "kelvin_domain": _check_vol(kelvin_vol, output_dir, kelvin_domain=True),
    }
    with ng.TaskManager():
        inventories = {
            "iron": _ngsolve_inventory(iron_vol),
            "kelvin_domain": _ngsolve_inventory(kelvin_vol),
        }
        reflection = {
            "iron": _reflection_inventory(iron_vol),
            "kelvin_domain": _reflection_inventory(kelvin_vol),
        }
        kelvin_identification = _kelvin_identification_inventory(kelvin_vol)
        kelvin_fes = _kelvin_fes_inventory(kelvin_vol)
        gap_inventory = _gap_inventory(kelvin_vol)
    for name, inventory in inventories.items():
        if inventory["volume_element_vertex_counts"] != [4]:
            raise RuntimeError(f"{name} is not a pure TET mesh: {inventory}")

    reflection_is_exact = all(
        row["missing_reflected_vertices"] == 0
        and row["maximum_reflected_vertex_error_m"] <= 1e-13
        and row["missing_reflected_volume_elements"] == 0
        for row in reflection.values()
    )
    kelvin_is_periodic = bool(
        kelvin_identification["pair_count"] > 0
        and kelvin_identification["maximum_pair_translation_error_m"] is not None
        and kelvin_identification["maximum_pair_translation_error_m"] <= 1e-12
        and kelvin_fes["slaved_free_dofs"] > 0
        and kelvin_fes["trace_norm_ratio"] is not None
        and abs(kelvin_fes["trace_norm_ratio"] - 1.0) <= 1e-12
    )

    iron_rows = checks["iron"]["payload"].get("materials", [])
    iron_volume = next(
        (float(row["ng_volume"]) for row in iron_rows if row.get("name") == "iron"),
        None,
    )
    volume_error = None if iron_volume is None else (
        (iron_volume - EXACT_IRON_VOLUME_M3) / EXACT_IRON_VOLUME_M3
    )
    # The historical gate accepted TWO elements across the gap.  It stays the
    # default so an existing family does not silently change verdict, but it
    # is now a named, recorded requirement that a caller can raise: two
    # elements through the gap is not a field-accuracy claim.
    required_across_gap = float(options.gap_elements_across)
    # The requirement is on the MEASURED crossing count, not on the indicator:
    # sizing a volume does not by itself guarantee a division count, so the
    # mesh that was actually produced has to be walked.
    # The longest traversal a line makes inside one element is the number a
    # coarse stretch hides behind a healthy count, so it is required too.  It
    # is tied to the division requirement rather than fixed: at N divisions
    # the nominal traversal is gap_height/N, and 1.2 times that is the mesh
    # DESIGN tolerance -- not an accuracy guarantee.  Raising N therefore
    # tightens both requirements together.
    profile = gap_inventory["line_profile"]
    nominal = gap_inventory["gap_height_m"] / max(required_across_gap, 1)
    segment_limit = float(options.gap_segment_factor) * nominal
    gap_is_resolved = bool(
        gap_inventory["elements"] > 0
        and profile["minimum_segments"] >= required_across_gap
        and profile["curved_minimum_segments"] >= required_across_gap
        and profile["clean_traversal"]
        and profile["maximum_segment_m"] <= segment_limit
        and profile["curved_maximum_segment_m"] <= segment_limit
    )
    result = {
        "schema": "radia.validation.c-type-cubit-meshes.v1",
        "passed": bool(
            volume_error is not None
            and abs(volume_error) <= 1e-8
            and gap_is_resolved
            and reflection_is_exact
            and kelvin_is_periodic
        ),
        "machine": platform.node(),
        "cad_authority": str(CAD_JOURNAL),
        "cad_sha256": sha256(CAD_JOURNAL),
        "cubit": str(cubit),
        "command_plugin_dir": (
            str(command_plugin_dir) if command_plugin_dir is not None else None
        ),
        "mesh_policy": {
            "hdiv_mmm": "exact iron-only Cubit/ACIS TET .vol; Coulomb open boundary",
            "reduced_a": "exact iron + locally refined physical air + Kelvin exterior",
            "mixed_total_reduced_omega": (
                "same periodic Kelvin .vol as reduced-A; H1 TOSCA "
                "mixed total/reduced Omega"
            ),
        },
        "kelvin_radius_m": float(options.kelvin_radius),
        "kelvin_physical_center_m": [0.0, 0.0, 0.0],
        "iron_size_m": float(options.iron_size),
        "air_size_m": float(options.air_size),
        "gap_size_m": float(options.gap_size),
        "kelvin_mesh_size_m": float(options.kelvin_mesh_size),
        "curve_order": int(options.curve_order),
        "gap_inventory": gap_inventory,
        "gap_elements_across_required": int(options.gap_elements_across),
        "gap_segment_limit_m": segment_limit,
        "gap_segment_factor": float(options.gap_segment_factor),
        "reflection_inventory": reflection,
        "kelvin_identification": kelvin_identification,
        "kelvin_fes": kelvin_fes,
        "exact_iron_volume_m3": EXACT_IRON_VOLUME_M3,
        "mesh_iron_volume_m3": iron_volume,
        "relative_iron_volume_error": volume_error,
        "artifacts": {
            "iron_vol": str(iron_vol),
            "iron_vol_sha256": sha256(iron_vol),
            "kelvin_domain_vol": str(kelvin_vol),
            "kelvin_domain_vol_sha256": sha256(kelvin_vol),
        },
        "inventory": inventories,
        "checks": checks,
        "runs": runs,
    }
    result_path = output_dir / "mesh_result.json"
    result_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if not result["passed"]:
        raise RuntimeError(f"Cubit C-type mesh contract failed; see {result_path}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cubit", type=Path, default=CUBIT_DEFAULT)
    parser.add_argument(
        "--command-plugin-dir",
        type=Path,
        help="Cubit command-plugin directory containing cubit_mesh_export.ccm",
    )
    parser.add_argument("--iron-size", type=float, default=0.010)
    parser.add_argument("--air-size", type=float, default=0.015)
    parser.add_argument("--gap-size", type=float, default=0.002,
                        help="pole-face SURFACE size; this refines the gap in "
                             "plane and does not by itself divide the gap "
                             "through its thickness")
    parser.add_argument("--gap-elements-across", type=int, default=2,
                        help="required elements a line must cross through the "
                             "FULL gap thickness, measured on the produced "
                             "mesh. The historical default is 2; a "
                             "field-accuracy family starts at 6 and refines "
                             "from there")
    parser.add_argument("--gap-segment-factor", type=float, default=1.2,
                        help="longest traversal inside one element, as a "
                             "multiple of the nominal gap_height/N. A mesh "
                             "design tolerance, not an accuracy guarantee")
    parser.add_argument("--kelvin-radius", type=float, default=0.22)
    parser.add_argument("--kelvin-mesh-size", type=float, default=0.025)
    parser.add_argument("--curve-order", type=int, choices=range(2, 6), default=2)
    options = parser.parse_args()
    if any(value <= 0.0 for value in (
        options.iron_size,
        options.air_size,
        options.gap_size,
        options.kelvin_radius,
        options.kelvin_mesh_size,
    )):
        raise ValueError("mesh sizes must be positive")
    build(options)


if __name__ == "__main__":
    main()
