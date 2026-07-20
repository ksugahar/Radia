"""Thread, sheet-metal, STEP-PMI, and repaired-BREP identity checks."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


THREAD = "screwthread_pitch_handedness_start_topology_volume_owner_identity"
SHEET = (
    "sheetmetal_thickness_bendallowance_neutralaxis_flatpattern_owner_identity"
)
PMI = "step_pmi_tolerance_datum_unit_product_owner_identity"
BREP = "brep_repair_tolerance_sewnshell_orientation_volume_owner_identity"
_LENGTH_UNITS = {"m", "mm", "cm", "inch"}


def _digest(value: object) -> bool:
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def _generation(row: Mapping[str, object], *names: str) -> bool:
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def _number(value: object, *, positive: bool = False, nonnegative: bool = False) -> bool:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        return False
    number = float(value)
    return (
        math.isfinite(number)
        and (not positive or number > 0.0)
        and (not nonnegative or number >= 0.0)
    )


def _close(left: object, right: object) -> bool:
    return _number(left) and _number(right) and math.isclose(
        float(left), float(right), rel_tol=1.0e-10, abs_tol=1.0e-12
    )


def _result(row: Mapping[str, object]) -> bool:
    return _digest(row.get("result_sha256")) and row.get(
        "accepted_result_sha256"
    ) == row.get("result_sha256")


def _topology(value: object) -> bool:
    if not isinstance(value, Mapping) or set(value) != {
        "solids",
        "shells",
        "faces",
        "edges",
        "vertices",
    }:
        return False
    if not all(
        isinstance(item, int) and not isinstance(item, bool) and item >= 0
        for item in value.values()
    ):
        return False
    return (
        value["solids"] == 1
        and value["shells"] == 1
        and value["faces"] >= 4
        and value["edges"] >= 4
        and value["vertices"] >= 2
        and value["vertices"] - value["edges"] + value["faces"] == 2
    )


def _thread_ok(row: Mapping[str, object]) -> bool:
    names = (
        "pitch_m",
        "handedness",
        "start_count",
        "thread_topology",
        "volume_m3",
        "shape_owner",
    )
    return (
        _generation(
            row,
            "pitch_generation",
            "handedness_generation",
            "start_generation",
            "topology_generation",
            "volume_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(row.get("pitch_m"), positive=True)
        and row.get("handedness") in {"left", "right"}
        and isinstance(row.get("start_count"), int)
        and not isinstance(row.get("start_count"), bool)
        and 1 <= int(row["start_count"]) <= 8
        and _topology(row.get("thread_topology"))
        and _number(row.get("volume_m3"), positive=True)
        and str(row.get("shape_owner") or "").startswith("shape:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _sheet_ok(row: Mapping[str, object]) -> bool:
    thickness = row.get("thickness_m")
    radius = row.get("inside_bend_radius_m")
    angle = row.get("bend_angle_rad")
    neutral_axis = row.get("neutral_axis_factor")
    bend_allowance = row.get("bend_allowance_m")
    expected_allowance = (
        abs(float(angle)) * (float(radius) + float(neutral_axis) * float(thickness))
        if all(_number(value) for value in (thickness, radius, angle, neutral_axis))
        else math.nan
    )
    names = (
        "thickness_m",
        "inside_bend_radius_m",
        "bend_angle_rad",
        "neutral_axis_factor",
        "bend_allowance_m",
        "flat_pattern_area_m2",
        "folded_surface_area_m2",
        "shape_owner",
    )
    return (
        _generation(
            row,
            "thickness_generation",
            "bend_generation",
            "neutralaxis_generation",
            "pattern_generation",
            "area_generation",
            "owner_generation",
            "result_generation",
        )
        and _number(thickness, positive=True)
        and _number(radius, nonnegative=True)
        and _number(angle)
        and 0.0 < abs(float(angle)) <= 2.0 * math.pi
        and _number(neutral_axis, nonnegative=True)
        and 0.0 <= float(neutral_axis) <= 1.0
        and _number(bend_allowance, positive=True)
        and _close(bend_allowance, expected_allowance)
        and _number(row.get("flat_pattern_area_m2"), positive=True)
        and _number(row.get("folded_surface_area_m2"), positive=True)
        and _close(row.get("flat_pattern_area_m2"), row.get("folded_surface_area_m2"))
        and str(row.get("shape_owner") or "").startswith("shape:")
        and all(row.get("result_" + name) == row.get(name) for name in names)
        and _result(row)
    )


def _pmi_ok(row: Mapping[str, object]) -> bool:
    tolerances = row.get("geometric_tolerances")
    datums = row.get("datum_frame")
    associations = row.get("product_association")
    if (
        not isinstance(tolerances, Sequence)
        or isinstance(tolerances, (str, bytes))
        or not tolerances
        or not isinstance(datums, Mapping)
        or not datums
        or not isinstance(associations, Mapping)
        or not associations
    ):
        return False
    if not all(
        isinstance(name, str)
        and name
        and isinstance(entity, str)
        and entity.startswith(("face:", "axis:", "plane:"))
        for name, entity in datums.items()
    ):
        return False
    features: list[str] = []
    for tolerance in tolerances:
        if not isinstance(tolerance, Mapping) or set(tolerance) != {
            "feature",
            "kind",
            "value",
            "datum_refs",
        }:
            return False
        feature = tolerance["feature"]
        refs = tolerance["datum_refs"]
        if (
            not isinstance(feature, str)
            or not feature.startswith("feature:")
            or not isinstance(tolerance["kind"], str)
            or not tolerance["kind"]
            or not _number(tolerance["value"], positive=True)
            or not isinstance(refs, Sequence)
            or isinstance(refs, (str, bytes))
            or not refs
            or len(refs) != len(set(refs))
            or not set(refs).issubset(datums)
        ):
            return False
        features.append(feature)
    associations_ok = set(associations) == set(features) and all(
        isinstance(product, str) and product.startswith("part:")
        for product in associations.values()
    )
    return (
        _generation(
            row,
            "tolerance_generation",
            "datum_generation",
            "unit_generation",
            "product_generation",
            "revision_generation",
            "owner_generation",
            "result_generation",
        )
        and len(features) == len(set(features))
        and associations_ok
        and row.get("length_unit") in _LENGTH_UNITS
        and all(
            row.get("replayed_" + name) == row.get(name)
            for name in (
                "geometric_tolerances",
                "datum_frame",
                "length_unit",
                "product_association",
                "document_revision",
                "document_owner",
            )
        )
        and str(row.get("document_revision") or "").startswith("document:")
        and str(row.get("document_owner") or "").startswith("document:")
        and _result(row)
    )


def _brep_ok(row: Mapping[str, object]) -> bool:
    shells = row.get("sewn_shells")
    orientations = row.get("face_orientations")
    if (
        not isinstance(shells, Sequence)
        or isinstance(shells, (str, bytes))
        or not shells
        or not isinstance(orientations, Mapping)
        or not orientations
    ):
        return False
    faces: list[str] = []
    shell_names: list[str] = []
    for shell in shells:
        if not isinstance(shell, Mapping) or set(shell) != {"shell", "faces", "closed"}:
            return False
        shell_name = shell["shell"]
        shell_faces = shell["faces"]
        if (
            not isinstance(shell_name, str)
            or not shell_name.startswith("shell:")
            or not isinstance(shell_faces, Sequence)
            or isinstance(shell_faces, (str, bytes))
            or len(shell_faces) < 4
            or len(shell_faces) != len(set(shell_faces))
            or not all(isinstance(face, str) and face.startswith("face:") for face in shell_faces)
            or shell["closed"] is not True
        ):
            return False
        shell_names.append(shell_name)
        faces.extend(shell_faces)
    orientation_ok = set(orientations) == set(faces) and all(
        isinstance(value, int) and not isinstance(value, bool) and value in {-1, 1}
        for value in orientations.values()
    )
    return (
        _generation(
            row,
            "tolerance_generation",
            "shell_generation",
            "orientation_generation",
            "volume_generation",
            "owner_generation",
            "result_generation",
        )
        and len(shell_names) == len(set(shell_names))
        and len(faces) == len(set(faces))
        and orientation_ok
        and _number(row.get("healing_tolerance_m"), positive=True)
        and isinstance(row.get("closed_volume_count"), int)
        and not isinstance(row.get("closed_volume_count"), bool)
        and row.get("closed_volume_count") == len(shells)
        and _number(row.get("volume_m3"), positive=True)
        and str(row.get("shape_owner") or "").startswith("shape:")
        and all(
            row.get("replayed_" + name) == row.get(name)
            for name in (
                "healing_tolerance_m",
                "sewn_shells",
                "face_orientations",
                "closed_volume_count",
                "volume_m3",
                "shape_owner",
            )
        )
        and _result(row)
    )


def _public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    rows: list[Mapping[str, object]] = []
    reference = payload.get("reference")
    if isinstance(reference, Sequence) and not isinstance(reference, (str, bytes)):
        rows.extend(item for item in reference if isinstance(item, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for family in measured.values():
            if isinstance(family, Sequence) and not isinstance(family, (str, bytes)):
                rows.extend(item for item in family if isinstance(item, Mapping))
    return rows


def _report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }


def validate_public_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping):
        return {}
    rows = _public_rows(payload)
    checks: dict[str, bool] = {}
    threads = [row.get(THREAD) for row in rows if THREAD in row]
    sheets = [row.get(SHEET) for row in rows if SHEET in row]
    if threads:
        checks["v55_thread_pitch_hand_start_topology_volume_owner"] = (
            len(threads) == len(rows)
            and all(isinstance(item, Mapping) and _thread_ok(item) for item in threads)
        )
    if sheets:
        checks["v55_sheet_thickness_bend_neutralaxis_pattern_owner"] = (
            len(sheets) == len(rows)
            and all(isinstance(item, Mapping) and _sheet_ok(item) for item in sheets)
        )
    return _report("build123d_v55_public_identity_v1", checks) if checks else {}


def validate_source_identity(payload: object) -> dict[str, object]:
    if not isinstance(payload, Mapping) or not isinstance(payload.get("replay_identity"), Mapping):
        return {}
    identity = payload["replay_identity"]
    checks: dict[str, bool] = {}
    if identity.get(PMI) is not None:
        checks["v55_step_pmi_tolerance_datum_product_owner"] = (
            isinstance(identity[PMI], Mapping) and _pmi_ok(identity[PMI])
        )
    if identity.get(BREP) is not None:
        checks["v55_brep_repair_shell_orientation_volume_owner"] = (
            isinstance(identity[BREP], Mapping) and _brep_ok(identity[BREP])
        )
    return _report("build123d_v55_source_identity_v1", checks) if checks else {}
