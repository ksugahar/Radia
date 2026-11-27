"""Guards: parameter regeneration and a shared STEP readback for either CAD source."""

import asyncio
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

pytest.importorskip("build123d")

from radia_mcp.build123d import server  # noqa: E402


# Public two-block model. The boxes have x widths 10 and 6 mm.
SCRIPT = """
def build(params):
    left = Box(10, 4, 3)
    right = Pos(8 + params["gap_mm"], 0, 0) * Box(6, 4, 3)
    left.label = "left_block"
    right.label = "right_block"
    return Compound(label="gap_sample", children=[left, right])
"""


@pytest.mark.parametrize("value", ["x", "2.0", True, float("inf"), 10**400])
def test_invalid_parameter_returns_input_error(value, tmp_path):
    export_dir = tmp_path / "runs"
    result = json.loads(asyncio.run(server.build123d_run(
        SCRIPT, {"gap_mm": value}, str(export_dir)
    )))
    assert (result["status"], result["stage"], result["kind"]) == ("error", "input", "input")
    assert "gap_mm" in result["error"]
    assert not export_dir.exists()


def test_gap_regeneration_and_standalone_step_readback(tmp_path):
    results = [
        json.loads(asyncio.run(server.build123d_run(
            SCRIPT, {"gap_mm": gap}, str(tmp_path)
        )))
        for gap in (2.0, 0.5, -0.5)
    ]
    paths = []
    for result, gap in zip(results, (2.0, 0.5, -0.5), strict=True):
        assert result["status"] == "ok", result
        assert result["params"] == {"gap_mm": gap}
        assert result["script_sha256"] == hashlib.sha256(SCRIPT.encode()).hexdigest()
        step = Path(result["exported"])
        paths.append(step)
        measured = result["readback"]
        assert measured["sha256"] == hashlib.sha256(step.read_bytes()).hexdigest()
        assert measured["solid_count"] == 2
        assert measured["labels_unique_and_present"]
        volumes = {solid["label"]: solid["volume_mm3"] for solid in measured["solids"]}
        assert volumes == pytest.approx({"left_block": 120.0, "right_block": 72.0})
        assert measured["clearance_mm"] == pytest.approx(max(gap, 0))
        assert measured["interference_volume_mm3"] == pytest.approx(max(-gap, 0) * 12)

    assert len(set(paths)) == 3
    imported_step = tmp_path / "external.stp"
    imported_step.write_bytes(paths[0].read_bytes())
    assert json.loads(server.build123d_probe(
        str(imported_step), "step_readback"
    )) == results[0]["readback"]

    diff = json.loads(asyncio.run(server.build123d_diff_shapes(
        str(imported_step), str(paths[1])
    )))
    assert diff["status"] == "ok", diff
    assert diff["inputs"]["before"]["sha256"] == results[0]["readback"]["sha256"]
    assert diff["inputs"]["after"]["sha256"] == results[1]["readback"]["sha256"]
    assert diff["unresolved"] == []
    parts = {part["before"]["label"]: part["comparison"] for part in diff["matched_parts"]}
    assert set(parts) == {"left_block", "right_block"}
    assert parts["left_block"]["symmetric_difference_volume_mm3"] == pytest.approx(0)
    for key in ("bbox_min_mm", "bbox_max_mm"):
        assert parts["left_block"]["delta"][key] == pytest.approx([0, 0, 0])
        assert parts["right_block"]["delta"][key] == pytest.approx([-1.5, 0, 0])
    for part in parts.values():
        assert part["delta"]["volume_mm3"] == pytest.approx(0)
        assert part["delta"]["topology"] == dict(solids=0, faces=0, edges=0, vertices=0)
    assert parts["right_block"]["symmetric_difference_volume_mm3"] == pytest.approx(36)
    assert diff["assembly"]["symmetric_difference_volume_mm3"] == pytest.approx(36)


def test_regenerated_step_preserves_regions_in_checked_netgen_vol(tmp_path):
    """Guards: STEP-to-Netgen preserves named volumes after a gap edit."""
    if any(importlib.util.find_spec(name) is None for name in ("netgen", "ngsolve")):
        pytest.skip("Netgen and NGSolve are required by the STEP-to-VOL pipeline")
    pytest.importorskip("radia.vol_check")
    from radia_mcp.build123d.pipeline import _stage_mesh_multi

    expected_volumes = {"left_block": 120.0, "right_block": 72.0}
    for gap_mm in (2.0, 0.5):
        run = json.loads(asyncio.run(server.build123d_run(
            SCRIPT, {"gap_mm": gap_mm}, str(tmp_path),
        )))
        assert run["status"] == "ok", run
        step = Path(run["exported"])
        readback = json.loads(server.build123d_probe(str(step), "step_readback"))
        assert readback["status"] == "ok", readback
        assert readback["clearance_mm"] == pytest.approx(gap_mm)
        assert readback["interference_volume_mm3"] == 0
        assert readback["labels_unique_and_present"], readback
        labels = [solid["label"] for solid in readback["solids"]]
        assert {solid["label"]: solid["volume_mm3"] for solid in readback["solids"]} == (
            pytest.approx(expected_volumes)
        )

    result = _stage_mesh_multi(step, labels, tmp_path, "edited_gap", maxh=1.0, units="mm")
    contract = json.loads(Path(result["label_contract"]).read_text(encoding="utf-8"))
    report = json.loads(Path(result["vol_check"]).read_text(encoding="utf-8"))
    assert report["passed"], report
    assert contract["strict_labels"] and contract["units"] == "m", contract
    assert sorted(contract["required"]["materials"]) == sorted(expected_volumes)

    assert report["labels"]["passed"], report["labels"]
    assert sorted(report["labels"]["materials"]) == sorted(expected_volumes)
    volumes = {row["name"]: row["ng_volume"] * 1e9 for row in report["materials"]}
    assert volumes == pytest.approx(expected_volumes)


def test_assembly_readback_reports_every_pair_and_label_issue(tmp_path):
    from build123d import Box, Compound, Pos, export_step

    a = Box(2, 2, 2)
    b = Pos(2, 0, 0) * Box(2, 2, 2)
    c = Pos(1, -1, 0) * Box(2, 2, 2)
    d = Pos(2, 2, 0) * Box(2, 2, 2)
    a.label, b.label, c.label = "anchor", "bracket", "Bracket"
    step = tmp_path / "assembly.stp"
    assert export_step(Compound(label="assembly_title", children=[a, b, c, d]), step)

    result = json.loads(server.build123d_probe(str(step), "step_readback"))
    assert result["status"] == "ok", result
    assert result["solid_count"] == 4
    ids = {solid["label"]: solid["id"] for solid in result["solids"]}
    assert set(ids) == {"anchor", "bracket", "Bracket", ""}
    assert result["label_audit"]["missing_solid_ids"] == [ids[""]]
    [duplicates] = result["label_audit"]["duplicate_solid_id_groups"]
    assert set(duplicates) == {ids["bracket"], ids["Bracket"]}
    assert not result["labels_unique_and_present"]

    pairs = {frozenset(row["solid_ids"]): row for row in result["pairs"]}
    assert len(result["pairs"]) == len(pairs) == 6
    face = pairs[frozenset((ids["anchor"], ids["bracket"]))]
    assert face["clearance_mm"] == pytest.approx(0)
    assert face["interference_volume_mm3"] == pytest.approx(0)
    assert face["shared_face_area_mm2"] == pytest.approx(4)
    assert face["has_face_contact"]

    edge = pairs[frozenset((ids["anchor"], ids[""]))]
    assert edge["clearance_mm"] == pytest.approx(0)
    assert edge["shared_face_area_mm2"] == pytest.approx(0)
    assert not edge["has_face_contact"]

    overlap = pairs[frozenset((ids["anchor"], ids["Bracket"]))]
    assert overlap["interference_volume_mm3"] == pytest.approx(2)
    assert overlap["clearance_mm"] == pytest.approx(0)
    assert not overlap["has_face_contact"]

    clear = pairs[frozenset((ids["Bracket"], ids[""]))]
    assert clear["clearance_mm"] == pytest.approx(1)
    assert clear["interference_volume_mm3"] == pytest.approx(0)
    assert not clear["has_face_contact"]

    c.label = "renamed_bracket"
    relabeled = tmp_path / "relabeled.step"
    assert export_step(Compound(label="assembly_title", children=[a, b, c, d]), relabeled)
    diff = json.loads(server._diff_step_files(str(step), str(relabeled)))
    assert diff["status"] == "ok", diff
    assert [part["before"]["label"] for part in diff["matched_parts"]] == ["anchor"]
    issues = {entry["label_key"]: entry for entry in diff["unresolved"]}
    assert issues[""]["issues"] == ["missing_label"]
    assert issues["bracket"]["issues"] == ["duplicate_label"]
    assert len(issues["bracket"]["before"]) == 2
    assert issues["renamed_bracket"]["issues"] == ["missing_match"]
    assert diff["assembly"]["before"]["volume_mm3"] == pytest.approx(28)
    assert diff["assembly"]["symmetric_difference_volume_mm3"] == pytest.approx(0)


def _pins(tmp_path, name, offsets):
    from build123d import Box, Compound, Pos, export_step

    base = Box(30, 6, 2)
    base.label = "base"
    pins = []
    for x, height in offsets:
        pin = Pos(x, 0, 1 + height / 2) * Box(2, 2, height)
        pin.label = "pin"
        pins.append(pin)
    step = tmp_path / f"{name}.step"
    assert export_step(Compound(label=name, children=[base, *pins]), step)
    return step


def test_repeated_labels_pair_by_unambiguous_placement(tmp_path):
    """Guards: same-label instances pair by centroid; edits stay attached to the right one."""
    before = _pins(tmp_path, "before", [(-10, 4), (10, 4)])
    after = _pins(tmp_path, "after", [(-10, 4), (10.5, 6)])
    diff = json.loads(server._diff_step_files(str(before), str(after)))
    assert diff["status"] == "ok", diff
    assert diff["unresolved"] == []
    pins = [part for part in diff["matched_parts"] if part["before"]["label"] == "pin"]
    assert [part["match"] for part in pins] == ["placement", "placement"]
    deltas = sorted(part["comparison"]["delta"]["volume_mm3"] for part in pins)
    assert deltas == pytest.approx([0, 8])
    moved = max(pins, key=lambda part: part["centroid_distance_mm"])
    assert moved["centroid_distance_mm"] == pytest.approx((0.5**2 + 1**2) ** 0.5)
    assert [part["match"] for part in diff["matched_parts"]
            if part["before"]["label"] == "base"] == ["label"]


def test_ambiguous_repeated_labels_stay_unresolved(tmp_path):
    """Guards: a pattern shifted by half its pitch is not paired by guesswork."""
    before = _pins(tmp_path, "before", [(-10, 4), (0, 4)])
    after = _pins(tmp_path, "after", [(-5, 4), (5, 4)])
    diff = json.loads(server._diff_step_files(str(before), str(after)))
    assert diff["status"] == "ok", diff
    [entry] = [row for row in diff["unresolved"] if row["label_key"] == "pin"]
    assert entry["issues"] == ["duplicate_label", "ambiguous_placement"]
    assert len(entry["centroid_distance_mm"]) == 4
    assert all(part["before"]["label"] != "pin" for part in diff["matched_parts"])
