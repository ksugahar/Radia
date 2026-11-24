"""Guards: STEP declaration and millimetre measurement across CAD workflows."""

import json
import re
from pathlib import Path

import pytest

pytest.importorskip("build123d")

from radia_mcp.build123d import server  # noqa: E402


@pytest.fixture
def published_steps(tmp_path):
    from build123d import Box, Compound, Pos, export_step
    from OCP.Interface import Interface_Static
    from OCP.STEPControl import STEPControl_Controller

    a = Box(25.4, 12.7, 6.35)
    b = Pos(38.1, 0, 0) * Box(25.4, 12.7, 6.35)
    a.label, b.label = "fixed_block", "spaced_block"
    model = Compound(children=[a, b])
    STEPControl_Controller.Init_s()
    key = "write.step.unit"
    previous = Interface_Static.CVal_s(key)
    files = {}
    try:
        for unit in ("MM", "INCH"):
            assert Interface_Static.SetCVal_s(key, unit)
            files[unit] = tmp_path / f"two_blocks_{unit}.step"
            assert export_step(model, files[unit])
    finally:
        Interface_Static.SetCVal_s(key, previous)
    return files["MM"], files["INCH"]


def test_mm_and_inch_steps_produce_the_same_geometry(published_steps, tmp_path):
    from PIL import Image, ImageChops

    mm, inch = published_steps
    for path, source_unit in ((mm, "millimetre"), (inch, "inch")):
        result = json.loads(server._inspect_step_readback(str(path)))
        assert result["status"] == "ok", result
        assert [unit.casefold() for unit in result["declared_length_units"]] == [source_unit]
        assert result["measurement_length_unit"] == result["length_unit"] == "mm"
        assert result["clearance_mm"] == pytest.approx(12.7)
        volumes = {part["label"]: part["volume_mm3"] for part in result["solids"]}
        assert volumes == pytest.approx(dict.fromkeys(
            ("fixed_block", "spaced_block"), 25.4 * 12.7 * 6.35
        ))

    image_root = str(tmp_path / "views")
    before = json.loads(server._render_step_views(str(mm), image_root, ["svg", "png"], "xy", 0))
    assert before["status"] == "ok", before
    after = json.loads(server._render_step_views(
        str(inch), image_root, ["svg", "png"], "xy", 0, frame=before["frame"]
    ))
    assert after["status"] == "ok", after
    assert before["declared_length_units"] == ["millimetre"]
    assert [name.casefold() for name in after["declared_length_units"]] == ["inch"]
    assert before["measurement_length_unit"] == after["measurement_length_unit"] == "mm"
    assert before["frame"] == after["frame"]
    for name, paths in before["views"].items():
        assert Path(paths["svg"]).read_bytes() == Path(after["views"][name]["svg"]).read_bytes()
        with Image.open(paths["png"]) as first, Image.open(after["views"][name]["png"]) as second:
            assert ImageChops.difference(first, second).getbbox() is None

    diff = json.loads(server._diff_step_files(str(mm), str(inch)))
    assert diff["status"] == "ok", diff
    assert diff["inputs"]["before"]["declared_length_units"] == ["millimetre"]
    assert [name.casefold() for name in diff["inputs"]["after"]["declared_length_units"]] == ["inch"]
    assert diff["measurement_length_unit"] == "mm"
    assert diff["unresolved"] == []
    assert len(diff["matched_parts"]) == 2
    for comparison in [diff["assembly"], *(p["comparison"] for p in diff["matched_parts"])]:
        assert comparison["symmetric_difference_volume_mm3"] == pytest.approx(0)
        assert comparison["delta"]["bbox_min_mm"] == pytest.approx([0, 0, 0])
        assert comparison["delta"]["bbox_max_mm"] == pytest.approx([0, 0, 0])
        assert comparison["delta"]["volume_mm3"] == pytest.approx(0)
        assert not any(comparison["delta"]["topology"].values())


def test_unknown_step_unit_fails_in_all_three_workflows(published_steps, tmp_path):
    mm, inch = published_steps
    data = inch.read_bytes()
    marker = b"CONVERSION_BASED_UNIT('INCH',"
    first, found, rest = data.partition(marker)
    assert found and marker in rest
    unknown = tmp_path / "unknown.step"
    unknown.write_bytes(first + found + rest.replace(
        marker, b"CONVERSION_BASED_UNIT('UNRECOGNIZED_UNIT',", 1
    ))
    output = tmp_path / "no-views"
    for response in (
        server._inspect_step_readback(str(unknown)),
        server._render_step_views(str(unknown), str(output), ["svg"], "xy", None),
        server._diff_step_files(str(mm), str(unknown)),
    ):
        result = json.loads(response)
        assert result["status"] == "error", result
        assert "Unknown STEP declared length unit" in result["error"]
    assert not output.exists()


def test_shape_context_without_declared_length_is_rejected(published_steps, tmp_path):
    _, inch = published_steps
    data, replacements = re.subn(
        rb"(GLOBAL_UNIT_ASSIGNED_CONTEXT\s*\(\()#\d+,", rb"\1", inch.read_bytes(), count=1
    )
    assert replacements == 1
    missing = tmp_path / "missing_unit.step"
    missing.write_bytes(data)
    result = json.loads(server._inspect_step_readback(str(missing)))
    assert result["status"] == "error", result
    assert "declared length unit" in result["error"]
