"""Guards: headless STEP projections, section placement and SVG/PNG output."""

import asyncio
import hashlib
import json
from pathlib import Path
from xml.etree import ElementTree

import pytest

pytest.importorskip("build123d")

from radia_mcp.build123d import server  # noqa: E402


_SVG = "{http://www.w3.org/2000/svg}"


def _layer_bounds(svg, layer):
    root = ElementTree.parse(svg).getroot()
    lines = root.findall(f"./{_SVG}g[@id='{layer}']/{_SVG}polyline")
    xs, ys = zip(*(tuple(map(float, point.split(","))) for line in lines
                   for point in line.attrib["points"].split()))
    return min(xs), min(ys), max(xs), max(ys)


@pytest.fixture
def step(tmp_path):
    from build123d import Box, Cylinder, Pos, export_step

    part = Pos(10, 20, 30) * (Box(12, 8, 4) - Cylinder(1, 6))
    path = tmp_path / "input.stp"
    assert export_step(part, path)
    return path


def test_step_views_write_both_formats_with_shared_scale(step, tmp_path):
    image_module = pytest.importorskip("PIL.Image")
    result = json.loads(asyncio.run(server.build123d_render_views(
        str(step), str(tmp_path), section_offset_mm=31
    )))
    assert result["status"] == "ok", result
    assert result["step_sha256"] == hashlib.sha256(step.read_bytes()).hexdigest()
    assert result["section"] == {"plane": "xy", "offset_mm": 31.0}
    assert set(result["views"]) == {"isometric", "xy", "xz", "yz", "section"}
    expected = {"xy": (12, 8), "xz": (12, 4), "yz": (8, 4), "section": (12, 8)}
    for name, paths in result["views"].items():
        root = ElementTree.parse(paths["svg"]).getroot()
        lines = root.findall(f"./{_SVG}g[@id='visible']/{_SVG}polyline")
        assert lines
        if name in expected:
            points = [tuple(map(float, pair.split(","))) for line in lines
                      for pair in line.attrib["points"].split()]
            xs, ys = zip(*points)
            scale = result["pixels_per_mm"]
            assert ((max(xs) - min(xs)) / scale, (max(ys) - min(ys)) / scale
                    ) == pytest.approx(expected[name], abs=0.001)
        if name == "section":
            assert any(len(line.attrib["points"].split()) > 10 for line in lines)
        with image_module.open(paths["png"]) as png:
            assert png.format == "PNG" and png.size == (1024, 768)
            assert png.getextrema() != ((255, 255),) * 3
        assert Path(paths["svg"]).parent == Path(paths["png"]).parent

    side = result["views"]["xz"]
    hidden = ElementTree.parse(side["svg"]).find(f"./{_SVG}g[@id='hidden']")
    assert hidden is not None and len(hidden) > 0
    assert hidden.attrib["stroke"] == "#808080"
    assert hidden.attrib["stroke-dasharray"] == "6 4"
    disabled = json.loads(server._render_step_views(
        str(step), str(tmp_path), ["svg", "png"], "xy", 31,
        show_hidden=False, frame=result["frame"]
    ))
    assert disabled["status"] == "ok", disabled
    plain = disabled["views"]["xz"]
    assert ElementTree.parse(plain["svg"]).find(f"./{_SVG}g[@id='hidden']") is None
    with image_module.open(side["png"]) as dashed, image_module.open(plain["png"]) as solid:
        x0, y0, x1, y1 = _layer_bounds(side["svg"], "visible")
        x = round((x0 + x1) / 2 - result["pixels_per_mm"])
        samples = range(round(y0) + 8, round(y1) - 8)
        marks = [dashed.getpixel((x, y)) for y in samples]
        assert any(r == g == b and 90 < r < 215 for r, g, b in marks)
        assert (255, 255, 255) in marks
        assert all(solid.getpixel((x, y)) == (255, 255, 255) for y in samples)

    from build123d import Pos, export_step, import_step

    shifted = tmp_path / "shifted.step"
    assert export_step(Pos(0.5, 0.25, 0) * import_step(step), shifted)
    moved = json.loads(server._render_step_views(
        str(shifted), str(tmp_path), ["svg", "png"], "xy", 31, frame=result["frame"]
    ))
    assert moved["status"] == "ok", moved
    assert moved["frame"] == result["frame"] and moved["section"] == result["section"]
    scale = result["pixels_per_mm"]
    for name in ("xy", "section"):
        before = _layer_bounds(result["views"][name]["svg"], "visible")
        after = _layer_bounds(moved["views"][name]["svg"], "visible")
        assert [b - a for a, b in zip(before, after)] == pytest.approx(
            [0.5 * scale, -0.25 * scale] * 2, abs=0.002
        )
    with image_module.open(result["views"]["xy"]["png"]) as before_png:
        with image_module.open(moved["views"]["xy"]["png"]) as after_png:
            def ink_bounds(png):
                return png.convert("L").point(lambda value: 255 if value < 100 else 0).getbbox()

            assert [b - a for a, b in zip(ink_bounds(before_png), ink_bounds(after_png))
                    ] == pytest.approx([0.5 * scale, -0.25 * scale] * 2, abs=1.1)


def test_outside_section_reports_error_without_writing(step, tmp_path):
    output = tmp_path / "no-output"
    result = json.loads(server._render_step_views(
        str(step), str(output), ["svg"], "xz", 40
    ))
    assert (result["status"], result["stage"]) == ("error", "render_views")
    assert "produces no contour" in result["error"]
    assert not output.exists()


def test_boolean_section_offset_is_rejected(step, tmp_path):
    output = tmp_path / "boolean-offset"
    result = json.loads(asyncio.run(server.build123d_render_views(
        str(step), str(output), ["svg"], section_offset_mm=True
    )))
    assert (result["status"], result["stage"], result["kind"]) == ("error", "input", "input")
    assert not output.exists()
