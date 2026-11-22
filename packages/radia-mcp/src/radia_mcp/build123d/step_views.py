"""Local orthographic drawings from imported STEP geometry."""

from math import floor, hypot
from pathlib import Path
from tempfile import mkdtemp


_CAMERAS = {
    "isometric": ((1, -1, 1), (0, 0, 1)),
    "xy": ((0, 0, 1), (0, 1, 0)),
    "xz": ((0, -1, 0), (0, 0, 1)),
    "yz": ((1, 0, 0), (0, 0, 1)),
}
_SECTION_AXES = {"xy": (0, 1, 2), "xz": (0, 2, 1), "yz": (1, 2, 0)}
_WIDTH, _HEIGHT, _MARGIN = 1024, 768, 24
_VISIBLE, _HIDDEN = "#20262e", "#808080"
_DASH, _GAP = 6, 4


def _edge_paths(edges, axes, deflection):
    from OCP.BRep import BRep_Tool
    from OCP.BRepAdaptor import BRepAdaptor_Curve
    from OCP.GCPnts import GCPnts_QuasiUniformDeflection

    paths = []
    for edge in edges:
        if BRep_Tool.Degenerated_s(edge.wrapped):
            continue
        samples = GCPnts_QuasiUniformDeflection(BRepAdaptor_Curve(edge.wrapped), deflection)
        if not samples.IsDone() or samples.NbPoints() < 2:
            raise ValueError("Could not sample an edge for drawing")
        points = []
        for index in range(1, samples.NbPoints() + 1):
            point = samples.Value(index).Coord()
            points.append((point[axes[0]], point[axes[1]]))
        paths.append(points)
    return paths


def _section_edges(shape, axis, offset):
    from build123d import Shape
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.gp import gp_Dir, gp_Pln, gp_Pnt

    origin = [0.0] * 3
    normal = [0.0] * 3
    origin[axis], normal[axis] = offset, 1.0
    operation = BRepAlgoAPI_Section(
        shape.wrapped, gp_Pln(gp_Pnt(*origin), gp_Dir(*normal))
    )
    if not operation.IsDone():
        raise ValueError("STEP plane intersection failed")
    cut = operation.Shape()
    return [] if cut.IsNull() else Shape.cast(cut).edges()


def _dashed_segments(path):
    phase = 0.0
    for (x0, y0), (x1, y1) in zip(path, path[1:]):
        dx, dy = x1 - x0, y1 - y0
        length = hypot(dx, dy)
        if not length:
            continue
        start, end = 0.0, 1.0
        for origin, delta, limit in ((x0, dx, _WIDTH), (y0, dy, _HEIGHT)):
            if delta:
                lo, hi = sorted(((-2 - origin) / delta, (limit + 2 - origin) / delta))
                start, end = max(start, lo), min(end, hi)
            elif not -2 <= origin <= limit + 2:
                end = -1.0
        if start < end:
            period = _DASH + _GAP
            first, last = start * length, end * length
            for index in range(floor((first + phase) / period), floor((last + phase) / period) + 1):
                dash = index * period - phase
                a, b = max(dash, first), min(dash + _DASH, last)
                if a < b:
                    yield [(x0 + dx * s / length, y0 + dy * s / length) for s in (a, b)]
        phase = (phase + length) % (_DASH + _GAP)


def _write_drawing(visible, hidden, base: Path, formats):
    artifacts = {}
    if "svg" in formats:
        lines = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{_WIDTH}" '
            f'height="{_HEIGHT}" viewBox="0 0 {_WIDTH} {_HEIGHT}">',
            '<rect width="100%" height="100%" fill="white"/>',
        ]
        for layer, paths, color, style in (
            ("hidden", hidden, _HIDDEN, f'stroke-dasharray="{_DASH} {_GAP}"'),
            ("visible", visible, _VISIBLE, 'stroke-linecap="round" stroke-linejoin="round"'),
        ):
            if paths:
                lines.append(f'<g id="{layer}" fill="none" stroke="{color}" stroke-width="1.5" {style}>')
                for path in paths:
                    points = " ".join(f"{x:.3f},{y:.3f}" for x, y in path)
                    lines.append(f'<polyline points="{points}"/>')
                lines.append("</g>")
        lines.append("</svg>")
        output = base.with_suffix(".svg")
        output.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
        artifacts["svg"] = str(output)
    if "png" in formats:
        from PIL import Image, ImageDraw

        oversample = 2
        image = Image.new("RGB", (_WIDTH * oversample, _HEIGHT * oversample), "white")
        draw = ImageDraw.Draw(image)
        for path in hidden:
            for segment in _dashed_segments(path):
                draw.line([(x * oversample, y * oversample) for x, y in segment],
                          fill=_HIDDEN, width=3)
        for path in visible:
            draw.line([(x * oversample, y * oversample) for x, y in path],
                      fill=_VISIBLE, width=3, joint="curve")
        output = base.with_suffix(".png")
        image.resize((_WIDTH, _HEIGHT), Image.Resampling.LANCZOS).save(output)
        artifacts["png"] = str(output)
    return artifacts


def render_views(shape, target: Path, formats, section_plane, section_offset_mm,
                 show_hidden, frame):
    from build123d import Vector

    bbox = shape.bounding_box()
    if not shape.edges() or bbox.diagonal <= 0:
        raise ValueError("STEP contains no drawable geometry")
    center = bbox.center() if frame is None else Vector(frame["center_mm"])
    deflection = bbox.diagonal / 4096
    views = {}
    for name, (direction, up) in _CAMERAS.items():
        eye = center + Vector(direction).normalized() * bbox.diagonal
        visible, hidden = shape.project_to_viewport(eye, up, look_at=center)
        outlines = _edge_paths(visible, (0, 1), deflection)
        if not outlines:
            raise ValueError(f"STEP projection produced no visible edges for {name}")
        views[name] = (outlines, _edge_paths(hidden, (0, 1), deflection) if show_hidden else [])

    x, y, axis = _SECTION_AXES[section_plane]
    coordinates = tuple(center)
    offset = coordinates[axis] if section_offset_mm is None else section_offset_mm
    section = _edge_paths(_section_edges(shape, axis, offset), (x, y), deflection)
    if not section:
        raise ValueError(f"Section {section_plane} at {offset:g} mm produces no contour")
    views["section"] = (
        [[(px - coordinates[x], py - coordinates[y]) for px, py in path] for path in section], []
    )

    def extent(layers, half_width, half_height):
        return max(max(abs(px) / half_width, abs(py) / half_height)
                   for paths in layers for path in paths for px, py in path)

    scale = frame["pixels_per_mm"] if frame is not None else 1 / max(
        extent(layers, _WIDTH / 2 - _MARGIN, _HEIGHT / 2 - _MARGIN) for layers in views.values()
    )
    extents = {name: extent(layers, _WIDTH / 2, _HEIGHT / 2) * scale
               for name, layers in views.items()}
    if max(extents.values()) > 1_000_000:
        raise ValueError("Fixed frame places geometry beyond the drawable coordinate range")
    target.mkdir(parents=True, exist_ok=True)
    target = Path(mkdtemp(prefix="build123d-views-", dir=target))
    output = {}
    clipped = []
    for name, layers in views.items():
        if extents[name] > 1:
            clipped.append(name)
        visible, hidden = (
            [[(round(_WIDTH / 2 + px * scale, 3), round(_HEIGHT / 2 - py * scale, 3))
              for px, py in path] for path in paths]
            for paths in layers
        )
        output[name] = _write_drawing(visible, hidden, target / name, formats)
    return {
        "views": output,
        "pixels_per_mm": scale,
        "frame": {"center_mm": list(center), "pixels_per_mm": scale},
        "clipped_views": clipped,
        "section": {"plane": section_plane, "offset_mm": offset},
    }
