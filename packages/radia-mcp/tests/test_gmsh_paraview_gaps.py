"""The five capabilities that still separated the gmsh lane from ParaView.

Each is either a real substitute for something gmsh lacks (volume
rendering, LIC) or a composition gmsh has no single verb for (shared
camera/scale panels, cross-file colour range, compound selection).  The
substitutes are named for what they DO, and these tests hold them to the
property that makes them honest:

  * ``volume_render`` composites N cut planes -- so the number of drawn
    slices must equal what was asked, and the result must say it is not
    ray casting.
  * ``flow_texture`` is evenly spaced streamlines -- so a higher density
    must actually produce a finer spacing, and the result must say it is
    not LIC.
  * ``render_panels`` must REFUSE to put two different quantities on one
    colour bar (a shared scale over T and A/m^2 means nothing).
  * ``select`` must reject an unknown name instead of evaluating to
    something plausible, and must carry field VALUES into the extraction.
"""

from __future__ import annotations

import pytest

from radia_mcp.gmsh.post_process import field_range, flow_texture, select
from radia_mcp.gmsh.render import render_panels, volume_render

pytest.importorskip("gmsh", reason="gmsh not installed")


# --------------------------------------------------------------------
# fixtures: a small analytic field on a box, written as a v4.1 .msh
# --------------------------------------------------------------------

def _write_msh(path, nodes, elements, views):
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat"]
    lo = [min(p[i] for p in nodes.values()) for i in range(3)]
    hi = [max(p[i] for p in nodes.values()) for i in range(3)]
    lines += ["$Entities", "0 0 0 1",
              "1 %.9g %.9g %.9g %.9g %.9g %.9g 0 0"
              % (lo[0], lo[1], lo[2], hi[0], hi[1], hi[2]),
              "$EndEntities"]
    lines += ["$Nodes", f"1 {len(nodes)} 1 {len(nodes)}",
              f"3 1 0 {len(nodes)}"]
    lines += [str(t) for t in nodes]
    lines += ["%.15e %.15e %.15e" % tuple(p) for p in nodes.values()]
    lines += ["$EndNodes"]
    lines += ["$Elements", f"1 {len(elements)} 1 {len(elements)}",
              f"3 1 4 {len(elements)}"]
    lines += [f"{t} " + " ".join(str(n) for n in ns)
              for t, ns in elements.items()]
    lines += ["$EndElements"]
    for name, ncomp, rows in views:
        lines += ["$NodeData", "1", f'"{name}"', "1", "0",
                  f"3", "0", str(ncomp), str(len(rows))]
        lines += [f"{t} " + " ".join("%.9e" % v for v in vals)
                  for t, vals in rows.items()]
        lines += ["$EndNodeData"]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _box_field(tmp_path, name, scale=1.0, n=5):
    """Regular tet-filled box with |v| = scale * (x + 1) on the nodes."""
    nodes, idx = {}, {}
    tag = 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                nodes[tag] = [i / (n - 1), j / (n - 1), k / (n - 1)]
                idx[(i, j, k)] = tag
                tag += 1
    # 5-tet decomposition of every cell: the box must be FILLED, or a
    # cut plane lands in empty space and the streamline tracer finds
    # nothing to follow (measured: a 1-tet-per-cell box gives 0 lines).
    elements = {}
    et = 1
    for i in range(n - 1):
        for j in range(n - 1):
            for k in range(n - 1):
                c = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                     idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                     idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                for a, b, cc, d in ((0, 1, 3, 4), (1, 2, 3, 6),
                                    (1, 3, 4, 6), (1, 4, 5, 6),
                                    (3, 4, 6, 7)):
                    elements[et] = [c[a], c[b], c[cc], c[d]]
                    et += 1
    rows = {t: [scale * (p[0] + 1.0)] for t, p in nodes.items()}
    vec = {t: [0.0, scale * (p[0] + 1.0), 0.0] for t, p in nodes.items()}
    return _write_msh(tmp_path / f"{name}.msh", nodes, elements,
                      [("f", 1, rows), ("vec", 3, vec)])


# --------------------------------------------------------------------
# 4. cross-file colour range
# --------------------------------------------------------------------

def test_field_range_unions_across_files(tmp_path):
    a = _box_field(tmp_path, "a", scale=1.0)     # f in [1, 2]
    b = _box_field(tmp_path, "b", scale=3.0)     # f in [3, 6]
    one = field_range([a], view="f")
    assert one["range"] == pytest.approx([1.0, 2.0])
    both = field_range([a, b], view="f")
    assert both["range"] == pytest.approx([1.0, 6.0])
    assert both["n_files"] == 2
    assert set(both["per_file"]) == {str(a), str(b)}
    bad = field_range([a], view="nope")
    assert not bad["ok"] and "no matching view" in bad["error"]


def test_field_range_component_vs_magnitude(tmp_path):
    a = _box_field(tmp_path, "a", scale=2.0)     # vec = (0, 2(x+1), 0)
    mag = field_range([a], view="vec")
    comp0 = field_range([a], view="vec", component=0)
    assert mag["range"] == pytest.approx([2.0, 4.0])
    assert comp0["range"] == pytest.approx([0.0, 0.0])


def test_field_range_handles_element_node_samples_individually(tmp_path):
    path = _write_msh(
        tmp_path / "element_node.msh",
        {1: [0, 0, 0], 2: [1, 0, 0], 3: [0, 1, 0], 4: [0, 0, 1]},
        {1: [1, 2, 3, 4]}, [])
    path.write_text(
        path.read_text(encoding="utf-8")
        + "$ElementNodeData\n1\n\"local\"\n1\n0\n3\n0\n1\n1\n"
          "1 4 -2 -1 3 4\n$EndElementNodeData\n",
        encoding="utf-8")
    result = field_range([path], view="local")
    assert result["ok"], result
    assert result["range"] == pytest.approx([-2.0, 4.0])


# --------------------------------------------------------------------
# 3. shared camera / zoom / colour panels
# --------------------------------------------------------------------

def test_render_panels_shares_the_range_and_frame(tmp_path):
    a = _box_field(tmp_path, "a", scale=1.0)
    b = _box_field(tmp_path, "b", scale=3.0)
    res = render_panels([{"path": str(a), "label": "A"},
                         {"path": str(b), "label": "B"}],
                        tmp_path / "panels.png", view="f",
                        camera_preset="+y", width=260, height=240)
    assert res.get("ok"), res
    assert res["shared_range"] == pytest.approx([1.0, 6.0])
    assert res["frame"] and res["shared_camera"]
    assert len(res["panels"]) == 2


def test_render_panels_refuses_to_share_across_quantities(tmp_path):
    a = _box_field(tmp_path, "a")
    other = _write_msh(tmp_path / "other.msh",
                       {1: [0, 0, 0], 2: [1, 0, 0], 3: [0, 1, 0],
                        4: [0, 0, 1]},
                       {1: [1, 2, 3, 4]},
                       [("totally_different", 1, {1: [1.0], 2: [2.0],
                                                  3: [3.0], 4: [4.0]})])
    with pytest.raises(ValueError, match="no view in common"):
        render_panels([str(a), str(other)], tmp_path / "x.png",
                      width=200, height=200)
    # explicitly opting out is allowed
    res = render_panels([str(a), str(other)], tmp_path / "ok.png",
                        share_color=False, width=200, height=200)
    assert res.get("ok"), res


def test_render_panels_requires_one_common_quantity(tmp_path):
    a = _box_field(tmp_path, "a")
    b = _box_field(tmp_path, "b")
    with pytest.raises(ValueError, match="multiple views in common"):
        render_panels([str(a), str(b)], tmp_path / "ambiguous.png",
                      width=200, height=200)


def test_render_panels_accepts_shared_merges_and_item_overrides(
        tmp_path, monkeypatch):
    import radia_mcp.gmsh.render as render_module
    from PIL import Image

    a = _box_field(tmp_path, "a")
    shared = tmp_path / "shared.geo"
    local = tmp_path / "local.geo"
    shared.write_text("Point(1) = {0, 0, 0, 1};\n", encoding="utf-8")
    local.write_text("Point(2) = {1, 0, 0, 1};\n", encoding="utf-8")
    calls = []

    def fake_render(path, png, **kwargs):
        calls.append(kwargs)
        Image.new("RGB", (20, 20), "white").save(png)
        return {"ok": True, "png": str(png)}

    monkeypatch.setattr(render_module, "render_png", fake_render)
    result = render_panels(
        [{"path": str(a), "width": 321, "camera_preset": "+x",
          "merge_files": [str(local)]},
         {"path": str(a), "label": "second"}],
        tmp_path / "overrides.png", share_camera=False, share_color=False,
        merge_files=[str(shared)])
    assert result["ok"], result
    assert calls[0]["width"] == 321
    assert calls[0]["camera_preset"] == "+x"
    assert calls[0]["merge_files"] == [str(shared), str(local)]


# --------------------------------------------------------------------
# 5. compound selection
# --------------------------------------------------------------------

def test_select_compound_expression_and_carry(tmp_path):
    a = _box_field(tmp_path, "a", scale=1.0)
    res = select(a, "f > 1.5 and x > 0.4", out_file=tmp_path / "s.msh",
                 carry="f", extract=False)
    assert res["ok"], res
    assert 0 < res["n_selected"] < res["n_elements"]
    assert "f" in res["available_names"] and "v0" in res["available_names"]
    assert res["carried_view"] == "f"
    lo, hi = res["carried_range"]
    assert lo > 1.5 - 1e-9 and hi <= 2.0 + 1e-9
    # the mask view really is in the output
    text = (tmp_path / "s.msh").read_text(encoding="utf-8")
    assert '"selection"' in text and '"selection_f"' in text


def test_select_all_and_none_are_consistent(tmp_path):
    a = _box_field(tmp_path, "a")
    every = select(a, "True", out_file=tmp_path / "all.msh", extract=False)
    none = select(a, "False", out_file=tmp_path / "none.msh", extract=False)
    assert every["n_selected"] == every["n_elements"]
    assert every["fraction"] == pytest.approx(1.0)
    assert none["n_selected"] == 0
    # a coordinate half-space takes roughly half of a symmetric box
    half = select(a, "x > 0.5", out_file=tmp_path / "half.msh",
                  extract=False)
    assert 0.2 < half["fraction"] < 0.8


def test_select_unknown_name_fails_loudly(tmp_path):
    a = _box_field(tmp_path, "a")
    bad = select(a, "nosuchfield > 1", extract=False)
    assert not bad["ok"]
    assert "not defined" in bad["error"]
    assert "f" in bad["available_names"]


def test_select_rejects_python_object_access(tmp_path):
    a = _box_field(tmp_path, "a")
    bad = select(a, "(1).__class__.__mro__", extract=False)
    assert not bad["ok"]
    assert "Attribute" in bad["error"]


def test_select_keeps_coordinate_names_when_a_view_collides(tmp_path):
    path = _write_msh(
        tmp_path / "collision.msh",
        {1: [0, 0, 0], 2: [1, 0, 0], 3: [0, 1, 0], 4: [0, 0, 1]},
        {1: [1, 2, 3, 4]},
        [("x", 1, {1: [1], 2: [1], 3: [1], 4: [1]})])
    result = select(path, "x < 0.8 and x_view > 0", extract=False)
    assert result["ok"], result
    assert "x_view" in result["available_names"]
    assert result["n_selected"] == 1


# --------------------------------------------------------------------
# 1. pseudo-volume rendering
# --------------------------------------------------------------------

def test_volume_render_stacks_the_requested_slices(tmp_path):
    a = _box_field(tmp_path, "a", n=6)
    res = volume_render(a, tmp_path / "v.png", view="f", n_slices=9,
                        axis="x", alpha=0.4, alpha_power=2.0,
                        color={"range": [1.0, 2.0]},
                        width=280, height=260)
    assert res.get("ok"), res
    assert res["n_slices"] == 9
    assert "NOT ray-cast" in res["method"]
    assert res["color_range"] == [1.0, 2.0]
    assert not res.get("blank_check", {}).get("looks_blank", False)


def test_volume_render_rejects_bad_arguments(tmp_path):
    a = _box_field(tmp_path, "a")
    with pytest.raises(ValueError, match="axis must be"):
        volume_render(a, tmp_path / "v.png", axis="w")
    with pytest.raises(ValueError, match="n_slices"):
        volume_render(a, tmp_path / "v.png", n_slices=1)
    with pytest.raises(ValueError, match="alpha"):
        volume_render(a, tmp_path / "v.png", alpha=0.0)
    with pytest.raises(ValueError, match="alpha_power"):
        volume_render(a, tmp_path / "v.png", alpha_power=-1.0)
    with pytest.raises(ValueError, match="integer"):
        volume_render(a, tmp_path / "v.png", n_slices=2.5)
    vector = volume_render(a, tmp_path / "v.png", view="vec")
    assert not vector["ok"] and "scalar view" in vector["error"]


# --------------------------------------------------------------------
# 2. flow texture (LIC alternative)
# --------------------------------------------------------------------

def test_flow_texture_density_sets_the_spacing(tmp_path):
    a = _box_field(tmp_path, "a", n=6)
    coarse = flow_texture(a, view="vec", plane="xy", density=8,
                          out_file=tmp_path / "c.pos")
    fine = flow_texture(a, view="vec", plane="xy", density=32,
                        out_file=tmp_path / "f.pos")
    assert coarse.get("ok") and fine.get("ok"), (coarse, fine)
    assert fine["d_sep"] == pytest.approx(coarse["d_sep"] / 4.0)
    assert "NOT line integral convolution" in fine["method"]
    # a finer texture must place strictly more curves
    assert fine["n_lines"] > coarse["n_lines"]


def test_flow_texture_rejects_bad_plane_and_density(tmp_path):
    a = _box_field(tmp_path, "a")
    with pytest.raises(ValueError, match="unknown plane"):
        flow_texture(a, plane="ab")
    with pytest.raises(ValueError, match="density must be"):
        flow_texture(a, density=0)
    scalar = flow_texture(a, view="f")
    assert not scalar["ok"] and "3-component" in scalar["error"]
