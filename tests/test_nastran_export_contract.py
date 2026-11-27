"""Contract of radia.nastran_export: refusals and card layout.

Guards: order-2 export refuses a straight mesh, a discontinuous element map,
an unsupported order and an existing output; cards stay in Nastran small/large
fields with '+'/'*' continuations; every element node is a written GRID.
Assumes: nothing about numerical accuracy; the Gmsh reader check in
validation_test/gmsh_post/test_nastran_export_gmsh_reader.py owns that.
"""
import re

import ngsolve as ng
import pytest
from netgen.occ import Box, OCCGeometry, Pnt

from radia.nastran_export import main, write_nastran_bdf

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")


def _box(curve=2):
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.5))
    if curve > 1:
        mesh.Curve(curve)
    return mesh


def _cards(text):
    """(name, [int fields]) per element card, continuation lines joined."""
    cards, cur = [], None
    for line in text.splitlines():
        assert len(line) <= 80, line
        if line.startswith("+"):
            cur[1].extend(int(line[i:i + 8]) for i in range(8, len(line.rstrip("+").rstrip()), 8))
        elif line.startswith("C"):
            body = line[:72]
            cur = (line[:8].strip(), [int(body[i:i + 8]) for i in range(8, len(body.rstrip()), 8)])
            cards.append(cur)
    return cards


def test_second_order_cards_reference_written_grids(tmp_path):
    out = tmp_path / "box.bdf"
    info = write_nastran_bdf(_box(), out, order=2)
    text = out.read_text()
    grids = {int(m) for m in re.findall(r"^GRID\*\s+(\d+)", text, re.M)}
    assert len(grids) == info["n_grid"] == info["n_vertices"] + info["n_mid_edge"]
    sizes = {"CTETRA": 10, "CTRIA6": 6}
    cards = _cards(text)
    assert {n for n, _ in cards} == set(sizes), {n for n, _ in cards}
    for name, fields in cards:
        assert len(fields) == 2 + sizes[name], (name, fields)
        assert set(fields[2:]) <= grids, (name, fields)
    assert sorted(f[0] for _, f in cards) == list(range(1, len(cards) + 1))
    assert text.rstrip().endswith("ENDDATA")


def test_order_two_needs_a_curved_mesh(tmp_path):
    with pytest.raises(ValueError, match="curve order 1"):
        write_nastran_bdf(_box(curve=1), tmp_path / "x.bdf", order=2)
    assert not (tmp_path / "x.bdf").exists()


def test_unsupported_order_and_existing_output(tmp_path):
    mesh = _box()
    with pytest.raises(ValueError, match="order must be 1 or 2"):
        write_nastran_bdf(mesh, tmp_path / "x.bdf", order=3)
    out = tmp_path / "y.bdf"
    out.write_text("keep")
    with pytest.raises(FileExistsError):
        write_nastran_bdf(mesh, out, order=2)
    assert out.read_text() == "keep"
    write_nastran_bdf(mesh, out, order=2, overwrite=True)
    assert out.read_text().startswith("$ Nastran BDF")


def test_deformed_mesh_is_refused(tmp_path):
    """The GRID points are the mesh vertices; a deformed map leaves them."""
    mesh = _box()
    gf = ng.GridFunction(ng.VectorH1(mesh, order=1))
    gf.Set(ng.CF((0, 0, 0.1 * ng.x)))
    mesh.SetDeformation(gf)
    try:
        with pytest.raises(ValueError, match="does not pass through its vertex"):
            write_nastran_bdf(mesh, tmp_path / "deformed.bdf", order=2)
    finally:
        mesh.UnsetDeformation()
    assert not (tmp_path / "deformed.bdf").exists()


def test_torn_shared_edge_is_measured():
    from radia.nastran_export import _merge_mid_edges
    a = ([0, 1, 2], ((0, 1), (1, 2)), [(0.5, 0, 0), (1, 0.5, 0)])
    b = ([1, 3, 0], ((0, 2), (0, 1)), [(0.5, 0, 0.002), (2, 0, 0)])
    nodes, mismatch = _merge_mid_edges([a, b])
    assert sorted(nodes) == [(0, 1), (1, 2), (1, 3)]
    assert mismatch == pytest.approx(0.002)


def test_torn_mesh_is_refused(tmp_path, monkeypatch):
    import radia.nastran_export as ne

    real = ne._merge_mid_edges

    def torn(entries):
        nodes, _ = real(entries)
        return nodes, 1e-3
    monkeypatch.setattr(ne, "_merge_mid_edges", torn)
    with pytest.raises(ValueError, match="shared mid-edge nodes disagree"):
        write_nastran_bdf(_box(), tmp_path / "torn.bdf", order=2)
    assert not (tmp_path / "torn.bdf").exists()


def test_cli_reports_refusals_with_exit_2(tmp_path, capsys):
    vol = tmp_path / "straight.vol"
    _box(curve=1).ngmesh.Save(str(vol))
    assert main([str(vol), str(tmp_path / "o.bdf"), "--order", "2"]) == 2
    assert "curve order 1" in capsys.readouterr().err
    assert main([str(vol), str(tmp_path / "o.bdf")]) == 0


def test_a_file_that_appears_late_is_not_replaced(tmp_path, monkeypatch):
    """The output is created exclusively, not checked once and replaced."""
    import radia.nastran_export as ne

    out = tmp_path / "late.bdf"
    real = ne.os.link

    def publish_after_competitor(stage, target):
        assert stage.parent == target.parent
        assert stage.read_bytes().endswith(b"ENDDATA\n")
        assert not target.exists()
        target.write_text("someone else")
        real(stage, target)

    monkeypatch.setattr(ne.os, "link", publish_after_competitor)
    with pytest.raises(FileExistsError):
        write_nastran_bdf(_box(), out, order=2)
    assert out.read_text() == "someone else"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["late.bdf"]


def test_labels_are_written_as_ascii_comments(tmp_path):
    from netgen.occ import Box, OCCGeometry, Pnt
    box = Box(Pnt(0, 0, 0), Pnt(1, 1, 1))
    box.mat("鉄心\nPSOLID")
    mesh = ng.Mesh(OCCGeometry(box).GenerateMesh(maxh=0.6))
    out = tmp_path / "u.bdf"
    write_nastran_bdf(mesh, out, order=1)
    text = out.read_bytes().decode("ascii")
    assert r"$ material \u9244\u5fc3\u000aPSOLID" in text, text[-300:]
    assert sum(line.startswith("PSOLID") for line in text.splitlines()) == 1
    assert sorted(p.name for p in tmp_path.iterdir()) == ["u.bdf"]


def test_failed_overwrite_keeps_the_old_file_and_no_stage(tmp_path, monkeypatch):
    import radia.nastran_export as ne

    out = tmp_path / "keep.bdf"
    out.write_text("old")

    def fail(*_args):
        raise OSError("replace failed")
    monkeypatch.setattr(ne.os, "replace", fail)
    with pytest.raises(OSError, match="replace failed"):
        write_nastran_bdf(_box(), out, order=2, overwrite=True)
    assert out.read_text() == "old"
    assert sorted(p.name for p in tmp_path.iterdir()) == ["keep.bdf"]


@pytest.mark.parametrize("overwrite", [False, True])
def test_cleanup_failure_preserves_publication_error(tmp_path, monkeypatch, overwrite):
    import errno
    from pathlib import Path

    import radia.nastran_export as ne

    out = tmp_path / "output.bdf"
    if overwrite:
        out.write_bytes(b"original")
    primary = OSError(errno.ENOSPC, "publication failed")

    def fail_publish(stage, target):
        assert stage.read_bytes() == b"complete deck"
        assert target == out
        raise primary

    def fail_cleanup(self, **kwargs):
        raise PermissionError("stage cleanup denied")

    monkeypatch.setattr(ne.os, "replace" if overwrite else "link", fail_publish)
    monkeypatch.setattr(Path, "unlink", fail_cleanup)
    with pytest.raises(OSError) as raised:
        ne._publish(out, b"complete deck", overwrite)

    assert raised.value is primary
    stage, = tmp_path.glob("*.stage")
    assert str(stage) in raised.value.__notes__[0]
    assert "stage cleanup denied" in raised.value.__notes__[0]
    if overwrite:
        assert out.read_bytes() == b"original"
    else:
        assert not out.exists()


@pytest.mark.parametrize("close_fails", [False, True])
def test_partial_write_only_touches_stage(tmp_path, monkeypatch, close_fails):
    import errno

    import radia.nastran_export as ne

    out = tmp_path / "output.bdf"
    real = ne.tempfile.NamedTemporaryFile
    primary = OSError(errno.ENOSPC, "injected short write")

    class FailingStream:
        def __init__(self, *args, **kwargs):
            self.stream = real(*args, **kwargs)
            self.name = self.stream.name

        def close(self):
            self.stream.close()
            if close_fails:
                raise OSError("stage close failed")

        def write(self, data):
            self.stream.write(data[:5])
            self.stream.flush()
            assert not out.exists()
            raise primary

    monkeypatch.setattr(ne.tempfile, "NamedTemporaryFile", FailingStream)
    with pytest.raises(OSError) as raised:
        ne._publish(out, b"complete deck", False)

    assert raised.value is primary
    if close_fails:
        assert "stage close failed" in raised.value.__notes__[0]
    assert list(tmp_path.iterdir()) == []
