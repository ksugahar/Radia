"""Cross-server probe contract: build123d_probe (CAD side) and
cubit_probe (mesh side) must emit the same core per-entity vocabulary
(PROBE_SOLID_CORE_KEYS / PROBE_FACE_CORE_KEYS) so per-body numbers can
be compared across the STEP -> Cubit handoff.

Also locks the assembly-Location regression: a solid placed at
Pos(3,0,0) must report its WORLD centroid, not local coordinates
(bug found by the 2026-08-05 compat E2E).
"""

import json

import pytest

from radia_mcp.common.server_hardening import (
    PROBE_FACE_CORE_KEYS,
    PROBE_SOLID_CORE_KEYS,
)


# ---------------------------------------------------------------------------
# Cubit side: probe_ops output carries the core keys (mock cubit)
# ---------------------------------------------------------------------------

def test_cubit_probe_entities_carries_core_keys():
    from radia_mcp.cubit import probe_ops
    from test_cubit_probe_ops import _MockCubit

    ent = probe_ops.op_probe(_MockCubit(), ["entities"])
    assert PROBE_SOLID_CORE_KEYS <= set(ent["volumes"][0])
    assert PROBE_FACE_CORE_KEYS <= set(ent["surfaces"][0])


# ---------------------------------------------------------------------------
# build123d side: real geometry (skipped in minimal-dep CI)
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def labeled_assembly_step(tmp_path_factory):
    b3d = pytest.importorskip("build123d")
    from build123d import Box, Compound, Cylinder, Pos, export_step

    box = Box(1, 1, 2)
    box.label = "yoke"
    cyl = Pos(3, 0, 0) * Cylinder(0.4, 1.0)
    cyl.label = "core"
    p = tmp_path_factory.mktemp("compat") / "asm.step"
    export_step(Compound(children=[box, cyl]), str(p))
    return p


def test_b3d_probe_entities_core_keys_and_world_placement(
        labeled_assembly_step):
    from radia_mcp.build123d.server import build123d_probe

    out = json.loads(build123d_probe(str(labeled_assembly_step),
                                     query="entities"))
    assert out["status"] == "ok"
    assert out["solid_count"] == 2
    for s in out["solids"]:
        assert PROBE_SOLID_CORE_KEYS <= set(s), s.keys()
        assert "label" in s          # build123d EXTENDS the core schema
    for f in out["faces"]:
        assert PROBE_FACE_CORE_KEYS <= set(f), f.keys()

    by_label = {s["label"]: s for s in out["solids"]}
    assert set(by_label) == {"yoke", "core"}
    # Location regression: 'core' sits at x=3 in WORLD coordinates.
    assert abs(by_label["core"]["centroid"][0] - 3.0) < 1e-6
    assert abs(by_label["core"]["volume"] - 0.502655) < 1e-4
    assert by_label["yoke"]["extent"] == [1.0, 1.0, 2.0]
    # FACES are world-placed too (verified 2026-08-05: shape.faces() on
    # the top-level compound applies child Locations): the cylinder at
    # x=3 contributes exactly 3 faces centred at x > 2.
    assert sum(1 for f in out["faces"] if f["center"][0] > 2.0) == 3


def test_b3d_probe_labels_audit(labeled_assembly_step):
    from radia_mcp.build123d.server import build123d_probe

    out = json.loads(build123d_probe(str(labeled_assembly_step),
                                     query="labels"))
    assert out["audit"]["passed"] is True
    assert [s["label"] for s in out["solids"]] in (
        [["yoke", "core"], ["core", "yoke"]][0],
        [["yoke", "core"], ["core", "yoke"]][1],
    )


def test_b3d_probe_unlabeled_solid_warns(tmp_path):
    pytest.importorskip("build123d")
    from build123d import Box, export_step
    from radia_mcp.build123d.server import build123d_probe

    p = tmp_path / "plain.step"
    export_step(Box(1, 1, 1), str(p))
    out = json.loads(build123d_probe(str(p), query="labels"))
    assert out["audit"]["passed"] is True     # warnings only
    assert any("unnamed" in w for w in out["audit"]["warnings"])


def test_b3d_labeled_solids_deep_nesting(tmp_path):
    """3-level assembly label inheritance, locked from measured behavior
    (2026-08-05, identical in-memory and after STEP round trip):

    * a part keeps its own label;
    * an unlabeled part inside a labeled sub-assembly inherits the
      SUB-ASSEMBLY label (nearest labeled ancestor);
    * parts under an unlabeled sub-assembly inherit the root label;
    * world placement is preserved at every level.
    """
    pytest.importorskip("build123d")
    from build123d import (Box, Compound, Cylinder, Pos, Sphere,
                           export_step, import_step)
    from radia_mcp.build123d.server import _labeled_solids

    m1 = Box(1, 1, 1)
    m1.label = "magnet_a"
    m2 = Pos(2, 0, 0) * Box(1, 1, 1)              # unlabeled
    sub1 = Compound(children=[m1, m2], label="rotor")
    c1 = Pos(0, 4, 0) * Cylinder(0.5, 1)
    c1.label = "pole"
    c2 = Pos(2, 4, 0) * Cylinder(0.5, 1)          # unlabeled
    sub2 = Compound(children=[c1, c2])            # unlabeled sub-assembly
    s = Pos(0, 8, 0) * Sphere(0.6)
    s.label = "sensor"
    root = Compound(children=[sub1, sub2, s], label="machine")

    def rows(shape):
        return sorted(
            (lab, round(sol.center().X, 1), round(sol.center().Y, 1))
            for lab, sol in _labeled_solids(shape))

    expected = sorted([
        ("magnet_a", 0.0, 0.0),
        ("rotor", 2.0, 0.0),          # inherited from labeled sub-asm
        ("pole", -0.0, 4.0),
        ("machine", 2.0, 4.0),        # inherited from root
        ("sensor", -0.0, 8.0),
    ])
    assert rows(root) == expected

    p = tmp_path / "nested3.step"
    export_step(root, str(p))
    assert rows(import_step(str(p))) == expected


def test_b3d_probe_summary_and_errors(tmp_path):
    pytest.importorskip("build123d")
    from build123d import Box, export_step
    from radia_mcp.build123d.server import build123d_probe

    p = tmp_path / "b.step"
    export_step(Box(1, 1, 1), str(p))
    smy = json.loads(build123d_probe(str(p), query="summary"))
    assert smy["solids"] == 1 and smy["faces"] == 6

    missing = json.loads(build123d_probe(str(tmp_path / "nope.step")))
    assert missing["status"] == "error" and missing["kind"] == "input"

    unknown = json.loads(build123d_probe(str(p), query="nope"))
    assert unknown["status"] == "error"
    assert "Unknown probe query" in unknown["error"]
