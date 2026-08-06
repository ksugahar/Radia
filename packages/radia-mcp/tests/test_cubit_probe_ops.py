"""Unit tests for the SHARED probe implementation (probe_ops.op_probe)
used by both Cubit session runners, exercised against a mock cubit
module (no Cubit license needed)."""

from radia_mcp.cubit import probe_ops


class _MockVolume:
    def __init__(self, vid):
        self._id = vid

    def centroid(self):
        return (1.0, 2.0, 3.0)

    def bounding_box(self):
        # (min_x, min_y, min_z, max_x, max_y, max_z) — verified layout
        return (-0.5, -0.5, -0.5, 0.5, 0.5, 1.5)

    def volume(self):
        return 2.0


class _MockSurface:
    def center_point(self):
        return (0.0, 0.0, 0.5)

    def bounding_box(self):
        return (-0.5, -0.5, 0.5, 0.5, 0.5, 0.5)

    def area(self):
        return 1.0


class _MockCubit:
    """Minimal cubit-module stand-in for probe queries."""

    def get_volume_count(self):
        return 1

    def get_surface_count(self):
        return 1

    def get_curve_count(self):
        return 0

    def get_vertex_count(self):
        return 0

    def get_node_count(self):
        return 8

    def get_hex_count(self):
        return 1

    def get_tet_count(self):
        return 0

    def parse_cubit_list(self, kind, expr):
        if kind == "volume":
            return (1,)
        if kind == "surface":
            return (1,)
        if kind == "hex":
            return (1,)
        return ()

    def get_quality_values(self, kind, ids, metric):
        return (0.9,) * len(ids)

    def get_volume_meshing_scheme(self, vid):
        return "sweep"

    def volume(self, vid):
        return _MockVolume(vid)

    def surface(self, sid):
        return _MockSurface()

    def get_block_id_list(self):
        return (1,)

    def get_exodus_entity_name(self, kind, eid):
        return "iron"

    def get_block_hexes(self, bid):
        # Real Cubit 2025.12: DIRECT element membership is empty for a
        # `block N add volume V` block; elements come via the volume.
        return ()

    def get_block_tets(self, bid):
        return ()

    def get_block_wedges(self, bid):
        return ()

    def get_block_pyramids(self, bid):
        return ()

    def get_block_tris(self, bid):
        return ()

    def get_block_faces(self, bid):
        return ()

    def get_block_volumes(self, bid):
        return (1,)

    def get_block_surfaces(self, bid):
        return ()

    def get_sideset_id_list(self):
        return ()

    def get_sideset_surfaces(self, sid):
        return ()


def test_summary_and_counts():
    cub = _MockCubit()
    assert probe_ops.op_probe(cub, ["volume_count"]) == 1
    smy = probe_ops.op_probe(cub, ["summary"])
    assert smy["volumes"] == 1 and smy["hexes"] == 1


def test_entities_extent_derived_from_minmax():
    ent = probe_ops.op_probe(_MockCubit(), ["entities"])
    v = ent["volumes"][0]
    assert v["bbox_min"] == [-0.5, -0.5, -0.5]
    assert v["bbox_max"] == [0.5, 0.5, 1.5]
    assert v["extent"] == [1.0, 1.0, 2.0]     # max - min, NOT bb[3:6]
    s = ent["surfaces"][0]
    assert s["extent"] == [1.0, 1.0, 0.0]     # flat face: zero normal extent


def test_labels_membership_and_audit():
    lab = probe_ops.op_probe(_MockCubit(), ["labels"])
    b = lab["blocks"][0]
    assert b["name"] == "iron" and b["volume_elems"] == 1
    assert b["surface_elems"] == 0
    assert lab["audit"]["passed"] is True


def test_per_volume_rows():
    rows = probe_ops.op_probe(_MockCubit(), ["per_volume"])
    assert rows == [{"id": 1, "scheme": "sweep", "hex": 1, "tet": 0,
                     "meshed": True}]


def test_unknown_query_lists_valid_queries():
    out = probe_ops.op_probe(_MockCubit(), ["nope"])
    assert "Unknown probe query" in out["error"]
    for q in ("summary", "quality", "per_volume", "entities", "labels"):
        assert q in out["error"]


def test_trimmed_traceback_drops_runner_frames():
    def user_code():
        raise ValueError("boom")

    try:
        user_code()
    except ValueError:
        text = probe_ops.trimmed_traceback(
            exclude_basenames=("test_cubit_probe_ops.py",))
        # All frames are in this test file -> trimming would leave
        # nothing -> falls back to the FULL traceback.
        assert "boom" in text and "user_code" in text
        text2 = probe_ops.trimmed_traceback(exclude_basenames=("other.py",))
        assert "user_code" in text2 and "boom" in text2
