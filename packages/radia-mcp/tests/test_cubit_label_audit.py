"""Tests for the cubit-free block/sideset label audit shared by the
Cubit daemon probe("labels") and the MCP server."""

from radia_mcp.cubit.label_audit import audit_label_records, is_strict_label


def _block(bid, name, volume_elems=0, surface_elems=0):
    return {
        "id": bid,
        "name": name,
        "volume_elems": volume_elems,
        "surface_elems": surface_elems,
    }


def test_clean_convention_passes():
    blocks = [
        _block(1, "iron", volume_elems=100),
        _block(2, "air", volume_elems=200),
    ]
    sidesets = [{"id": 1, "name": "source"}, {"id": 2, "name": "kelvin_ext"}]
    audit = audit_label_records(blocks, sidesets)
    assert audit["passed"] is True
    assert audit["errors"] == []
    assert audit["warnings"] == []


def test_mixed_block_is_error():
    blocks = [_block(1, "mixed", volume_elems=50, surface_elems=10)]
    audit = audit_label_records(blocks, [])
    assert audit["passed"] is False
    assert len(audit["errors"]) == 1
    assert "LOST" in audit["errors"][0]


def test_unnamed_block_and_sideset_warn():
    audit = audit_label_records([_block(3, "", volume_elems=1)],
                                [{"id": 7, "name": ""}])
    assert audit["passed"] is True  # warnings only
    assert len(audit["warnings"]) == 2
    assert all("unnamed" in w for w in audit["warnings"])


def test_casefold_collision_is_error():
    blocks = [
        _block(1, "iron", volume_elems=1),
        _block(2, "Iron", volume_elems=1),
    ]
    audit = audit_label_records(blocks, [])
    assert audit["passed"] is False
    assert any("casefold collision" in e for e in audit["errors"])


def test_non_snake_case_warns_but_passes():
    audit = audit_label_records([_block(1, "Iron_Yoke", volume_elems=1)], [])
    assert audit["passed"] is True
    assert any("strict naming" in w for w in audit["warnings"])


def test_strict_label_rules():
    assert is_strict_label("iron")
    assert is_strict_label("kelvin_ext")
    assert is_strict_label("sym_bn=0_x")
    assert is_strict_label("sym_ht=0_z")
    assert is_strict_label("GND")
    assert not is_strict_label("Iron")
    assert not is_strict_label("volume 1")
    assert not is_strict_label("gnd_")
