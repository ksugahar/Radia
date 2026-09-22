"""A recorded hash that is recomputed on every edit records nothing.

Four stored result artifacts were edited after their runs, to take a product
name out of a formulation label.  The certificate that indexes them calls
``result_sha256`` a content identity, so that edit had to be made auditable
rather than absorbed: the pre-rename hash stays recorded, and putting the old
label back has to reproduce it exactly.

This performs that audit from the repository alone.  It is the whole point of
keeping the old label in the record, and it fails if a future edit to those
artifacts changes anything besides the label -- including an edit that then
updates the hashes to match itself.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CERTIFICATES = sorted(
    (ROOT / "validation_test").rglob("*_certificate.json"))


def _with_rename():
    for path in CERTIFICATES:
        data = json.loads(path.read_text(encoding="utf-8"))
        if "label_rename" in data:
            yield path, data


@pytest.mark.parametrize(
    "case", list(_with_rename()),
    ids=lambda c: c[0].name if isinstance(c, tuple) else str(c))
def test_restoring_the_old_label_reproduces_the_recorded_hash(case):
    """The audit the record promises, actually carried out."""
    path, data = case
    rename = data["label_rename"]
    old, new = rename["from"].encode(), rename["to"].encode()
    assert old != new
    assert rename["artifacts"], "a rename record with no artifacts proves nothing"

    for row in rename["artifacts"]:
        artifact = path.parent / row["result"]
        assert artifact.exists(), f"{row['result']} is recorded but missing"
        current = artifact.read_bytes()

        # As it stands now.
        assert hashlib.sha256(current).hexdigest() == row["result_sha256"], (
            f"{row['result']} no longer matches its current recorded hash")
        # And as the run produced it.
        restored = current.replace(new, old)
        assert restored != current, (
            f"{row['result']} does not contain the renamed label at all, so "
            f"the record does not describe it")
        assert (hashlib.sha256(restored).hexdigest()
                == row["pre_rename_result_sha256"]), (
            f"{row['result']} differs from the artifact the run produced by "
            f"more than the label; the rename record is no longer true")


@pytest.mark.parametrize(
    "case", list(_with_rename()),
    ids=lambda c: c[0].name if isinstance(c, tuple) else str(c))
def test_the_level_rows_agree_with_the_rename_record(case):
    """The index and the record cannot drift apart silently."""
    path, data = case
    recorded = {row["result"]: row["result_sha256"]
                for row in data["label_rename"]["artifacts"]}
    seen = 0
    for key in ("level_results", "convergence_levels"):
        for row in data.get(key, []):
            if not isinstance(row, dict):
                continue
            name = row.get("result")
            if name in recorded:
                assert row["result_sha256"] == recorded[name], (
                    f"{name}: the level row and the rename record disagree")
                seen += 1
    assert seen == len(recorded), (
        "the rename record names artifacts the certificate does not index")


def test_a_renamed_label_carries_no_product_name():
    """The rename has to have removed something, not merely moved it."""
    for path, data in _with_rename():
        new = data["label_rename"]["to"]
        assert "TOSCA" not in new.upper(), path.name
