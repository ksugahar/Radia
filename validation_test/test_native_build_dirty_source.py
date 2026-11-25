"""Validate provenance rejection when an already dirty checkout changes."""

from tests.test_native_build_provenance import _assert_change_rejected, _quote, _repo


def test_dirty_source_change_during_build_is_still_rejected(tmp_path):
    repo, binary = _repo(tmp_path)
    (repo / "source.txt").write_text("dirty at start\n", encoding="utf-8")
    mutation = f"[IO.File]::WriteAllText({_quote(repo / 'source.txt')}, 'dirtier')"
    _assert_change_rejected(repo, binary, mutation)
