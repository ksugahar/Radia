from zipfile import ZipFile

from radia_mcp.presentation import tools


def _write_notes_pptx(path, payload: bytes) -> None:
    with ZipFile(path, "w") as archive:
        archive.writestr("ppt/notesSlides/notesSlide1.xml", payload)


def test_presentation_check_pptx_notes_encoding_accepts_utf8_and_bom(tmp_path):
    pptx = tmp_path / "notes.pptx"
    notes = tmp_path / "notes.md"
    _write_notes_pptx(pptx, "<notes>日本語の発表ノート</notes>".encode("utf-8"))
    notes.write_bytes(b"\xef\xbb\xbf" + "発表セリフ".encode("utf-8"))

    report = tools.presentation_check_pptx_notes_encoding(str(pptx), str(notes))

    assert report["passed"] is True
    assert report["notes_xml_checked"] == 1
    assert report["exported_notes"]["utf8_bom"] is True


def test_presentation_check_pptx_notes_encoding_flags_missing_bom(tmp_path):
    pptx = tmp_path / "notes.pptx"
    notes = tmp_path / "notes.md"
    _write_notes_pptx(pptx, "<notes>日本語</notes>".encode("utf-8"))
    notes.write_text("発表セリフ", encoding="utf-8")

    report = tools.presentation_check_pptx_notes_encoding(str(pptx), str(notes))

    assert report["passed"] is False
    assert any(item["issue"] == "missing_utf8_bom" for item in report["findings"])
