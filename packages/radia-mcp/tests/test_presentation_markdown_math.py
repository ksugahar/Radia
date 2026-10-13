"""Real PPTX packaging with a fake equation backend; no font/Office activity."""
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys
import types
import zipfile
import pytest

pytest.importorskip("pptx")
from pptx import Presentation
from lxml import etree

SPEC = spec_from_file_location("markdown_math_contract", Path(__file__).parents[1] /
    "src/radia_mcp/presentation/_markdown_equations.py")
assembly = module_from_spec(SPEC)
SPEC.loader.exec_module(assembly)


@pytest.fixture
def backend(monkeypatch):
    calls = []
    def export(tex, output, timeout_s=30):
        calls.append(tex)
        root = etree.Element(f"{{{assembly.MATH_NS}}}oMathPara", nsmap={"m": assembly.MATH_NS, "a": assembly.DRAWING_NS})
        math = etree.SubElement(root, f"{{{assembly.MATH_NS}}}oMath")
        run = etree.SubElement(math, f"{{{assembly.MATH_NS}}}r")
        etree.SubElement(run, f"{{{assembly.DRAWING_NS}}}t").text = tex
        Path(output).write_bytes(etree.tostring(root))
        return Path(output)
    fake = types.SimpleNamespace(__version__="3.1.4", export_powerpoint_math=export)
    monkeypatch.setitem(sys.modules, "eqnedit64", fake)
    return calls, fake


def test_inline_display_notes_and_code_survive_pptx_roundtrip(tmp_path, backend):
    calls, _ = backend
    source = "# Slide\n\nBefore $x$ after.\n\n$$x^2$$\n\n> Say $x$ and\n> $$y$$\n\n- Item $x$\n\n```text\n$not_math$\n```\n"
    path = tmp_path / "deck.pptx"
    result = assembly.presentation_markdown_to_pptx(source, str(path), title="Deck")
    assert result["body_math"] == 3 and result["notes_math"] == 2
    assert result["slides"] == 2 and calls == ["x", "x^2", "y"]
    reopened = Presentation(path)
    assert len(reopened.slides) == 2
    again = tmp_path / "roundtrip.pptx"
    reopened.save(again)
    for candidate in (path, again):
        with zipfile.ZipFile(candidate) as archive:
            slide = etree.fromstring(archive.read("ppt/slides/slide2.xml"))
            assert len(list(slide.iter(f"{{{assembly.A14_NS}}}m"))) == 3
            assert len(list(slide.iter(f"{{{assembly.MC_NS}}}AlternateContent"))) == 1
            assert not list(slide.iter(f"{{{assembly.MC_NS}}}Fallback"))
            text = "".join(slide.itertext())
            assert "Before " in text and " after." in text and "$not_math$" in text
            notes_name = next(n for n in archive.namelist() if n.startswith("ppt/notesSlides/notesSlide") and n.endswith(".xml"))
            notes = etree.fromstring(archive.read(notes_name))
            assert len(list(notes.iter(f"{{{assembly.A14_NS}}}m"))) == 2


@pytest.mark.parametrize("source", ["# x\n\n$x", "# x\n\n| a | b |", "# x\n\n**bold**", "# x\n\n- a\n  - b", "```\nno close", "# $x$", "# x\n\n[x](url)"])
def test_unsupported_markdown_raises_and_preserves_output(tmp_path, backend, source):
    path = tmp_path / "deck.pptx"
    path.write_bytes(b"previous")
    with pytest.raises(ValueError):
        assembly.presentation_markdown_to_pptx(source, str(path))
    assert path.read_bytes() == b"previous"


@pytest.mark.parametrize("version", ["3.1.3", "", "invalid"])
def test_old_or_unknown_backend_fails_loudly(tmp_path, backend, version):
    _, fake = backend
    fake.__version__ = version
    with pytest.raises(RuntimeError, match="required"):
        assembly.presentation_markdown_to_pptx("# x\n$x$", str(tmp_path / "deck.pptx"))


def test_missing_backend_fails_without_images(tmp_path, monkeypatch):
    monkeypatch.setitem(sys.modules, "eqnedit64", None)
    with pytest.raises(RuntimeError, match="Install eqnedit64"):
        assembly.presentation_markdown_to_pptx("# x\n$x$", str(tmp_path / "deck.pptx"))


def test_plain_code_escapes_and_backend_error_are_not_silently_flattened(tmp_path, backend):
    calls, fake = backend
    assert list(assembly._pieces(r"Price \$5; `**$code$**`")) == [("text", "Price $5; "), ("literal", "**$code$**")]
    def fail(*args, **kwargs):
        raise RuntimeError("unsupported TeX")
    fake.export_powerpoint_math = fail
    with pytest.raises(RuntimeError, match="unsupported TeX"):
        assembly.presentation_markdown_to_pptx("# x\n$bad$", str(tmp_path / "deck.pptx"))
    assert not (tmp_path / "deck.pptx").exists()
