from __future__ import annotations

import json
from pathlib import Path

from PIL import Image
from pptx import Presentation
from pptx.util import Inches

from radia_mcp.presentation.tools import (
    presentation_check_embedded_figure_text_size,
)


def _deck_with_picture(tmp_path: Path) -> tuple[Path, str]:
    image_path = tmp_path / "figure.png"
    Image.new("RGB", (1000, 500), "white").save(image_path)
    prs = Presentation()
    slide = prs.slides.add_slide(prs.slide_layouts[6])
    picture = slide.shapes.add_picture(
        str(image_path), Inches(1), Inches(1), width=Inches(10)
    )
    path = tmp_path / "deck.pptx"
    prs.save(path)
    return path, picture.name


def _write_manifest(
        tmp_path: Path, shape_name: str, bbox_height_px: int) -> Path:
    path = tmp_path / "ocr.json"
    payload = {
        "pictures": [
            {
                "slide": 1,
                "shape": shape_name,
                "words": [
                    {
                        "text": "HACApK",
                        "confidence": 0.99,
                        "bbox": [100, 100, 300, 100 + bbox_height_px],
                    }
                ],
            }
        ]
    }
    path.write_text(json.dumps(payload), encoding="utf-8")
    return path


def test_embedded_figure_text_at_24pt_passes(tmp_path: Path) -> None:
    deck, shape_name = _deck_with_picture(tmp_path)
    manifest = _write_manifest(tmp_path, shape_name, bbox_height_px=24)

    result = presentation_check_embedded_figure_text_size(
        str(deck), ocr_backend="manifest", ocr_manifest_path=str(manifest)
    )

    assert result["passed"] is True
    assert result["violation_count"] == 0
    assert result["pictures"][0]["minimum_estimated_font_pt"] == 24.0


def test_embedded_figure_text_below_20pt_is_a_violation(
        tmp_path: Path) -> None:
    deck, shape_name = _deck_with_picture(tmp_path)
    manifest = _write_manifest(tmp_path, shape_name, bbox_height_px=14)

    result = presentation_check_embedded_figure_text_size(
        str(deck), ocr_backend="manifest", ocr_manifest_path=str(manifest)
    )

    assert result["passed"] is False
    assert result["violation_count"] == 1
    assert result["pictures"][0]["minimum_estimated_font_pt"] == 14.0


def test_source_size_evidence_uses_actual_paste_width(tmp_path: Path) -> None:
    deck, shape_name = _deck_with_picture(tmp_path)
    manifest = tmp_path / "source-evidence.json"
    manifest.write_text(json.dumps({
        "pictures": [{
            "slide": 1,
            "shape": shape_name,
            "source_evidence": {
                "minimum_source_font_pt": 24,
                "source_width_in": 12,
            },
        }]
    }), encoding="utf-8")

    result = presentation_check_embedded_figure_text_size(
        str(deck), ocr_backend="manifest", ocr_manifest_path=str(manifest)
    )

    assert result["passed"] is True
    assert result["pictures"][0]["evidence_type"] == "source_size"
    assert result["pictures"][0]["embed_scale"] == 0.8333
    assert result["pictures"][0]["minimum_estimated_font_pt"] == 20.0


def test_source_size_evidence_below_floor_fails(tmp_path: Path) -> None:
    deck, shape_name = _deck_with_picture(tmp_path)
    manifest = tmp_path / "source-evidence-small.json"
    manifest.write_text(json.dumps({
        "pictures": [{
            "slide": 1,
            "shape": shape_name,
            "source_evidence": {
                "minimum_source_font_pt": 18,
                "source_width_in": 10,
            },
        }]
    }), encoding="utf-8")

    result = presentation_check_embedded_figure_text_size(
        str(deck), ocr_backend="manifest", ocr_manifest_path=str(manifest)
    )

    assert result["passed"] is False
    assert result["pictures"][0]["minimum_estimated_font_pt"] == 18.0


def test_unverified_picture_does_not_pass(tmp_path: Path) -> None:
    deck, _ = _deck_with_picture(tmp_path)

    result = presentation_check_embedded_figure_text_size(str(deck))

    assert result["passed"] is False
    assert result["unresolved_count"] == 1


def test_visually_confirmed_textless_picture_is_exempt(tmp_path: Path) -> None:
    deck, shape_name = _deck_with_picture(tmp_path)

    result = presentation_check_embedded_figure_text_size(
        str(deck), confirmed_textless_shapes=[f"1:{shape_name}"]
    )

    assert result["passed"] is True
    assert result["pictures"][0]["status"] == "confirmed_textless"
