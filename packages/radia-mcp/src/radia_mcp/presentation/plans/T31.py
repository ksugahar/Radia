"""Presentation T31: TTS audio embedding for PowerPoint decks.

Generate per-slide MP3 narration from a Markdown script and embed it into
a PPTX using PowerPoint COM automation.  The first slide can be started by
click, while later slides can start with the previous animation, matching the
common "audio narration deck" workflow for conference backups.
"""

from __future__ import annotations

import asyncio
import os
import pathlib
import re
import shutil
import tempfile
from typing import Any


def _parse_slide_markdown(script_md_path: pathlib.Path) -> dict[int, str]:
    text = script_md_path.read_text(encoding="utf-8-sig")
    blocks = re.split(r"^##\s+Slide\s+(\d+)\b.*$", text, flags=re.M)
    notes: dict[int, str] = {}
    for i in range(1, len(blocks), 2):
        slide_no = int(blocks[i])
        body = _clean_tts_text(blocks[i + 1])
        if body:
            notes[slide_no] = body
    if notes:
        return notes

    # Fallback for TTS files written as "Slide N. Title." paragraphs.
    pattern = re.compile(r"^Slide\s+(\d+)\b.*$", re.M)
    matches = list(pattern.finditer(text))
    for idx, m in enumerate(matches):
        slide_no = int(m.group(1))
        start = m.end()
        end = matches[idx + 1].start() if idx + 1 < len(matches) else len(text)
        body = _clean_tts_text(text[start:end])
        if body:
            notes[slide_no] = body
    return notes


def _clean_tts_text(text: str) -> str:
    """Convert Markdown-ish script text to plain TTS-friendly prose."""
    lines: list[str] = []
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        if line == "---" or line.startswith("> **References"):
            continue
        if re.match(r"^\|?\s*:?-{3,}:?", line):
            continue
        if line.startswith("|"):
            # Tables are usually redundant with spoken prose.
            continue
        if line.startswith(">"):
            continue
        line = re.sub(r"^[-*+]\s+", "", line)
        line = re.sub(r"^\d+\.\s+", "", line)
        line = re.sub(r"`([^`]+)`", r"\1", line)
        line = re.sub(r"\*\*([^*]+)\*\*", r"\1", line)
        line = re.sub(r"\*([^*]+)\*", r"\1", line)
        line = re.sub(r"\[([^\]]+)\]\([^)]+\)", r"\1", line)
        line = re.sub(r"\s+", " ", line)
        if line:
            lines.append(line)
    return " ".join(lines).strip()


async def _edge_tts_save(slide_no: int, text: str, out_path: pathlib.Path,
                         voice: str, rate: str) -> tuple[int, str, int]:
    try:
        import edge_tts  # type: ignore
    except ImportError as exc:  # pragma: no cover - dependency guard
        raise RuntimeError(
            "edge-tts is required. Install with: pip install edge-tts"
        ) from exc
    out_path.parent.mkdir(parents=True, exist_ok=True)
    await edge_tts.Communicate(text, voice, rate=rate).save(str(out_path))
    return slide_no, str(out_path), out_path.stat().st_size


def _estimate_duration_seconds(text: str, wpm: int = 145) -> float:
    words = re.findall(r"[A-Za-z0-9]+|[一-鿿ぁ-んァ-ヴー]+", text)
    return max(1.5, len(words) / max(1, wpm) * 60.0)


def _media_duration_seconds(shape: Any, fallback_text: str) -> float:
    try:
        length_ms = float(shape.MediaFormat.Length)
        if length_ms > 0:
            return length_ms / 1000.0
    except Exception:
        pass
    return _estimate_duration_seconds(fallback_text)


def _set_audio_trigger(slide: Any, shape: Any, trigger_type: int) -> int:
    set_count = 0
    try:
        seq = slide.TimeLine.MainSequence
        for j in range(1, seq.Count + 1):
            eff = seq.Item(j)
            try:
                if eff.Shape.Id == shape.Id:
                    eff.Timing.TriggerType = trigger_type
                    set_count += 1
            except Exception:
                continue
    except Exception:
        pass
    return set_count


def _remove_existing_audio_shapes(slide: Any) -> int:
    removed = 0
    for i in range(slide.Shapes.Count, 0, -1):
        shape = slide.Shapes.Item(i)
        try:
            media_type = int(shape.MediaType)
        except Exception:
            media_type = 0
        try:
            name = str(shape.Name)
        except Exception:
            name = ""
        if media_type != 0 or name.startswith("radia_mcp_audio_"):
            try:
                shape.Delete()
                removed += 1
            except Exception:
                pass
    return removed


def _run_blocking(coro):
    """Run a coroutine to completion whether or not an event loop is already
    running in this thread.  Inside the MCP server the tool handler is itself
    awaited, so asyncio.run() raises "cannot be called from a running event
    loop"; a worker thread with its own loop is the way through."""
    try:
        asyncio.get_running_loop()
    except RuntimeError:
        return asyncio.run(coro)
    import concurrent.futures

    with concurrent.futures.ThreadPoolExecutor(max_workers=1) as pool:
        return pool.submit(asyncio.run, coro).result()


def presentation_embed_tts_audio_in_pptx(
    pptx_path: str,
    script_md_path: str,
    output_pptx_path: str | None = None,
    audio_dir: str | None = None,
    voice: str = "en-US-AndrewMultilingualNeural",
    rate: str = "+0%",
    first_slide_on_click: bool = True,
    following_slides_with_previous: bool = True,
    auto_advance_after_audio: bool = False,
    advance_padding_seconds: float = 0.35,
    overwrite_output: bool = True,
    remove_existing_audio: bool = True,
    update_speaker_notes: bool = True,
) -> dict:
    """Generate per-slide TTS MP3 and embed it into a PowerPoint deck.

    Args:
        pptx_path: Source PPTX path.
        script_md_path: Markdown script. Preferred section format is
            ``## Slide N — title``. TTS style ``Slide N. title.`` is also
            accepted.
        output_pptx_path: Destination PPTX. Defaults to ``*_audio.pptx``.
        audio_dir: Directory for generated MP3 files. Defaults to
            ``<pptx_dir>/audio``.
        voice: edge-tts voice name.
        rate: edge-tts rate string, e.g. ``+0%`` or ``-5%``.
        first_slide_on_click: Trigger slide 1 audio on click.
        following_slides_with_previous: Trigger slides 2..N with previous.
        auto_advance_after_audio: Set slide transition to advance after the
            embedded audio duration plus padding. Uses PowerPoint media length
            when available, otherwise estimates from script length.
        advance_padding_seconds: Extra seconds added to auto-advance time.
        overwrite_output: Allow replacing an existing destination file.
        remove_existing_audio: Delete existing audio/media shapes before adding.
        update_speaker_notes: Copy narration text into speaker notes.

    Returns:
        dict with output path, generated audio files, and per-slide status.
    """
    src = pathlib.Path(pptx_path)
    script = pathlib.Path(script_md_path)
    if not src.exists():
        return {"error": f"pptx not found: {pptx_path}"}
    if not script.exists():
        return {"error": f"script markdown not found: {script_md_path}"}

    out = pathlib.Path(output_pptx_path) if output_pptx_path else (
        src.with_name(src.stem + "_audio.pptx")
    )
    audio_root = pathlib.Path(audio_dir) if audio_dir else src.parent / "audio"
    notes = _parse_slide_markdown(script)
    if not notes:
        return {"error": f"no slide narration blocks found in {script_md_path}"}
    if out.exists() and not overwrite_output:
        return {"error": f"output exists: {out}"}

    async def gen_all() -> list[tuple[int, str, int]]:
        tasks = []
        for slide_no, text in sorted(notes.items()):
            audio_path = audio_root / f"slide_{slide_no:02d}_en.mp3"
            tasks.append(_edge_tts_save(slide_no, text, audio_path, voice, rate))
        return await asyncio.gather(*tasks)

    audio_rows = _run_blocking(gen_all())
    audio = {n: pathlib.Path(path) for n, path, _size in audio_rows}

    try:
        import win32com.client  # type: ignore
    except ImportError as exc:
        return {
            "error": "pywin32 is required for PowerPoint COM automation.",
            "detail": str(exc),
        }

    ppt = win32com.client.Dispatch("PowerPoint.Application")
    ppt.Visible = True
    work = pathlib.Path(tempfile.gettempdir()) / f"radia_mcp_audio_{os.getpid()}.pptx"
    shutil.copyfile(src, work)
    prs = ppt.Presentations.Open(str(work), WithWindow=False)

    # PowerPoint trigger constants:
    # msoAnimTriggerOnPageClick = 1, msoAnimTriggerWithPrevious = 2
    per_slide: list[dict] = []
    try:
        n_slides = int(prs.Slides.Count)
        for slide_no in range(1, n_slides + 1):
            slide = prs.Slides(slide_no)
            row: dict[str, Any] = {"slide_no": slide_no, "has_audio": False}
            if remove_existing_audio:
                row["removed_audio_shapes"] = _remove_existing_audio_shapes(slide)

            if slide_no in notes and update_speaker_notes:
                try:
                    slide.NotesPage.Shapes(2).TextFrame.TextRange.Text = notes[slide_no]
                    row["notes_updated"] = True
                except Exception as exc:
                    row["notes_updated"] = False
                    row["notes_error"] = str(exc)

            if slide_no in audio:
                shape = slide.Shapes.AddMediaObject2(
                    str(audio[slide_no].resolve()), False, True, -100, -100, 0, 0
                )
                try:
                    shape.Name = f"radia_mcp_audio_{slide_no:02d}"
                except Exception:
                    pass
                ps = shape.AnimationSettings.PlaySettings
                ps.PlayOnEntry = True
                ps.HideWhileNotPlaying = True
                ps.StopAfterSlides = 1

                if slide_no == 1 and first_slide_on_click:
                    trigger = 1
                    trigger_name = "on_click"
                elif slide_no >= 2 and following_slides_with_previous:
                    trigger = 2
                    trigger_name = "with_previous"
                else:
                    trigger = 1
                    trigger_name = "on_click"
                row["trigger_effects_set"] = _set_audio_trigger(slide, shape, trigger)
                row["trigger"] = trigger_name
                row["has_audio"] = True
                row["audio_file"] = str(audio[slide_no])

                if auto_advance_after_audio:
                    duration = _media_duration_seconds(shape, notes[slide_no])
                    advance = max(1.0, duration + float(advance_padding_seconds))
                    try:
                        slide.SlideShowTransition.AdvanceOnTime = True
                        slide.SlideShowTransition.AdvanceTime = advance
                        row["advance_seconds"] = round(advance, 2)
                    except Exception as exc:
                        row["advance_error"] = str(exc)
            per_slide.append(row)

        if out.exists() and overwrite_output:
            out.unlink()
        prs.SaveAs(str(out.resolve()))
    finally:
        try:
            prs.Close()
        except Exception:
            pass
        try:
            work.unlink()
        except OSError:
            pass
        # Do not quit PowerPoint; the user may have other decks open in the
        # same application instance.

    return {
        "output_pptx_path": str(out),
        "source_pptx_path": str(src),
        "script_md_path": str(script),
        "audio_dir": str(audio_root),
        "n_script_blocks": len(notes),
        "n_audio_files": len(audio),
        "voice": voice,
        "rate": rate,
        "auto_advance_after_audio": auto_advance_after_audio,
        "per_slide": per_slide,
    }
