---
name: pdf2ppt-pdfgear
description: "Convert a PDF into a PowerPoint .pptx automatically using the locally-installed PDFgear desktop app on Windows -- no manual clicking. Use this whenever the user wants to turn a PDF into PPT/PowerPoint slides VIA PDFGEAR SPECIFICALLY -- e.g. 'use PDFgear to convert this report.pdf to pptx', 'automate the PDFgear PDF->PPT conversion', 'batch-convert these PDFs to PowerPoint with PDFgear', 'PDFgearでPDFをパワポに変換して' -- or when an earlier part of the conversation already settled on PDFgear as the converter. It drives PDFgear's pdfconverter.exe in PDFToPPT mode through Windows UIAutomation (pywinauto): writes PDFgear's conv_docs.json, launches the converter pre-loaded (no editor), fires the 変換/Convert button via the InvokePattern (no mouse), waits for the cloud result, and moves the .pptx to the requested path. Do NOT use it for: authoring or editing a PowerPoint deck (use the pptx skill), building slides from scratch, converting PPTX->PDF (the reverse), or any offline/confidential conversion -- with 上級モード (Advanced Mode, the skill's default) PDFgear emits an EDITABLE pptx (text boxes + vector shapes), or with --no-advanced one full-page IMAGE per slide; either way it UPLOADS the PDF to PDFgear's cloud (apiw.pdfgear.com), so for an OFFLINE/private result prefer a local PyMuPDF->python-pptx render. Windows-only; requires PDFgear installed plus `pip install pywinauto`, and PDFgear able to convert without an interactive sign-in."
---

# pdf2ppt-pdfgear — Automated PDF → PPTX via PDFgear

Fully automated PDF→PowerPoint using the locally installed **PDFgear**
desktop app. No GUI clicking by hand — the skill drives PDFgear through
UIAutomation.

## Usage

```powershell
python "C:\Users\Administrator\.claude\skills\pdf2ppt-pdfgear\pdf2ppt_pdfgear.py" "input.pdf"
python ".../pdf2ppt_pdfgear.py" "input.pdf" -o "S:\path\out.pptx"
python ".../pdf2ppt_pdfgear.py" "input.pdf" --timeout 240
```

- No `-o`: result is left at `<PDFgear 出力先>\<stem> conv.pptx`
  (default `%USERPROFILE%\OneDrive\PDFgear\<stem> conv.pptx`).
- With `-o`: the result is moved to that path.
- Prints the final `.pptx` path on success; exits non-zero with a clear
  message on failure (timeout / sign-in required / UI changed).

## How it works (reverse-engineered 2026-06-02)

PDFgear's standalone converter accepts positional args and opens
pre-loaded in PDF→PPT mode (this is exactly what the editor's
"PDFからPPTへ" button does under the hood):

```
pdfconverter.exe "app1" "PDFToPPT" "<conv_docs.json>"
```

`conv_docs.json` schema:

```json
{"ConvertType": "PDFToPPT",
 "FilesPath": [{"FilePath": "C:\\abs\\in.pdf", "OpenDate": "<.NET ticks>", "Password": null}]}
```

The skill: writes that JSON → launches the converter → finds its window by
**pdfconverter.exe PID** (the title carries the document name, so don't match
on title) → ensures the **`上級モード` (Advanced Mode) CheckBox is ON** via the
UIA TogglePattern (this is what makes the output EDITABLE; default on) → finds
the `変換` button → fires it via the UIA **InvokePattern** (works even though
the button reports `visible=False`; no mouse movement) → polls the output dir
for `<stem> conv.pptx` → moves it to `--output`.

Key automation facts:
- The local `pdfconverter.exe` menu (結合/分割/圧縮/X→PDF) does **NOT** list
  PDF→PPT; that direction only exists via the `app1 PDFToPPT` launch (cloud).
- `MouseOverPdfToPPTButton` is the editor-side entry button (aid), if you ever
  need the editor route instead of the direct JSON launch.
- Output dir is a PDFgear *setting* (出力先), not in the JSON; default is
  `OneDrive\PDFgear`. Override with `--outdir` if the user changed it.
- PrintWindow with `PW_RENDERFULLCONTENT (2)` screenshots PDFgear reliably
  even when other windows overlap (BitBlt screen-grab does not).

## Caveats — read before recommending this skill

1. **Output mode = 上級モード (Advanced Mode) checkbox.** The skill checks it
   by default (`--no-advanced` to skip):
   - **上級モード ON → EDITABLE** pptx: text boxes + vector shapes, **0 raster**
     — titles / body / math / figures all editable, high visual fidelity
     (verified on the 14-slide CEFC beamer deck: 315 text boxes, 0 pictures).
   - **上級モード OFF → image-per-slide** (one `PICTURE` per slide), the same
     result as a local `PyMuPDF`→`python-pptx` render.
   (The earlier "always image-only, not editable" claim was WRONG — it was
   measured with 上級モード OFF, which the skill did not used to set.)  Even in
   editable mode, figures become many small vector shapes and beamer math is
   re-laid-out per glyph, so fidelity is high but not byte-identical.
2. **Cloud upload.** The conversion POSTs the PDF to `apiw.pdfgear.com`.
   Do not use for confidential or unpublished documents.
3. **Sign-in.** On a machine where PDFgear's free cloud convert needs an
   interactive login, the unattended click will stall → the skill times out
   with a clear message. Pre-sign-in in the PDFgear GUI once.
4. **Live desktop.** Uses InvokePattern (no mouse hijack), but the converter
   window does briefly appear; PDFgear also pops an Explorer window at the
   output dir (its "変換後エクスプローラーで表示" setting).

## When to use something else

- **Offline / identical image quality / no cloud** → local render:
  `fitz` (PyMuPDF) page → PNG → `python-pptx` full-bleed picture on a 16:9
  slide. Deterministic, instant, private. (This is what produced
  `cefc2026_oral.pptx`.)
- **Genuinely editable PowerPoint** → 上級モード ON (this skill's default) gives
  editable text boxes + vector figures — good enough to tweak text / move
  objects — but beamer **math is re-laid-out per glyph** (messy to re-edit) and
  figures become many tiny shapes. For clean editable math, the Beamer `.tex`
  is the real editable master (edit + recompile); re-author in PowerPoint only
  if you must leave LaTeX.

## Dependencies

`pip install pywinauto` (uses `comtypes`). Windows + PDFgear installed
(`C:\Program Files\PDFgear\pdfconverter.exe`).
