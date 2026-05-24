# EasyOCR Queue — Lab Library Books Awaiting OCR

Generated from `catalog.json` (PyMuPDF batch index) on 2026-05-21.

The PDF batch indexer flagged 62 entries as `is_scanned: true`.  Of
these, the **14 priority books below** have no usable text layer
(`abstract_excerpt < 100 chars` despite `page_count >= 100`) and
should be passed through EasyOCR before they can contribute to
`bibliography_index_knowledge.py`.

The remaining 48 entries are mostly **single-page figure PDFs**
(eps-converted-to outputs from NMR/MR papers) which do not need
OCR.

---

## Priority 1 — Major Japanese 電磁気学 textbooks (0-60 chars OCR)

| Pages | Title | Path |
|-------|-------|------|
|  473  | 電磁気学_紀伊国屋書店_砂川重信 | `W:\03_文献・論文\00_電磁界解析\01_教科書\10_電磁気学・物理基礎\洋書・専門\` |
|  329  | 電磁気学_岩波書店_砂川重信 | (same dir) |
|  288  | 電磁気学_コロナ社_(稲垣) | `…\01_教科書\10_電磁気学・物理基礎\和書\` |
|  242  | 電磁気学_初めて学ぶ人のために | (same dir) |
|  226  | 電磁気学_コロナ社_(石井) | (same dir) |
|  197  | 工科の物理3_電磁気学 | (same dir) |
|  139  | 電磁気学とは何か | (same dir) |

After OCR these classify into `radia_ngsolve` (general EM) and may
also extend `differential_forms` if the Sunakawa books include
exterior-calculus treatments (likely; check after OCR).

---

## Priority 2 — Specialty references (high page count, no OCR)

| Pages | Title | Subpackage hint |
|-------|-------|-----------------|
|  821  | Analysis of Multiconductor Transmission Lines | `peec` |
|  685  | NMR Spectroscopy Explained | `nmr_mri` |
|  268  | アマチュアのV・UHF技術 | `_misc` → re-classify after OCR |
|  254  | 原子力材料 | (probably out-of-scope; skip) |
|  232  | Magnetic_Hysteresis | `electromagnet` |
|  223  | Systems With Hysteresis (Ikhouane Faycal) | `electromagnet` |
|  193  | 経済学のための最適化理論入門 | `topology_optimization` / `pinn` |

`Analysis of Multiconductor Transmission Lines` (Paul) is *the* PEEC
reference and will dramatically improve `radia_mcp.peec` knowledge.

`Systems With Hysteresis` (Ikhouane-Rodellar) is the Bouc-Wen model
canonical reference and complements the existing Play / Energy model
material in `radia_mcp.differential_forms.forces_knowledge`.

---

## EasyOCR command template

```powershell
# Example: OCR one priority book and save text alongside PDF
$src = "W:\03_文献・論文\00_電磁界解析\01_教科書\10_電磁気学・物理基礎\洋書・専門\電磁気学_紀伊国屋書店_砂川重信.pdf"

python -c "
import easyocr, pdf2image, sys, json
reader = easyocr.Reader(['ja','en'], gpu=True)
pages = pdf2image.convert_from_path(sys.argv[1], dpi=200)
out = []
for i, img in enumerate(pages):
    lines = reader.readtext(img, detail=0, paragraph=True)
    out.append({'page': i+1, 'text': '\n'.join(lines)})
    print(f'page {i+1}/{len(pages)} done')
import json
open(sys.argv[1] + '.ocr.json', 'w', encoding='utf-8').write(json.dumps(out, ensure_ascii=False, indent=2))
" "$src"
```

(GPU EasyOCR runs ~2-3 s/page; the 685-page NMR book takes ~30 min.)

After OCR, re-run the batch indexer with EasyOCR output integrated
to repopulate `abstract_excerpt`, then re-run
`catalog_to_knowledge.py` to refresh the per-subpackage indexes.

---

## Lower-priority small scans (1-50 pages)

Mostly NMR/MRI figure files (1-page PDF per figure) — these are
decorative outputs, not text-bearing.  Skip.

Two non-figure low-priority entries:

| Pages | Title | Subpackage |
|-------|-------|------------|
|   56  | 信号処理_fpga-tuner_seminar | `radia_ngsolve` (slide deck) |
|   33  | 電磁場解析解析で役に立つ解析積分公式集 | `radia_ngsolve` (analytical formulas — relevant!) |
|   33  | PDF_2022_11_14_JTC_配布資料② | `_misc` |
|   36  | 宇宙科学研究所報告177_200102 | `_misc` |
|   17  | 基礎電磁気学演習(第15回) | `_textbooks` |
|    4  | 反磁場係数 | `_misc` (likely demagnetization-factor table — short but useful) |
