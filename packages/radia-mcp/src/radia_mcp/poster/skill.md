# poster — 学会ポスター発表の作文・組版支援

学会ポスター (A1/A0、縦・横レイアウト) 特化の lint + テンプレ + コンパイル。
スライド向けの [`presentation`](../presentation/skill.md) とは別軸 — ポスター
は「1m 離れて読む」ことが固有の制約なので、フォントサイズ・図優位・密度
の基準が大きく異なる。

## 起源と進化

- v1: COMSOL Conference 2025 Tokyo (2025-12-05) で発表した `Kelvin.tex` を雛形
  promote (5 tools: template / lint / compile / figures_audit)
- v2 (本版): web 調査 (Mike Morrison #betterposter, Colin Purrington,
  Editage, BetterPosters blog 2025, WCAG 2.1, ScienceUX eye-tracking
  pilot) を反映して 31 tools / 7 Tier に拡張

## ポスター lint 基準 (presentation と異なる軸)

| 軸 | スライド (presentation) | ポスター (poster) |
|---|---|---|
| 視聴距離 | 3-5m | 1m |
| 最小フォント | 18pt | 24pt 本文 / 20pt 参考文献 (Editage) |
| 1 枚あたりの主張 | 1 主張 1 slide | 1 ポスターで 3-5 章 + 単一 finding |
| 図優位 | 図 ≥ 本文 | 60% 図 + caption / 40% 本文 (Purrington) |
| 動線 | 時系列 (slide N → N+1) | Z 字 / F 字 (左上 → 右下) |
| 参考文献 | 末尾 1 slide | 各 contentbox 末尾に分散配置 |
| 語数 | (slide ごと) | ≤1000 語 (Purrington), 5-10 分で完読 |

## Tools (31 個 / 7 Tier)

### Tier 0 — Template + Compile (v1 から継続)

| Tool | 用途 |
|---|---|
| `poster_template_kelvin(out_path=None, replace=None)` | A1 縦・bxjsarticle ポスター雛形 (`Kelvin.tex` ベース) |
| `poster_lint(tex_path)` | レガシー 5 軸 lint (新規 `poster_health_report` を推奨) |
| `poster_compile(tex_path, engine="platex")` | platex + dvipdfmx で PDF 化、PDF lock 自動解放 |
| `poster_figures_audit(tex_path)` | `\\includegraphics` 参照と figures/ 実体の整合 |
| `poster_skill_doc()` | この skill.md を取得 |

### Tier 1 — Lint 拡充

| Tool | チェック内容 | 出典 |
|---|---|---|
| `poster_word_count` | 本文 ≤1000 語 (EN words + JP chars × 0.5) | Purrington |
| `poster_line_length` | 文長 ≤35 EN words / ≤80 JP chars | Purrington + Editage |
| `poster_caption_self_contained` | 図 caption に「軸+単位+凡例+(統計)」 | Editage |
| `poster_jp_font_check` | `\\sfdefault` / `\\gtdefault` 必須、Mincho 警告 | Editage |
| `poster_fontsize_by_distance` | Editage 公式 `h_mm = d_m/0.25` で逆算 | Editage |
| `poster_typography_lints` | no double-space / Title Case / underline / ALL CAPS 等 9 規則 | Purrington |

### Tier 2 — アクセシビリティ

| Tool | チェック内容 | 出典 |
|---|---|---|
| `poster_color_contrast_wcag` | text/bg ペアで WCAG AA 3:1 / 4.5:1 判定 | WCAG 2.1 |
| `poster_colorblind_hint` | deuteranopia 模擬で衝突ペア検出 | ACS Chem Health |
| `poster_color_count_321` | `\\definecolor` ≤3 色 (70/25/5 ルール) | Editage |

### Tier 3 — BetterPoster (Morrison 派)

| Tool | 用途 |
|---|---|
| `poster_template_betterposter(out_path=None, replace=None)` | A0 横 / 25-50-25 ゾーン雛形 |
| `poster_betterposter_billboard_lint(tex_path)` | 中央 finding ≥120pt / ≤15 語 / ≤2 jargon 検証 |
| `poster_zone_balance_check(tex_path)` | minipage 幅が 0.25 / 0.50 / 0.25 ± 0.05 |

### Tier 4 — QR コード

| Tool | 用途 |
|---|---|
| `poster_qr_inject(tex_path, url, label, size_cm)` | `qrcode` package + tikz overlay で右下に挿入 |
| `poster_qr_audit(tex_path, check_url=False)` | サイズ ≥2.5cm / 具体ラベル / 任意 URL HEAD |

### Tier 5 — 印刷準備

| Tool | 用途 |
|---|---|
| `poster_print_readiness_audit(tex_path)` | 用紙宣言 + 各 figure の DPI (vector PDF は OK) |
| `poster_font_embed_check(pdf_path)` | `pdffonts` で全フォント埋め込み確認 |

### Tier 6 — Intelligence Layer

| Tool | 用途 |
|---|---|
| `poster_health_report(tex_path)` | Tier 1-2 集約 + 0-100 スコア |
| `poster_root_cause_diagnosis(tex_path)` | 5 patterns (title-heavy / figure-starved / color-chaotic / billboard-buried / typography-noisy) |
| `poster_next_5_actions(tex_path)` | impact × effort で top-5 + 修正例 |
| `poster_rewrite_suggest(target, slot_values)` | title / subtitle / billboard / abstract / takeaway / qr_label の 3-4 案 |
| `poster_adaptive_health_report(tex_path, conference)` | COMSOL / IEEJ / IEEE_Magnetics / Compumag / 応物 で severity 調整 |
| `poster_run_full_workflow(tex_path, conference, rewrite_target)` | 1 コール chain (intake → adaptive → root cause → top-5 → rewrite) |

### Tier 7 — 入口の広さ

| Tool | 用途 |
|---|---|
| `poster_from_pptx(pptx_path, out_tex_path)` | PowerPoint → A1 縦 .tex 骨格 (python-pptx) |
| `poster_from_paper_tex(paper_tex, out_tex_path, title)` | 論文 .tex → poster .tex (Intro/Method/Result/Conclusion 抽出) |
| `poster_elevator_pitch_generate(tex_path)` | 3 分発表台本 (30s hook + 90s main + 60s implication) |
| `poster_qa_anticipation_list(tex_path)` | 想定 Q&A 10 項目 (reviewer-type 分類付き) |

## クイック起動例

```python
# 1-call で health report + 修正計画 + リライト候補をまとめて取得
poster_run_full_workflow(
    tex_path="public-safe curated corpus",
    conference="IEEJ",
    rewrite_target="title",
)
```

```python
# QR コードを挿入 + audit で BetterPosters ルールに照らす
poster_qr_inject(tex_path, url="https://...", label="Scan for full paper", size_cm=4.0, in_place=True)
poster_qr_audit(tex_path)
```

```python
# 既存の論文 .tex から poster .tex を生成
poster_from_paper_tex("paper.tex", out_tex_path="poster_draft.tex")
```

## 出典 (web 調査 2026-05)

- Mike Morrison, "How to design a better research poster" (2019); OSF
  `#betterposter` template — https://osf.io/ef53g/
- ScienceUX, "Eye Tracking Pilot: #betterposter Design v2" (2020) —
  faster comprehension than dense-traditional layout
- BetterPosters Blog, "QR codes on conference posters: Some scan, some
  don't" (2025-05) — ~50% scan rate, prominence + concrete label rules
- Colin Purrington, "Designing conference posters" (2025) — ≤1000 語、
  45-65 char/line、no double-space、no Title Case 等 9 規則
- Editage『学会ポスター発表完全ガイド6』(2025) — 視認距離公式、
  Gothic 必須、3 色 70/25/5 ルール
- WCAG 2.1 (W3C) — text contrast 3:1 (AA-Large) / 4.5:1 (AA)
- ACS Chemical Health & Safety — "Beyond the Visual: Achieving
  Accessibility in Scientific Figures" (赤緑回避)
- K. Sugahara, "Electromagnetic Analysis of Eddy Current Testing With
  Kelvin Transformation," IEEE Trans. Magn., vol.58, no.9, 2022.
- K. Sugahara, "Extended Kelvin Transformation for Solving Radiating
  Electromagnetic Fields," IEICE Trans. Electron., vol.E108-C, no.4,
  pp.189-194, 2025.
