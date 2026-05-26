# Paper Writing Skill — Journal 論文を通すための作文技術

Journal paper 執筆ガイド。IEEE Transactions / IEEJ 論文誌 / APS / Elsevier 系を対象に、IMRAD 構造・contribution 明示・reviewer 対応の実戦ノウハウを統合。

---

## 🛠️ 診断ツール (Virtual Tools → Actual Tool Mapping)

v0.12.0 以降、Plan B Tier 1 + Tier 2 の実測ツールで多くがカバーされた。
各 virtual tool が「どの実装に接続されているか」を明示。

### T1. Contribution 明確度 → ✅ `paper_writing_contribution_clarity_score` (v0.12.0)
Introduction 末尾の Contribution ブロックを検出し、claim verbs (propose /
present / demonstrate 等) 数、item count (3-5 理想)、parallel POS、hedge
不在、定量要素の有無を総合スコアリング。
- **NG**: "In this paper, we discuss..." (discuss だけで貢献が不明)

### T2. IMRAD バランス → ✅ `paper_writing_check_imrad_balance`
IMRAD 構造の字数バランス。
- 目安: Intro 15%, Method 25%, Result 30%, Discussion 20%, Concl 5%, Abs 5%
- Discussion < Method は要警戒 (結果を並べて考察が浅い)

### T3. Reviewer-2 triggers → 🔗 複合 (既存 + 新規)
査読者が確実に突っ込む語彙を検出。以下のツールを組み合わせ:
- 弱気修飾語: `paper_writing_count_weak_expressions`
- 定量化されていない主張: **`paper_writing_claim_quantification` (v0.13.0)**
- 未定義略語の初出: `paper_writing_find_undefined_acronyms` (v0.13.0 cross_lint)
- English red flags: `paper_writing_check_english_redflags`

### T4. Figure caption 品質 → ✅ `paper_writing_check_figure_caption_showing`
Figure caption 単独で図が理解できるか。
- 目標: caption だけ読めば図の主張が伝わる
- 推奨: "Fig. 3. <主張の動詞句>: <セットアップ>, <キー数値>"

### T5. Citation density → ✅ `paper_writing_related_work_density` (v0.13.0)
Introduction 内の引用密度・自己引用比率・年度分布を診断。
- 目安: Intro に 20-40 件 (研究動向 20, 直接比較対象 5-10)
- 自己引用比率 < 20%、直近 5 年分 40% 以上が目安

### T6. 制約検証 → ✅ `paper_writing_validate_pdf_pages` / `paper_writing_validate_abstract_length` / `paper_writing_check_overfull_hbox` / `paper_writing_check_pdf_edge_overflow`
IEEE/IEEJ/APS 固有の制約検証:
- ページ数 (IEEE Trans: extended OK / PRL 4 pages / IEEJ 和文 10 pages)
- Abstract (IEEE 200 words / IEEJ 400 字 / APS 250 words)
- Overfull hbox ゼロ (LaTeX log 経由 = `check_overfull_hbox`)
- **紙面端へのはみ出しゼロ (PDF を直接検査 = `check_pdf_edge_overflow`)**
  log には出ないケース (\\sloppy + underfull、長い URL、長い数式、figure*
  忘れ、table 列数オーバー) も検出。pymupdf で各 span の bbox を取得し、
  全 span 右端の 99 パーセンタイルを「通常の右端」として、それを 1pt 以上
  超える span を flag。`page_boundary_violation = True` のものは MediaBox
  の外に出ているので投稿前に必ず修正。
- **明らかな体裁エラー 7 カテゴリ一括 (`check_pdf_obvious_errors`)**:
  unresolved_references (`??`/`[?]`), bibkey_residues (`[Foo2024]`),
  raw_latex_residues (`\\frac` 等), **pseudo_math_residues (ASCII pseudo-math
  by IEEEtran bibstyle render: `(R/rho)^2 modulated nu` → `(R/rho)2modulatednu`)**,
  blank_pages, font_issues (Type 3 reject), image_bbox_violations.
  Score 0-10 (10 = 全カテゴリクリーン)。投稿前最終チェック必須。
- **Reviewer-2 高頻度異常 8 カテゴリ (`check_pdf_advanced_anomalies`)**:
  obvious_errors 後の更に深い検査。citation_order (IEEE 初出順 [1][2][3]) /
  citation_space (`word[N]` NG → `word [N]`) / acronym_first_use (二重定義・
  未定義使用) / caption_position (表は上・図は下) / equation_gaps (式番号スキップ・
  重複・未参照) / heading_capitalization (Title Case 揃い) / orphan_widow
  (1-2 行取り残し段落) / image_overlap (図 vs 図、図 vs body)。Score 0-10。

### T7. Cover letter 診断 → 🟡 部分的 `paper_writing_generate_cover_letter` (生成のみ、診断はまだ mental)
Cover letter が editor にアピールできているか。
- **必須**: why this journal, what's new, why important NOW
- **NG**: "Please consider our paper" の一文だけ

### T8. Abstract 再構成 → ✅ `paper_writing_abstract_strength` (v0.12.0)
4 要素 (problem / method / quantitative-result / impact) + strong-verb +
banned-phrase + background-ratio の total score。
1. 問題 (1-2 文): 既存研究の limitation
2. 手法 (1-2 文): 本研究の approach
3. 結果 (1-2 文): 定量的 key number
4. 意義 (1 文): implication / contribution

### 新規 (v0.13.0 Plan B Tier 2) — skill.md に無かった追加ツール

- `paper_writing_limitation_statement_presence` — Discussion 内の
  limitation 段落の存在・位置 (mid 推奨)・充実度 (skill.md §O / Wallwork §9.12)
- `paper_writing_figure_referencing_coverage` — 全 \\label が本文で
  \\ref される回数を集計、unreferenced / singleton を flag (Tufte)
- `paper_writing_given_new_ordering` — Wallwork §3.4-3.6 Given-New 配置
  の段落内 heuristic 評価

### 和文 lint (v0.13.0 cross_lint re-exports)

grant_writing 由来の 8 本を paper_writing 名前空間で露出:
notation_variants / find_undefined_acronyms / acronym_usage_audit /
check_kanji_ratio / lint_bedrock / check_misuse_japanese /
suggest_redundancy_fixes / check_subject_predicate_distance

### Plan B Tier 3 (v0.19.0) — paper 特有の追加診断 6 ツール

grant_writing にはない paper 特有の観点 (truth-claim 軸 / IMRAD 後処理 /
journal fit) で追加された 6 ツール:

#### `paper_writing_title_abstract_conclusion_triangle(text)` — 三位一体整合性
論文の **入口 (Title) と要約 (Abstract) と出口 (Conclusion)** が同じ物語を
語っているかを **三角形 keyword overlap** + 数値 claim 整合性で診断。
Title 過大広告 / Conclusion-Abstract 段差は reviewer の信頼を一発で失う。
Harmonic mean を score 化するので、3 辺のうち最も弱い辺が全体を引っ張る。
出典: Day & Gastel (7th ed) Ch.5 + Wallwork Ch.12。

#### `paper_writing_reproducibility_open_science_check(text, journal_required)` — 6 軸
近年 IEEE / Nature / Science / PLoS が必須化する Open Science 6 軸:
code_availability / data_availability / preregistration /
methods_replicability / random_seed_disclosure / competing_interests。
`journal_required="code,data,coi"` のようにカンマ区切りで journal 必須項目を
渡すと、不在時に critical 扱いで減点。
出典: TOP Guidelines (Nosek 2015 Science) + NeurIPS Reproducibility Checklist
+ ICMJE recommendations。

#### `paper_writing_statistical_reporting_compliance(text)` — p+effect+CI+n
**p 値単独報告は reviewer #2 トリガー筆頭**。各 p 値 occurrence の周辺
(default 250 字以内) で 4 要素 (p / effect size / CI / sample size) の
充足を per-occurrence 評価。多重比較補正 (Bonferroni / FDR) の言及不足も
flag。理論系論文 (p 値 0 件) では適用外として score=10 を返す。
出典: APA 7th + NEJM Statistical Reporting + CONSORT 2010 + STROBE。

#### `paper_writing_citation_health_4_axes(bib, author_last_names, current_year)` — 4 軸
.bib 入力で 4 軸診断: **recency** (median age >10年で warning) /
**author_concentration** (top1 share >20% で echo chamber) /
**geographic_diversity** (1 region 偏重を flag) /
**self_citation** (>30% で manipulation 疑い)。
既存 `check_self_citation_ratio` と相補的に細分化。
出典: Bornmann & Daniel (2008) + ICMJE + DORA。

#### `paper_writing_discussion_structure_4_elements(text)` — Discussion 4 要素
Discussion section を抽出し 4 要素を評価:
**interpretation_implications** / **limitations** / **future_work** /
**generalization_caveats**。Discussion を Results 再記述で終わらせがちな
典型 NG を防ぐ。well_developed = 同要素 hit ≥2 件。既存 T4
(limitation_statement_presence) は paper 全体だが本ツールは Discussion
section 内のみ厳格評価。
出典: Day & Gastel (7th ed) Ch.13 + Wallwork Ch.13。

#### `paper_writing_journal_fit_assessment(text, journal_name, aims_and_scope)` — Pre-submission
target journal の **Aims & Scope** を caller が貼り付けて、論文の
title / abstract / keywords との F1-style overlap を診断。score < 4 は
desk reject リスク高。Web で取得した aims_and_scope を渡す前提。
共通 keywords / 不足 critical terms / paper-only terms を提示し、
abstract 微修正の手がかりを返す。
出典: editor desk reject 統計 + journal manuscript central guidance。

#### `paper_writing_check_typography_hacks(tex_path)` — フォント縮小・余白拡張の自動検出
**Sugahara 2026-05-21 "fontを小さくするはだめ" を automate**。`.tex` ソースを
scan して page-limit hack を検出:

- **CRITICAL**: `\fontsize{X}{Y}\selectfont` で X<10、本文中の `\scriptsize` / `\tiny`
- **MAJOR**: 本文中の `\small` / `\footnotesize` (caption/table/figure 外)、
  `\linespread{<0.95}`、`baselineskip<12pt`、`\setlength{\textheight/textwidth/...}`、
  `\geometry{...}` (journal class margin の上書き)
- **INFO** (score に影響しない): `\vspace{-Xmm}` 等の負 vspace。
  Sugahara 2026-05-21 「\vspace は、見た目が良くなるならありにしよう」 —
  `\vspace` は layout tool として OK、報告のみ (audit 用)。
  「見た目改善」 vs「内容詰め込み」は context 判断。

Context-aware: caption / table / figure 環境内の `\small` 等は IEEE 標準なので
flag しない。

Score 0-10、検出ごとに critical=-3.0, major=-1.5, warning=-0.5。
CRITICAL 検出時は recommendation が「全部消せ、その上で `suggest_concept_drops`」
に切り替わる。
出典: IEEE TMag author guidelines (10pt 最低)。

#### Page-limit overflow への対処 — **タイポグラフィ ハックは厳禁**
ページ制限を超えた時の **bad responses 階層** (上ほど絶対禁止):

0. **Typography hacks — 厳禁** (Sugahara 2026-05-21):
   - 本文フォント縮小 (`\small` / `\footnotesize` / `\fontsize{9}`) — IEEE/IEEJ
     の 10pt 最低基準を破る。
   - 行間圧縮 (`\linespread{0.9}` 以下、`baselineskip` 12pt 以下)。
   - マージン / textheight 拡大 (journal class default より広く取る)。
   - **`\vspace` は厳禁ではない** (Sugahara 2026-05-21 refinement,
     「見た目が良くなるならあり」): layout tool として使用 OK、detector は
     INFO 報告のみで score に影響させない。

   font/spacing/margin の3つは desk-reject トリガー。reviewer は
   journal class の標準と違うタイポグラフィを一瞥で見抜く — 「正規
   スペースで収まらなかった」ことを露呈する。`\vspace` は通常 LaTeX
   layout の一部 — 全体のタイポが標準に見えるかどうかが基準。
1. **Prose 圧縮過剰** — 下記 `check_prose_density` で検出される anti-pattern
   (nominalisation, em-dash chain, jargon クラスタ).
2. **概念を1つ削る** — 正解。`suggest_concept_drops` で候補を提案。

#### `paper_writing_check_prose_density(text)` — 圧縮 anti-pattern 検出
語数・ページ制限に達して **悪い圧縮** を始めた draft を検出する per-sentence
診断。5 軸: nominalisation (`we augment` → `the augmentation`) /
em-dash + semicolon clause chaining / jargon クラスタ / 40 語超の長文 /
acronym 密集。各文 score 0-4、>=2 で flag。
**flag 率 >30% の場合**、recommendation が `rewrite` ではなく
**`reduce_content`** に切り替わる — 圧縮には床がある。これ以上圧縮しても
読めなくなるだけで、本来は **概念を 1 つ落とすべき** という政策判断。
全体 score 0-10。表面 detector (font / citation / hbox) が全部通っても
読みにくい draft の最後の砦。

#### `paper_writing_suggest_concept_drops(text)` — どの概念を落とすか
`check_prose_density` の companion。compression floor を超えた時、
**具体的に何を捨てれば自然な文章に戻るか** を pattern match で提案する。
5 droppable パターン (低コスト順):

1. **CITATION_ONLY_PARENTHETICAL** (cost 1): `(see [N])` / `(cf.\cite{X})`
   — 本文に同じ引用があれば削除可。
2. **PARALLEL_CITATION** (cost 2): `X \cite{A} and Y \cite{B}` で
   両方が同じ修辞役割。2 つ目を削除。IGTE digest で manual に発見した
   パターン (DtN \cite{GivoliKeller1989} を Warburg \cite{Randles1947}
   の隣で drop) — 11 語節約 + prose_density 6.5→10.0。
3. **EM_DASH_INTERPOLATION** (cost 2): `A --- X --- B` の X が
   3-15 語。査読者要望の遅い追加で sentence flow を切る典型。
4. **PARENTHETICAL_ASIDE** (cost 3): 非引用カッコ補足。
5. **TRAILING_SEMICOLON_CLAUSE** (cost 3): `main; X also Y`。

出力: ranked candidate list (drop_text / words_saved / rewrite_preview /
cost / rationale)。recommendation field に「#1 を drop すれば 1-page
制約を満たすかも」のヒント。
**注意**: 候補は提案であり機械的に削るな — 投稿先 (digest / journal /
extended abstract) ごとに残す/削るは判断が違う。
出典: IGTE 2026 digest 2026-05-21 debugging trace。

---

## 📋 Workflow Phases

### Phase 1: 📖 **intake** — ドラフト全体を把握
- Abstract と Conclusion を先に読む
- 図表を章立てなしで並べて物語になっているかチェック
- IMRAD の各 section がどこから始まるか確認

### Phase 2: 🔍 **diagnose** — 全体診断
- `paper_writing_analyze_sentences(text)` で平均文長
- `paper_writing_count_weak_expressions(text)` で弱気修飾語
- `paper_writing_check_citation_usage(tex, bib)` で引用整合
- `paper_writing_check_english_redflags(text)` で英文 red flag
- Contribution が 1 段落内に明示されているか目視

### Phase 3: ✍️ **rewrite** — 重要 section 優先で書き直し
優先順位:
1. **Abstract + Introduction 末尾** (rejection の過半はここで決まる)
2. **Figure captions** (ページをめくる reviewer に最初に読まれる)
3. **Contribution 宣言 (intro 末尾)**
4. **Discussion** (method との対比で主張)
5. Method (最後。最も書きやすいので後回し)

### Phase 4: ✅ **validate** — 制約検証
- `paper_writing_validate_pdf_pages(pdf, limit)` でページ数
- `paper_writing_validate_abstract_length(text)` で abstract 字数
- `paper_writing_check_overfull_hbox(log)` で組版警告
- 査読用 checklist (journal 公式) を順に潰す

### Phase 5: 📝 **submit** — 投稿
- Cover letter (1 ページ)
- Supplementary material 整理
- reference style journal 規定にフォーマット
- ORCID / 著者 affiliation 確認

### Phase 6: 🔄 **revise** — reviewer 対応
- Response letter は query を明示、line-by-line で対応
- 変更箇所を colored text で本文にも示す
- 反論するなら **データで**。感情論は禁。

---

## 📏 定量的閾値 (Quantitative Targets)

| 項目 | 目標 | 警告 |
|---|---|---|
| Abstract 語数 (EN) | 150-200 words | >250 |
| Abstract 字数 (JP) | 300-400 字 | >500 |
| 平均文長 (JP) | ≤60 字 | ≥100 字 |
| 弱気修飾語 | 0 (Abs/Contribution) | ≥3/page |
| Overfull hbox | 0 | ≥1 |
| 下線 | 0 | ≥5/page (journal では稀) |
| 図表数 | 6-10 (原著) | >12 (PRL は 5-6 上限) |
| 引用総数 | 30-50 | >80 (review 論文以外は過剰) |
| 自己引用率 | <20% | ≥30% |
| 未使用 bib entry | 0 | ≥5 |

---

## 📚 基本理論

### A. IMRAD の各 section の機能

| Section | 役割 | キーフレーズ例 |
|---|---|---|
| Abstract | 論文全体の 30 秒要約 | "We propose ... and demonstrate ..." |
| Introduction | なぜ今これが重要か + 本論文の貢献 | "However, prior work has not addressed ..." |
| Method | 再現性を担保する記述 | "The measurement setup consists of ..." |
| Result | 生データ + 代表的な figure | "Fig. 3 shows ... at X mT" |
| Discussion | 結果の物理的解釈 + 限界 | "This behavior can be explained by ..." |
| Conclusion | 貢献の再確認 + 今後の展望 | "We have demonstrated ..." |

### B. Introduction の 4 段構造 (Swales CARS model)

1. **Territory** — この分野は重要
2. **Niche** — 既存研究にはギャップがある
3. **Occupy** — 本論文はそのギャップを埋める
4. **Announce** — 具体的な貢献 (contribution list)

### C. Contribution list の書き方

```
The main contributions of this paper are:
1. We propose <手法名> for <解く問題>.
2. We demonstrate <key result> with <quantitative metric>.
3. We show <implication> for <broader field>.
```

各項目で **動詞を明示**。"We discuss / We consider" は弱い。

### D. Reviewer 2 対策

"Reviewer 2" とは厳しい査読者の俗語。典型的 trigger:
- **Contribution が曖昧**: "What's new compared to [17]?" → Intro 末尾で明示
- **弱気修飾**: "might / could / seems" → 断定的に ("We show", "We demonstrate")
- **数値なし**: "significantly faster" → "3.2× faster (Fig. 7)"
- **fair comparison なし**: 自分の手法だけ最適化、既存は default → 条件を揃える
- **Limitation 未記述**: discussion で limit を 1 段落書く。正直に。

### E. Figure の書き方 (Tufte principles)

- データ・インク比を最大化 (背景色・凡例枠は最小)
- Caption は「図の要約」、本文を読まなくても主張が伝わる
- 単位・凡例・軸ラベルは必ず (学振審査員ですら凡例なし図で落とす)
- 白黒印刷でも区別できる (hatching / line style)
- figures/ ディレクトリに pdf/svg と source (pptx/draw.io) を並置

---

## 🗂️ OK/NG 文例

### Contribution 宣言

**NG**:
> In this paper, we discuss the effect of eddy current on the induction heating workpiece.

**OK**:
> In this paper, we **propose** a hybrid surface-integral + reduced-order solver for the target problem class. We **demonstrate** a 4× speedup over the reference finite-element pipeline with <3% accuracy loss.

### Abstract opening

**NG**:
> Induction heating is a widely used technology.  (generic truism)

**OK**:
> Resolving the 0.1 mm skin depth in a 50 kHz eddy-current problem requires a fine mesh whose cost forbids parameter sweeps over more than a handful of operating points.  (problem stated with numbers)

### Discussion の深さ

**NG**:
> Fig. 5 shows that the power loss increases with frequency.

**OK**:
> Fig. 5 shows an f^{0.5} dependence of the power loss, consistent with the classical skin-effect scaling P ∝ δ^{-1}, where δ = sqrt(2/ωμσ). The 3% deviation at high frequency arises from ...

### 英文の hedging

**NG**:
> Our method might be useful for this kind of problem.

**OK**:
> Our method reduces the simulation time from 6 h to 1.5 h (Table II), making parametric sweeps over 50 operating points feasible within a single day.

---

## 🚫 NG パターン 10 (journal 特有)

1. **Contribution が discussion 調** ("We discuss..." で終わる) — 動詞で宣言
2. **Abstract に将来形 / 提案形** ("We will show") — 現在形 ("We show")
3. **fair comparison なし** — 比較対象の条件を揃える (表で明示)
4. **負の結果を隠す** — 正直に書く方が accept されやすい
5. **Related work が自己引用だらけ** — 競合他者を 5 件以上
6. **Figure caption が 1 行** — 図の主張を caption 単独で伝える
7. **pragmatic な数値なし** — "significantly" だけ → "3.2x (Fig. 7)"
8. **cover letter で journal への fit が未言及** — "why this journal"
9. **response letter が箇条書きなし** — reviewer comment を line-by-line で
10. **supplementary material に本論の main result** — 本文に移す

---

## 💡 実例集 (汎用 template)

### Case 1: 国内和文論文誌 (IEEJ / 電気学会 論文誌 クラス)
- 構造: 和文 8 ページ、図 6 枚、引用 25 件
- Hook: 1 段落で「既存手法の限界 → 提案」を定量的に示す
- Contribution list: 3 件を intro 末尾に箇条書き
- Reviewer 対応: 反論 figure を必ず 1 枚追加する余裕を持つ

### Case 2: 国際英文誌 (IEEE Trans / IEICE Trans / IOP 系)
- 構造: 英文 6 ページ double-column、図 8 枚、引用 35 件
- Abstract: 200 words、problem→method→result→impact の 4 文構造
- Cover letter: "This work addresses a long-standing tradeoff in [domain]"

### Case 3: APS PRB long form / Elsevier JMMM
- 構造: 英文 12 ページ、図 10、引用 50
- 式番号 eq.(1), (2) を本文と密に結ぶ
- appendix を積極活用（derivation は appendix、本文は結果重視）

## ✅ 完成検証チェックリスト

### Abstract
- [ ] 4 要素 (問題/手法/結果/意義) が揃っている
- [ ] 単独で論文の主張が分かる
- [ ] 語数 (EN 200 / JP 400) 以内
- [ ] 弱気修飾語ゼロ
- [ ] 略語は初出で定義 (ただし広く知られた "FFT" 等は例外)

### Introduction
- [ ] CARS 4 段構造 (Territory / Niche / Occupy / Announce)
- [ ] Contribution list が末尾に 3 項目以内で箇条書き
- [ ] 関連研究を 15-25 件引用 (自己引用 < 20%)
- [ ] 本論文の構成が末尾に 1 文で宣言 ("Section II describes ...")

### Method
- [ ] 再現可能な level の記述 (パラメータ・設定値すべて)
- [ ] 新規性がある部分を明示 (既知手法との差分)
- [ ] 図で setup 図示

### Result
- [ ] 各 figure は独立して主張が伝わる caption
- [ ] 定量的数値 (誤差・標準偏差) 明記
- [ ] 比較対象との **fair comparison**

### Discussion
- [ ] 結果の **物理的 / 数学的解釈**
- [ ] Limitation を 1 段落 (正直に)
- [ ] 今後の方向性を示唆

### Conclusion
- [ ] Contribution を再確認
- [ ] 数値で締める (定性的な "effective" ではなく "4× speedup")

### Format
- [ ] 投稿規定のページ・語数制限内
- [ ] Overfull hbox ゼロ
- [ ] References が journal style に沿っている
- [ ] Figure は 300 dpi 以上、ベクター可
- [ ] Supplementary material が必要なら整理済

### 英文 (英語論文のみ)
- [ ] 時制統一 (Method 過去形、Result 過去形、Discussion 現在形)
- [ ] 冠詞 (initial reference "a", subsequent "the")
- [ ] Active voice を discussion / contribution で優先
- [ ] Native proof (IEEE Author Tools / DeepL Write / 校閲業者) 済

---

## 📖 推奨参考

- Swales, J., *Academic Writing for Graduate Students*. IMRAD + CARS 原典
- Tufte, E., *The Visual Display of Quantitative Information*. Figure 設計
- Strunk & White, *The Elements of Style*. 英文簡潔化
- Williams, *Style: Lessons in Clarity and Grace*. 英文 paragraph 構成
- IEEE Author Tools (editor 自動チェック) https://journals.ieeeauthorcenter.ieee.org/
- 木下是雄『理科系の作文技術』 中公新書 (和文 paper の古典)

---

## 🎯 アクセプト術 (Wallwork + 佐藤) 深掘りセクション

『日本人研究者のためのアクセプト術』(Wallwork, 514 p) + 『なぜあなたは論文が書けないのか』(佐藤雅昭, 182 p) の濃縮エッセンス。基本原則を超えた実戦ノウハウ。

### A. 執筆順序 (Wallwork §1.5, 佐藤 Q21)

```
figures → methods → results → intro → discussion → abstract → title
```

- 図は論文の骨格 → 最初に確定。欠けている panel = 追加実験
- Title は最後。本文が固まるまで論文の背骨は決まらない
- Abstract は penultimate (最後から 2 番目)。本文で主張が動けば書き直し

### B. "第一の創造" Phase 0 (佐藤 Q17)

本文を書く前の必須 artifacts:

1. **核となる datum** の特定 (Q15)。この結果があるから論文が存在する、という 1 図
2. **Figure-kamishibai** (紙芝居) を PowerPoint で構築 (Q17)。欠けているスライド = 追加実験
3. 研究に関わっていない 3 人以上に kamishibai を見せる (Q18)。内輪では当然の前提が他人には見えない
4. 関連文献 30 件以上 + モデル論文 3 本 (Q19-Q20)
5. Word テンプレ + 文献管理ソフト連携 を**先に**作る (Q21)

この Phase 0 を飛ばして書き始めた draft は棄却対象。

### C. Paragraph-level 規律 (Wallwork §14.9)

| 項目 | 目標 | 警告 | 致命 |
|---|---|---|---|
| 1 パラグラフ語数 (EN) | 75-175 words | >200 | >300 |
| 1 パラグラフ字数 (JP) | 150-400 字 | >500 | >800 |
| Introduction 総パラグラフ数 | 4-8 | <3 or >10 | 1 paragraph |

**長い Intro は弱い論文の signal** — 貢献が薄い著者ほど Intro を膨らませて隠す。Reviewer はこれを知っている。

### D. Given-New 情報配置 (Wallwork §3.4-3.6)

各文内で **既知情報を文頭、新情報を文末**。読者注意は文末で peak、文頭で renewal。重要な数値や発見を文中に埋めるな。

段落内では逆: **新情報を段落最初の文の最初、既知情報はあとで**。段落冒頭で読者を引き込む。

### E. Reviewer 2 対応 4 段階 (佐藤 Q40)

コメントを**返信前に分類**:

- **Difficulty A** (typo / 明快提案): そのまま直して "We thank the reviewer. This has been corrected (P11, L21)."
- **Difficulty B** (書き直し + 引用追加): 同意 → 根拠 → 改訂文を inline で示す
- **Difficulty C** (追加実験要求): paper が survive できるか判定。必要なら editor に延長申請
- **Difficulty D** (技術的に不可能な要求): 最も危険。必ず agree first → 理由 → 代替提案

### F. "Agree first, pivot" テンプレ (佐藤 p.167)

```
We agree with the reviewer in that [reviewer の懸念を寛容に restate].
Unfortunately, [理由 X, 技術的制約 Y].
Thus, we introduced [代替 Z] and demonstrated [結果].
The manuscript has been revised accordingly (Figure 8).
```

絶対に言うな: "The reviewer misread..." → 自分の説明不足として謝る:
```
We apologize for the confusion; what we meant in [section] is not 〜as the reviewer pointed,
but …. We have revised this part to carry a clearer message (P_, L_).
```

### G. 文法時制 3 部構成 (discussion) (佐藤 Q32 + Wallwork §14.8)

- **仮説** → 現在形: "air pollution **is** a stimulant"
- **自分の結果** → 過去形: "the distance **was** associated"
- **確立された背景** → 現在完了形: "X **has been shown** to ..."

この混同は**日本人著者最大のエラー**。

### H. Showing not telling (Wallwork §17.9)

**NG** (describe): "Figure 4 shows the relationship between A and B."
**OK** (claim): "The abundances of A and B were inversely related (Fig. 4)."

図番は**主語ではなく括弧内引用**。動詞は主張の核を表す。

### I. 禁断の Introduction 冒頭 (Wallwork §14.11)

以下のフレーズで始めたら書き直し必須:
- "Recent advances in..."
- "The last few years have seen..."
- "In this paper..."
- "X is an important field..."

**冒頭 5 語以内**に domain-specific noun を置くこと。

### J. Abstract 内 Background 配分 (Wallwork §13.16)

Background が abstract 全体の **25% 以下**。50% 超は Japanese-style で editor が即座 reject する典型パターン。

### K. Self-plagiarism の落とし穴 (Wallwork §11.10, §16.6)

過去論文の Method section を copy-paste すると自己剽窃。解決:
```
"Full details are given in our previous paper [ref].
In brief, [2-3 sentences summary + modifications]."
```

### L. 強調副詞の予算 (Wallwork §8.12)

`remarkably / notably / importantly / interestingly / novel / innovative / cutting-edge / crucial` の使用上限:
- 全体で **2 個以内** (Abstract + Intro で 1-2 個 max)
- 5 個超 → reviewer が "hype" と判定、信頼感低下

弱気修飾の禁止 + 強気修飾の節約 = プロの語彙コントロール。

### M. 数値 OK/NG ペア (Wallwork 実例)

**NG**: "Recent advances in the field of fragmentation..."
**OK**: "The physical process of fragmentation is relevant to several areas of science and technology. Because different physical phenomena are at work during the fragmentation of a solid body, it has mainly been studied from a statistical viewpoint [1-5]."

**NG** (past-tense 全般化): "Steel plates are treated with an open-flame burner."
**OK** (our specific action): "Steel plates were treated with an open-flame burner (Fig. 2)."

**NG** (意見を先):"The large difference between populations C and D is particularly interesting."
**OK** (事実→解釈): "While the mean size generally varies by only a few cm, populations C and D differed by 25 cm. Two hypotheses could account for this..."

### N. Discussion の A+B→C パターン (佐藤 Q33)

1 段落 = A (自分の結果の再要約) + B (関連文献) + C (2 つを合わせた洞察)。3-6 段落、各 100-150 words。

**NG**: "WH516 inhibits lung cancer. We showed this in our experiments."
**OK**:
```
WH516 has been demonstrated to inhibit gastric and breast cancer cell lines
in vitro and in vivo [refs]. Similar to these studies, we demonstrated that
WH516 inhibits a lung cancer cell line in a dose-dependent manner (Fig. 1).
These results support the reported direct inhibitory effect of WH516 on
cancer cells.
```

A+B→C が揃って初めて Discussion として機能する。

### O. Limitation 配置は Discussion の末尾ではなく中間 (Wallwork §9.12)

末尾 limitation = 防衛的。推奨: **主解釈段落の後、結論段落の前**。limitation を future work への橋渡しに使う。

### P. 反直感的な Tips

- **Long Introduction = 低新規性論文** (p.292)。Intro 積極的短縮は自らの contribution への自信の signal
- **Bullet は文法を揃えよ** — 一項目は動詞で始まり、もう一項目は名詞で始まるのは copy-editor flag (§8.5)
- **陰性結果も書け** (Goldacre, §17.7): 「あなたの方法が期待したものを測っていないと分かった」こと自体が findings
- **First draft 完了後のバーンアウト罠** (佐藤 p.155): supervisor 遅延を言い訳にするな。週次でリマインド、締切自己設定
- **学会発表の増加 = 論文の停滞** (佐藤 Q4): 会議の締切 → slides → 安堵 → 論文書かず ループ
- **Structured Discussion 見出し** (Wallwork §18.4): 医学系で普及しつつある "Principal findings / Strengths and limitations / Comparison / Meaning / Unanswered questions" は高水準 journal でも歓迎傾向
- **盲検 review での self-reference リーク** (Wallwork §7.11): "in our previous paper" は de-anonymize 最有力。"Doe et al. (2017) demonstrated" に書き換え、acceptance 後に revert

### Q. 拡張 チェックリスト

上記の定量的 target を統合:

- [ ] 1 パラグラフ 75-175 words / 150-400 字
- [ ] Introduction 4-8 段落、単一段落は失格
- [ ] Abstract 冒頭が禁断フレーズ 4 つに該当しない
- [ ] Abstract 内 background ≤ 25%
- [ ] Discussion の段落が A+B→C 構造
- [ ] Limitation が Discussion 末尾ではない
- [ ] 強調副詞合計 ≤ 2
- [ ] Figure caption が主張の動詞句で始まる
- [ ] Methods の tense が過去形 (我々の行為)
- [ ] Discussion の hypothesis が現在形、自分の結果が過去形
- [ ] Bullet の文法が並列
- [ ] Phase 0 artifacts 5 点揃い済

### R. 使える診断ツール (paper_writing_*)

既存 + 新規で計画されている:
- `paper_writing_check_subject_verb_distance` — 主述の物理距離
- `paper_writing_check_paragraph_length` — 段落字数の範囲
- `paper_writing_check_abstract_background_ratio` — abstract 内 background 比率
- `paper_writing_check_tense_consistency` — discussion の 3 部時制
- `paper_writing_check_figure_caption_showing` — caption が showing vs telling か
- `paper_writing_check_strong_adjective_budget` — 強調副詞の過剰
- **`paper_writing_check_word_repetition`** (v0.8.0) — 同一単語の近接障害 (中島・塚本)
- **`paper_writing_check_sentence_ending_variety`** (v0.8.0) — 文末表現の単調さ (中島・塚本)

---

## 📚 中島利勝・塚本真也『知的な科学・技術文章の書き方』(v0.8.0 で取込)

コロナ社 1996 年刊。日本工学教育協会賞・文部科学省特色 GP 採択教材。実験リポートから学術論文構築までの和文技術文章の教科書。**grant / paper 横断で適用可**だが、特に paper 執筆の「推敲時に気づかない悪癖」を潰すのに強い。

### §1.3.5 同一単語の近接障害 (本書の白眉)

**原則**: 狭い窓 (40 字程度) 内で同一単語が繰り返されると **幼稚な印象** を与える。語彙不足ではなく **意識の欠如** に原因がある。

**本家の例** (中島・塚本 p.37):
```
悪い: 最近の為替相場が円高傾向へ [大きく] 変化し、
      産業界特に自動車産業ではこの影響は [大きく]、
      圏内生産体制や部品調達システムなどが従来方式と比べて
      [大きく] 変化している。
```
→ 「急速に / 絶大 / 大規模」など類語で書き換え。

**もう一つの例**: 「...を測定することによって...の位置を測定する測定器」 → 「検出 / 決定 / 調べる」で置換。

**「の」 の使用回数**: 連続は **2 回まで**。`X の Y の Z の W の最大値` → `X と Y における Z の W の最大値` 等に分解。

**検出**: `paper_writing_check_word_repetition(text, window_chars=40)`。
**hint として使う** — 技術用語 (ESIM, POD 等) の反復は避けられないケースもある。全件機械置換は NG。

### §1.3.2 文末表現の変化

**原則**: 「...である。」「...である。」「...である。」と連続すると単調で幼稚。文末を **変化に富む場所へ落とす** ことで、格調の高い歯切れよい文章になる。

**本家の例** (p.26):
```
悪い: ...で「ある」。地球環境問題としては...が「ある」。
      ...予測を示す研究報告が「ある」。それは、...という報告で「ある」。
      ...最大の手段は、...代替フロンの研究開発で「ある」。
```

同じく「...した。...した。...した。」の過去形一辺倒も避ける。**実験方法の全文過去形 ≠ 再現性のある実験**。現在形・過去形・受動態を混ぜる。

**検出**: `paper_writing_check_sentence_ending_variety(text)` が histogram + Shannon entropy + 連続 3 回以上の critical_runs を返す。
**hint として使う** — 論理的に同じ形が必要な文もある。

### §1.3.1 接続詞の厳密な選択

接続詞の用法は 5 種類:

| 用法 | 接続詞の例 |
|---|---|
| 選択 | あるいは、それとも、または |
| 添加 | また、かつ、さらに、しかも、そして |
| 並列 | ならびに、かつ、一方 |
| 順接 | したがって、すなわち、それで、よって |
| 逆接 | しかし、だが、けれども、ところが |

**接続詞移動技法** (p.24): 「しかし」「そこで」を文頭でなく **文中** に挿入すると、歯切れの良い格調高い文章になる。ただし頻用禁止、単発的に。

### §6.4 校閲委員対応 (reviewer 返信)

査読者に対する回答文の作成法。
- **掲載可が予想される肯定的照会** → 素直に補足・拡張で応答
- **掲載否が予想される否定的照会** → 反駁より「前提の再確認」から入る (agreement first)

既存の `paper_writing_generate_response_letter` 等で運用。

---

## 関連 MCP

- `grant-writing` - 申請書 (科研費 / JSPS / KDDI / パワーアカデミー)
- `presentation` - 学会発表スライド (IEEJ SA / IEEE conference / 社内セミナー)

三者は共通の診断 tool (overfull hbox / sentence length / weak expressions) を持ちつつ、それぞれの文脈で閾値と推奨を変えている。
