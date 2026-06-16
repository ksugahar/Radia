# Presentation Skill — 学会発表スライドを通すための作文技術

研究発表スライドの作文ガイド。国内 (IEEJ 系 / 応物 / 機械学会) + 国際 (IEEE conference / APS March Meeting / Compumag 等) + 社内セミナー共通の設計原則。

---

## 🗣️ 良いスライドの定義 — しゃべりやすさ (speakability)

**良いスライド = 台本（セリフ）がスラスラ出てくるスライド。** 見た目の綺麗さ
ではなく、**しゃべり手が流暢に語れるか**で良し悪しを判定する。条件は 3 つ:

1. **セリフがスラスラ出る** — そのスライドを前に、言葉に詰まらず自然に話し
   始められる。詰まる = 主張が曖昧 / 情報過多のサイン。
   （木下 10章: 録音して「つまる箇所 = 筋が悪い」）
2. **繰り返しが少ない** — 前スライドで言ったことを再び言わされない。各スライド
   が **固有メッセージ 1 つ**（宮野 S4「1 枚 1 メッセージ」）。重複は「さっき
   言ったのに…」という失速を生み、最もしゃべりにくい。
3. **次スライドへのつながりが明確** — 今の締めが次の "つかみ" になり、送る瞬間に
   「では次は—」が自然に出る。タイトルを並べたら目次として筋が通る。

### 運用

- **判定**: tool を閉じて **声に出してリハーサル**。詰まった / 同じ事を 2 回言った /
  「えーっと次は」と言い淀んだ箇所 = 直すべきスライド（Phase 4 rehearse の中核）。
- **直し方**:
  - 詰まる → 主張を 1 文に絞る（情報を次スライドへ分割 or 削除）
  - 繰り返す → 重複スライドの一方に **別メッセージ** を持たせる（例: 概要スライド
    で手法を述べたら、表スライドは「誰が何を作ったか（分担）」に振る）
  - つながらない → スライド順を入れ替える / 各スライド末尾に次への 1 行ブリッジ
- **関連 tool**: `presentation_suggest_redundancy_fixes`（繰り返し検出）/
  `presentation_single_message_per_slide_semantic`（1 主張）/
  `presentation_slide_titles_outline_coherence`（隣接つながり・outline）/
  `presentation_script_vs_slide_coverage`（台本↔スライド対応）。

> **要点**: スライドは「見せる」より「**しゃべる**」もの。スラスラ語れる deck が
> 良い deck。見た目が綺麗でも、しゃべって詰まるなら作り直す。

---

## 🛠️ 診断ツール (Virtual Tools → Actual Tool Mapping)

v0.12.0 以降、Plan B Tier 1/Tier 2 の実測ツールで段階的にカバー。
各 virtual tool が「どの実装に接続されているか」を明示。

### T1. Opening hook 診断 → ✅ `presentation_opening_hook_strength` (v0.12.0)
最初 2 枚から Hook 強度 (数値/対立語/現場語/疑問投げかけ/汎用冒頭) を測定。
- **目標**: 「今日の話、聞く価値がある」と 30 秒で思わせる
- **NG**: "今日は○○について発表します" だけ (聞く理由がない)
- **OK**: 予期せぬ数値 / ビジュアル / ビフォーアフター を最初に出す

### T2. One-slide-one-message → 🟡 部分的 `presentation_check_bullet_count_per_slide` + 新規 T3 チャート検出
bullet 数で 1 主張違反を粗く判定 (bullet >= 4 で warn)。スライドタイトルの
動詞明示は `presentation_check_slide_title_verb` を補助利用。
- 目標: スライドタイトルを並べたら目次にして筋が通る
- **NG**: "結果と考察" (2 topic)、"○○について" (主張不明)
- **OK**: "Method X is 4× faster than FEM" (主張が動詞で明示)

### T3. Slide density → ✅ `presentation_check_slide_density` + v0.13.0 追加
`presentation_slide_density_balance` で deck 全体のバランス評価。
`presentation_visual_text_ratio_score` で画像/文字の面積比を算出。
- 高橋メソッド: 大字・1 主張・20 字以下
- Garr Reynolds: 画像重視・本文は 40-80 字

### T4. Visual flow → 🟡 部分的 `presentation_arrow_usage` + `presentation_check_image_text_ratio`
矢印の過剰使用を検出。視線誘導の完全自動評価は難しいが、矢印密度と
画像/文字比で粗く近似。
- 図表の色・矢印・枠で視線を導くか
- 複雑な図に appear animation で段階提示するか (pptx)

### T5. Speaking time → ✅ `presentation_estimate_speaking_time`
原稿から発表時間を逆算。slot 整合は手で確認。
- 日本語: 300 字/分 (学会), 400 字/分 (速め)
- 英語: 130 wpm (standard), 150 wpm (native fast)
- 持ち時間 20 分 → 6000 字 or 2600 words が上限

### T6. Q&A defense → ✅ `presentation_check_qa_backup_slides` + v0.12.0 takehome
Q&A 対応 slides の有無、takehome slide の質を評価。
- key claim に対する「なぜそう言える?」の data
- Method の limitation (正直に言及した方が好印象)

### T7. 制約検証 → ✅ `presentation_validate_pdf_pages` + `presentation_count_slides` + `presentation_check_time_13_rule`
発表 slot / 時間・ページ数の実測。

### T8. Print survivability → 🟡 部分的 `presentation_check_color_accessibility`
色の区別性検査。hatching / marker 形への変更推奨は mental。
- 色のみ区別 (赤/緑) は NG
- 凡例なし NG (プリントで色消える)

### 新規 (v0.13.0 Plan B Tier 2) — 追加の実測ツール

- `presentation_visual_text_ratio_score` — 画像面積 vs 文字面積
- `presentation_speaker_note_ratio` — speaker note 準備度
- `presentation_font_consistency` — フォント種類統一度 (宮野 S11)
- `presentation_arrow_usage` — 矢印過剰使用 (宮野 S10)
- `presentation_check_underline_in_pptx` — 下線壁紙化 pptx 版 (宮野 S14)
- `presentation_slide_density_balance` — deck 内 char 分布の Gini

### 既存 Plan B Tier 1 (v0.12.0)

- `presentation_opening_hook_strength` (T1 実装)
- `presentation_takehome_strength` (木下 p.230 / 宮野 S14)
- `presentation_check_pie_3d_charts` (宮野 S8 NG chart)
- `presentation_check_logo_on_every_slide` (宮野 S13)
- `presentation_check_progress_indicator` (宮野 S15)

### 和文 lint (v0.13.0 cross_lint re-exports)

台本・スライド原稿用:
notation_variants / find_undefined_acronyms / acronym_usage_audit /
check_kanji_ratio / lint_bedrock / check_misuse_japanese /
suggest_redundancy_fixes

### Plan B Tier 3 (v0.21.0) — outline / message / 欧米 vs 日本 / 構造

skill.md 既述だが未実装だった「タイトルだけ並べたら筋が通る」「タイトルに
最強主張」を実装 + 欧米 PPT text-heavy への日本理系 adaptation。

#### `presentation_slide_titles_outline_coherence(pptx_path)` — タイトル並びの論理整合
全 slide の title を抽出し並べた時に **目次として筋が通るか** を 5 軸で診断
(title_density / claim_density / 4 phase 出現 / 隣接 overlap / outline 完全性)。
extracted_outline で人間が音読判定可能。

#### `presentation_title_body_alignment_check(pptx_path)` — title が本文の最強主張か
既存 `check_slide_title_verb` が title の verb 形式のみ check に対し、本ツールは
**title vs body の内容対応** (token overlap + 数値裏付け + 5W1H) を per-slide で評価。

#### `presentation_single_message_per_slide_semantic(pptx_path)` — 1 主張 semantic
既存 `check_bullet_count_per_slide` (bullet 数のみ) を超えて、
**shape 種別の多様性** (text / image / chart / table の並存数) で多テーマ slide を検出。

#### `presentation_text_density_per_slide_western_style(pptx_path)` — 欧米 text-heavy 検出
**欧米 PPT 文字過多** を日本理系基準で flag。1 slide 文字数 + 段落型 textbox +
bullet 入れ子深さ + font サイズ の 4 signal で western-heavy slide を識別。
日本聴衆向けは ≤50 字/slide 推奨。

#### `presentation_mini_imrad_structure_check(pptx_path)` — 理系 mini-IMRAD 構造
Duarte sparkline (corporate / TED) ではなく、**理系プレゼンの mini-IMRAD**
(title→background→problem→methods→results→discussion→conclusion) 7 phase の
出現と順序を診断。section divider slide があると検出が安定。

#### `presentation_rikei_minimalism_score(pptx_path)` — 理系 minimalism (TED とは違う)
Reynolds Presentation Zen (TED 風) ではなく、**理系プレゼン特化の minimalism**
5 軸: figure or chart / axis label + unit / 定量 evidence / 出典 / bullet ≤4 + font ≥20pt。
domain peer 前提なので citation / 軸単位 必須。

#### `presentation_chart_simplification_check(pptx_path)` — Knaflic style
既存 `check_pie_3d_charts` (pie 限定) を多 chart 種に拡張。
3D bar/column/line + 凡例 ≥7 + 1 slide 3+ charts + chart title 不在 を flag。

### Plan B Tier 4 (v0.21.0) — 理系プレゼン特化 (TED 用ではない)

論文発表 / 学会発表 / 研究室セミナー想定の理系特化評価:

#### `presentation_equation_slide_compliance(pptx_path)` — 数式 slide 4 要素
**記号定義 + 物理意味 + 次元/単位 + 出典** (paper の数式記述と同じ厳格性)。
LaTeX 数式マーカー / formula 記号 / OMath を検出して数式 slide を識別、
4 要素 keyword detection で compliance 評価。

#### `presentation_figure_slide_compliance(pptx_path)` — figure slide 5 要素
**軸ラベル + 単位 + 凡例 + caption + 出典** (corporate プレゼンでは省略可だが
理系学会では必須)。chart 内部の axis title を python-pptx で取得試行、
image figure には text 推定で評価。

#### `presentation_results_slide_statistical_evidence(pptx_path)` — Results 統計報告 4 要素
paper T12 (statistical_reporting_compliance) の slide 版。
**error bar + sample size n + 統計検定 + effect size** が Results slide で
articulate されているか診断。理系 Q&A 予防の core。

---

## 📋 Workflow Phases

### Phase 1: 📖 **storyboard** — スライドを作る前にストーリー
A4 紙に手書きで筋を書く:
- Title: 聞き手が一番知りたいこと
- Why: なぜこの問題を選んだ
- Problem: 既存手法の限界
- Idea: 提案のひらめき
- Demo: 最もかっこいい結果 (Hook 用)
- Results: 定量比較
- Impact: 誰にどう役立つ
- Next: 次にやること

ここで 10 分以上悩む = まだ主張が固まっていない → 書き直し。

### Phase 2: ✍️ **draft** — スライド化
- タイトル行: 動詞句で主張 ("PEEC is 4× faster", "Fig shows X")
- 本文: 画像 > グラフ > 表 > bullet > 本文 の優先順位
- 1 枚 1 主張。要素が 2 つ以上あるなら分ける。

### Phase 3: 🔍 **diagnose** — 全体診断
- `presentation_count_slides(path)` で total 数
- `presentation_estimate_speaking_time(script)` で時間
- `presentation_check_slide_density(text)` で密度
- `presentation_count_weak_expressions(text)` で弱気修飾
- `presentation_check_overfull_hbox(log)` (beamer のみ)

### Phase 4: 🎤 **rehearse** — 声に出して練習
- ストップウォッチで時間測定 (初稿は 120% 時間かかる)
- 録音して聞き返す (つまる箇所 = 筋が悪い)
- 他人に見せる (1 人ラボの場合は鏡 or ペット)

### Phase 5: 🛡️ **prepare Q&A** — 想定問答
- 主張ごとに「なぜそう言える?」を書き出し backup スライド
- Limitation は先手で言及 (反撃を誘う)
- 「この論文知ってる?」の可能性 → 主要 3 本 + 比較 1 枚

### Phase 6: 🎬 **deliver** — 本番
- 最初 30 秒は話せることを暗記
- Slide 送り遅れ防止: clicker 左手、lab pointer 右手
- 時間超過は **致命的** (座長の評価が下がる)
- Q&A では **質問を復唱** → 答え → "Does that answer your question?"

---

## 📏 定量的閾値 (Quantitative Targets)

| 項目 | 目標 | 警告 |
|---|---|---|
| **スライド数 / 分 (学会, 木下)** | **0.7-0.8 / min** (10min → 7-8 slides) | **>=1.2 / min** |
| スライド数 / 分 (standard, セミナー) | 1-1.5 / min | >=2 / min |
| スライド数 / 分 (高橋) | 2-3 / min | >=5 / min |
| 1 スライド字数 (standard) | <=80 字 | >=150 字 |
| 1 スライド字数 (高橋) | <=20 字 | >=40 字 |
| 1 スライド箇条書き | <=5 項目 | >=7 項目 |
| 1 スライド図表 | 1 個 | >=3 個 |
| フォントサイズ (本文) | >=20 pt | <18 pt |
| フォントサイズ (タイトル) | >=32 pt | <28 pt |
| 文中の弱気修飾 | 0 | >=2 / slide |
| Overfull hbox (beamer) | 0 | >=1 |
| 発表原稿の字数 (20 min ja 学会, 木下) | **4800-5000 字** (240 字/min) | >=6000 字 |
| 発表原稿の字数 (NHK 基準) | ~6000 字 (300 字/min) | >=7000 字 |
| 発表原稿 words (20 min en) | ~2600 words (130 wpm) | >=3000 |
| **1 パラグラフ字数 (木下)** | **200-300 字** | >=400 字 |
| **1 スライド縦横比 (木下)** | **縦 2 : 横 3** (壁面プロジェクタ標準) | — |
| **1 スライド最大行数 (横)** | 8 行 | >=10 行 |
| **1 スライド最大行数 (縦)** | 12 行 | >=15 行 |

---

## 📚 基本理論

### A. 三大流派

| 流派 | 代表 | 特徴 |
|---|---|---|
| **高橋メソッド** | 高橋征義 (ruby 関西) | 超大字、1 スライド 1 word。資料不要の "漫才型"。 |
| **Zen プレゼン** | Garr Reynolds | 画像中心、本文最小、余白重視。TED 型。 |
| **横徹流** | 横山徹 (学会資料系) | グラフ + 結論 1 行。学会発表の標準。 |

研究発表では横徹流が無難。招待講演や企業セミナーでは Zen 寄り。
学生向け教育セミナーや TEDx では高橋メソッドが効く。

### B. スライドの情報階層 (top-down)

1. **タイトル** (画面上部): このスライドの **主張 1 文**
2. **キービジュアル** (中央): 図 / グラフ / 式。面積 50%以上
3. **説明文** (下部): 20-40 字で何が重要か補足
4. **footer** (画面下部): 発表題 / 発表者 / ページ番号 (印刷時 必須)

### C. 20-minute talk の黄金配分

| Section | 時間 | スライド数 | 目的 |
|---|---|---|---|
| Title + Self intro | 0.5 min | 1 | 発表者の背景 |
| Hook / Motivation | 1.5 min | 1-2 | 掴み |
| Problem / Related | 2.5 min | 2-3 | 既存の限界 |
| Proposed method | 5 min | 3-4 | 提案 (ここが主) |
| Result | 5 min | 3-4 | 定量データ |
| Discussion | 2 min | 1-2 | 解釈 / 限界 |
| Summary | 1.5 min | 1 | Take-home message |
| Q&A | 2 min (別枠) | - | 議論 |

### D. Key slide の 3 枚

どんな発表でも「この 3 枚だけ覚えて帰ってほしい」を設計:
1. **Problem slide**: 解きたい問題が 1 行で
2. **Result slide**: 定量的な key number が 1 つ
3. **Take-home slide**: 全体の要約 1 行 (slide 末尾)

これら 3 枚は特に推敲する。他の slide はこの 3 枚を補強する beam。

### E. IEEE / 国際学会固有の作法

- **Language**: English-only (日本語混在 NG)
- **Self-contained figures**: caption で図単独で意味が分かる
- **Eye contact**: スライド見すぎ NG、audience を見る
- **Accent**: 気にしない。スピードは遅めでOK (聞き取りやすさ優先)
- **Q&A**: 質問を復唱 ("So the question is..."), 答えは 30 秒以内

### F. 動画収録対応 (online / hybrid)

- レーザーポインター使用不可 → pptx の highlight / beamer の \alert
- アニメーション は控えめに (video 圧縮で崩れる)
- 音声: 外付けマイク必須 (ノート PC 内蔵はノイズ)
- 背景: 無地または blur
- 画面共有: 1920×1080 で書き出し

---

## 🗂️ OK/NG スライド例

### Title slide

**NG**:
> 「誘導加熱の数値解析について」 (主張なし)

**OK**:
> 「提案手法: メッシュ非依存な渦電流 simulation」 (手法 + 特長)

### Problem slide

**NG**:
> ・従来手法は遅い
> ・精度に問題がある
> ・改善が必要

**OK**:
> 「高周波 eddy current 問題を FEM で解くと 1 ケース 6 時間。
>  パラメータ sweep 50 点で 2 週間。それを半日にしたい。」
> (数値で痛みを具体化)

### Result slide

**NG**:
> 結果 (グラフタイトル)
> [軸ラベルなし折れ線グラフ]

**OK**:
> 「提案手法は FEM 比 4× 高速、精度差 <3%」
> [FEM vs 提案の実時間棒グラフ、赤矢印で 4× 明示]

### Summary slide

**NG**:
> ・手法を提案した
> ・評価を行った
> ・今後の課題がある

**OK**:
> 「**Take-home**: 提案手法で解析を 4× 加速。
>   製品設計の DoE が 1 日で回せる」
> [主要 Figure のミニ版]

---

## 🚫 NG パターン 10 (発表特有)

1. **"今日は...について発表します" から始まる** — 退屈。結果や数値から開始
2. **スライドにフルセンテンス 3 行** — 読まれる / 聞かれない
3. **同じ bullet が 7 項目以上** — 覚えられない (miller's 7±2)
4. **凡例なし grpah** — 色消える / 意味不明
5. **Method 詳細を 10 slide** — 聞き手は method に興味なし。結果で勝負
6. **弱気語 on key slide** — 信頼感低下
7. **Summary が箇条書き 5 項目** — 1 行 take-home に絞る
8. **Q&A backup slide がない** — 準備不足を露呈
9. **動画が再生できない** — 必ず現場の PC で事前確認
10. **時間超過** — 座長評価が最も下がる。練習で 90% 時間に収める

---

## 💡 実例集 (汎用 template)

### Case 1: 国内学会 20-min talk
- slot: 20 min (15 slide + 5 Q&A)
- 全 15-18 slide が木下推奨上限 (0.7-0.8 slide/min)
- Key slides: Slide 2 (Hook = 数値的 Problem), 中盤 (4× 改善 等の key number), 末尾 (Take-home 1 行)
- 原稿 4800-5000 字 (240 字/min × 20)、リハで 10% 余裕

### Case 2: IEEE conference 12-min talk
- slot: 12 min (10 min talk + 2 min Q&A)
- 全 10-12 slide (タイトル + 結論含む)
- English-only、Hook は figure でインパクト重視
- Q&A backup 5 slide (method detail, failure mode, cost, scalability, future)

### Case 3: 招待講演 (60 min)
- slot: 60 min (45 min talk + 15 min discussion)
- 全 30-40 slide (Zen 寄り、図 70% / 本文 30%)
- Take-home を冒頭と末尾で同じ 1 行

## ✅ 完成検証チェックリスト

### Storyboard
- [ ] 1 行 take-home message が slide 末尾に 1 つ
- [ ] Problem slide で数値的痛み
- [ ] Result slide で数値的 key

### Structure
- [ ] Title / Hook / Problem / Method / Result / Discussion / Summary
- [ ] 各 section のスライド数が 20-min 配分と整合
- [ ] 1 スライド 1 主張

### Typography
- [ ] タイトル >= 32 pt, 本文 >= 20 pt
- [ ] 太字 / 下線 は key number のみ
- [ ] 弱気修飾語 ゼロ (contribution / result)

### Visuals
- [ ] 図表に軸ラベル / 凡例 / 単位
- [ ] 白黒印刷で読める (hatching / marker 形)
- [ ] 視線誘導 (矢印 / 枠) が設計されている

### Timing
- [ ] スライド数 x 1.5 min <= slot
- [ ] 原稿字数 <= 持ち時間 x 300 字/分 (ja) or 130 wpm (en)
- [ ] リハーサル 1 回は実施

### Q&A
- [ ] Backup slide 3-5 枚 (method detail / limitation / cost / future)
- [ ] Limitation を主発表で 1 枚言及
- [ ] 英語発表なら Q&A の出だし "Thank you for your question" 暗記

### Format
- [ ] PDF export で崩れない (フォント埋込)
- [ ] 発表 PC (会場) のフォント確認
- [ ] ページ番号 footer 表示

### Delivery
- [ ] 最初 30 秒暗記済
- [ ] Take-home を冒頭で 1 回、末尾で 1 回
- [ ] 時間 10% 余裕で終わる
- [ ] Clicker / pointer 動作確認

---

## 📖 推奨参考

- 高橋征義 『高橋メソッド』http://www.rubycolor.org/takahashi/
- Garr Reynolds, *Presentation Zen*. Zen style の原典
- Nancy Duarte, *slide:ology*. 視覚設計
- 木下是雄 『理科系の作文技術』 中公新書 (普遍)
- Patrick Winston's "How to Speak" MIT lecture (YouTube) — 必見
- IEEE Presentation Guidelines: https://www.ieee.org/conferences/

---

## 🎨 作図力学 (Sakuzu Rikigaku) — 学術図の作り方

『作図力学』(54 p) から抽出した、発表図の視覚設計規則。論文図にも適用可。

### 視覚の weight 階層 (pecking order)

1. **Marker 迫力**: `●` > `■` > `▲` > `◆` > `○` > `□` > `△` > `×`
   - 黒ベタ > 面積大 > サイズ大、ほど迫力
2. **Line 迫力**: solid thick > solid thin > dashed > dotted
3. **Frame は plot 線より太くするな** — 幼稚な印象
4. **Curve-frame 比で主張を encode**:
   - 理論と実測一致 → 理論曲線を太く
   - 不一致 → 理論曲線を細く
   - 複数手法比較で best 強調 → その曲線だけ太く

### 軸の選択 (Sakuzu p.4-7, p.27-28)

- **曲線 vs 折れ線**: 縦軸量の変化が連続 かつ 散布度小 → 曲線。それ以外 → 折れ線
- **dot plot (線なし)**: 密度が高く点群だけで trend 見える場合のみ。点は接するくらい、overlap はしない
- **3D は最後の手段** — 定量読み取り困難、shading 必須
- **縦横比**: 縦長は散布度を誇張。横長は連図 (panel series) の場合のみ可
- **原点マーク 0 は単一** (x=y=0 が共通の時、角に 1 つだけ)

### ラベルの rule

- **量記号はイタリック**、subscript は真の下付き (K_a、Ka ではない)
- **単位の書き方**: `J/(kg·K)` は OK、`J/kg/K` は禁止 (両義)
- **y 軸ラベル**: 下から上に書く (横向き文字を縦に配列)
- **実験条件は plot 内の空白に** — caption 下ではない。右上角に整列

### 系列の区別 (color は最後)

優先順位: **marker 内形状** (○ ● ◐ ⊕) → **shape** (○ □ △) → **line style** (solid/dash) → **color**

- 旗手の系列に `×` 使用禁止 — 粗雑な印象
- 平滑曲線に `□` 使用禁止 — 上下端の平面が曲線の連続性を壊す
- 最終印刷サイズで marker size 決定。卒論の図を学会で再利用するな (縮小率が違う)

### Grid / scatter 表示

- **Grid は原則 OFF**。「締まりのない」と感じた時だけ点線 grid を使う
- 散布度大時の trick: frame 太く + smoothing 曲線厚 + envelope 細 → 「散布度の圧縮」
- 生データ ("raw") 図を意図的に 1 枚入れる = 信頼性アピール

---

## 🖼️ 図・表・キャプション 追加規則 (mined from 8 書籍, 2026-06)

既存の「🎨 作図力学」節（marker/line/frame weight 階層・量記号イタリック・
`J/(kg·K)`・y 軸下→上・実験条件を余白へ・色は最後・縦横比・原点 0 単一・
生データ図・3D 最後の手段）と、`radia_mcp.figure` モジュール（title 禁止・
単位は括弧・side-by-side 4 cm・10pt@8cm・TikZ 出力・凡例 frame/overlap 禁止・
colorblind-safe・font 埋込）でカバー済みの規則はここでは繰り返さない。
以下は **それらに無い、または大幅に拡張する規則** のみ。論文図・学会図・
スライド図に共通で適用する。

---

### A. 数値・単位のタイポグラフィ（強調と SI）

- **単位は数字の約 6 割サイズにする** — 単位を数字と同サイズで書くと主役の
  数値が埋もれる。`300%` なら `%` を数字の ~60% に縮小し、数字本体は太めの
  書体にして主役を明確化する。（スライドデザイン p.102）
- **強調したい数値だけ太字にし、必要なら付加記号で意図を添える** — 数値は
  メッセージの核。`1000万kW` のような重要数値のみ太字化、危機感を出すなら
  `供給力不足!!` のように記号を必要時だけ付す（乱用しない）。（スライドデザイン p.102）
- **数値と単位の間は 1 字あけ、分母が積なら必ずカッコでくくる** — `3.0 kg`
  (○) / `3.0kg` (×)、`W/(m²·K)` (○) / `W/m²·K`・`W/m/K` (×)。くっつくと
  一語に見えて誤読される。（木下 §9.3.1 / 作図 §4.3.2 ⑪）
- **単位記号は SI を用い、固有名詞由来は記号 1 文字目だけ大文字・単位名は
  全小文字** — `Hz・N・Pa・W`（記号）/ `hertz・newton・watt`（名）。`sec` は
  不可で `s`、`Torr・kgf・cal・μ(ミクロン)` は避ける。（木下 §9.3.1）
- **軸ラベルは「量／単位」の無名数方式で書く** — スケールの数字が
  (量 ÷ 単位) のディメンションのない数になるよう、`V / (10³ m·h⁻¹)`・`t / h`
  のように書く（旧来の「軸上部に ×10³」方式より明確）。
  ※IEEE/IEEJ では `f (Hz)` の括弧方式が標準（`figure` モジュール参照）、
  無名数 `f / Hz`(ISO 80000) は別系統 — 投稿先の慣習に従い混在させない。（木下 §9.6）
- **スケール数字に不要な有効桁を付けない** — 整数で割れるなら `5,10,15`。
  `2.5,7.5` 等の小数を必要なときだけ併記する。（作図 §4.3.2 ⑤）

### B. 図番号・キャプションの体裁

- **図番号＋キャプションは図の下に置き、枠で囲まない** — 図内タイトルは不要
  （title 禁止規則のとおり）、説明は図の下が標準。英文は `Fig.` を用い先頭
  大文字。（作図 §4.3.1 / 知的論文 §4.3）
- **複数行キャプションは 2 行目以降をキャプション本文の行頭にそろえる** —
  図番号の位置にぶら下げない（ずれて読みにくい）。表のタイトルも同様。
  （作図 §4.3.1 12) / §5.1）
- **キャプション中に量記号を入れて図の内容を説明する** — 例
  「図3 比熱 *C* と熱膨張量 *Ka* を考慮したときの関係」。本文を読まずに何の
  関係を示す図か分かり、図の「独り歩き（転載）」に耐える。（作図 §4.3.1 12)）
- **座標軸線に矢印を付けない（矢印つきデカルト線図にしない）** — 高校数学・
  純粋数学の慣習であり、理工系の卒論/修論/学術論文では初心者っぽく見える。
  軸は枠線（長方形）で囲む。例外は枠が描けない生データ図のみ。
  （作図 §4.3.1 / 知的論文 §4.3.1）

### C. キャプション／凡例の自己完結・showing not telling

- **図表とキャプション（legend）は本文なしで単独理解できるよう、タイトル並みに
  濃い情報を載せる** — グラフが示す結果・研究対象・処理の背景・場所・解釈に
  必要な説明・条件（温度/培地）・サンプルサイズと統計検定結果まで。
  Nature 系は本文語数制限のため方法の大半を legend / 脚注に書く。
  （Accept 英語論文 §8.7 / §17.11 / なぜ論文 Q28）
- **legend には「図を見れば分かること」でなく解釈・意義を述べる（showing not
  telling）** — `Figure 4 shows the relationship between A and B` のような図示
  内容の反復は付加価値ゼロ。`The abundances of A and B were inversely related
  (Figure 4)` と書く。図表参照番号は括弧に入れ文末へ。（Accept 英語論文 §17.9-17.10）
- **本文で図表に言及するときは能動態・簡潔に、`Fig.` は略す/`Table` は略さない**
  — `This figure shows X`（○）/ `X is shown in this figure`（×）。
  `graphically/schematically` 等の冗語や `can be seen` を避ける。legend 内では
  Fig も Table も略さない。（Accept 英語論文 §5.16 / §17.11）
- **論文に載せた図表はすべて本文で言及し、示す結論・傾向を述べる** — 言及されない
  図表は宙に浮く。`Table 1 shows the summary` のように意識を向けるだけの文では
  なく、傾向・差を述べる。（Accept 英語論文 §17.11）
- **同じデータを図と表の両方で示さない／図表は重要なものだけ載せ相互に重複させ
  ない** — 図表枚数の上限を投稿先で確認し、各図表が異なる情報を持つようにする。
  「重要なものだけ載せるほどメッセージはクリア」（Nature）。（木下 §9.6 / Accept §8.7 / §17.14）
- **legend に矢印の意味・スケール・略語・統計値を入れ、本文との重複は一方に寄せる**
  — `arrows indicate ...`、scale bar、略語定義、`r=0.98, P<0.001` 等。図を見れば
  分かる差は本文か legend の一方だけに。（なぜ論文 Q28）
- **回覧用に図表＋legend を 1 対 1 で 1 ファイルにまとめる** — 本文・図・legend が
  別々だと読み手（指導教員/Reviewer）が行き来する負担を負う。1 ファイル化は
  読んでもらえる確率を上げる。（なぜ論文 Q38）
- **図表内の数値・単位・通貨を比較可能にそろえる** — 年収（$）と月収（€）の混在、
  通貨記号の有無不統一、非対称な母集団比較は読者を混乱させる。（Accept 英語論文 §17.13）
- **データの価値は「興味深い」と言わず具体値で示す** — `The large difference is
  interesting`（×）→ `populations C and D differed by 25 cm`（○）。読者自身が結論に
  達せるようにする。（Accept 英語論文 §17.8）

### D. 投稿・印刷向けの出力仕様

- **図はカラム幅に合わせて設計する** — 1 カラム ≈ 幅 8.7 cm (3.4 in)、2 カラム
  またぎ ≈ 18.1 cm (7.1 in)、高さ最大 ≈ 25 cm。作る前に投稿先のカラム幅を確認
  する。（なぜ論文 Q28）※lab の `figure` profile（IEEE 88.9/182 mm）と同主旨。
- **投稿規定の解像度・形式を満たす（典型 300 dpi・TIFF）** — PowerPoint 既定
  96 dpi のままでは不足。スライド幅を `必要dpi×幅cm÷96` に拡大して出力する
  （300 dpi・1 カラム → 約 27.2 cm）。線画 1200 dpi / 写真 300 dpi / 混在 600 dpi
  等の指定に従う。グラフは Excel/PPT 切貼りより**イラストソフトでトレースし直す**
  と線のガタつきが消える。（なぜ論文 Q28）
  > ⚠️ ただし IEEE/IEEJ/IGTE のエンジニアリング系は **vector PDF / TikZ** が原則
  > （`figure` モジュール「camera-ready に raster PNG を埋めない」）。300 dpi TIFF
  > 指定は主にライフサイエンス/Nature 系。**venue で要求形式が分かれる** — 投稿規定
  > を最優先。
- **元画像（写真等）は最初の保存段階で解像度を最大限上げておく** — 後から図全体の
  dpi を上げても元素材の質以上には改善できない。（なぜ論文 Q28）
- **カラー画像は RGB/CMYK 指定に従う** — 画面 RGB・印刷 CMYK。規定が CMYK なら
  イラストソフトで CMYK プレビューしつつ変換する（ポスターも同様）。（なぜ論文 Q28）
- **写真にはスケールを書き込み、刷り上がりの 1〜1.5 倍の印画を用意する** — 印刷所では
  1.5 倍印画が最もきれいに仕上がる。線画・本文は刷り上がりの 2〜2.5 倍に作図し、線は
  太め・字は大きめにする（縮写後にちょうどよくなる）。（木下 §9.6）
- **図に完璧を求めすぎない（誌の grade × 自分の重要度でかける手間を調整）** — 図は
  凝りだすとキリがない。全てに完璧を求めると生産性が落ちる。（なぜ論文 Q28）
- **普通紙印刷の写真は印画紙より解像度が落ちるので、掲載写真で確認できない微細形態を
  本文で過剰に論述しない** — 読者は掲載写真しか見られない。（作図 §4.2.1 実写図）

### E. 装置図・解析モデル図・可視化（信頼性 = 認識時間の短さ）

- **作図法の極意は「認識時間の短い図」を作ること（図の信頼性は認識時間に反比例）** —
  一見で理解できる図ほど信頼・評価される。図種選択・作図力学・可視化・矢印ラベルを、
  すべて読者の認識時間を最短化する観点で評価する。（作図 §4.5.1）
- **装置図・解析モデル図は線画イラストで描き、不要部分を省略し主張部分を明示する
  （写真にしない）** — 装置図は読者が最初に精査し論文全体の印象を決める。写真では
  不要情報を省けず主張部位も明示しにくい。（作図 §4.2.1 / 知的論文 §4.5.2）
- **装置図の具備条件は (1) アピール性 (2) 現実感 (3) 明確さ** — 陳腐な平面長方形
  構図でなく、独創的構図で「良い結果が出そう」と期待させる。単なる白抜き長方形で
  並べない。（作図 §4.5.2 / 知的論文 §4.5.2）
- **図の一部に写真を使うなら他の装置部も同等の詳密さで描き、リアリティをそろえる** —
  隣が単純な四角枠だと写真の現実感まで損なう。（作図 §4.5.2）
- **スクリーンシート（網掛け）は装置図には多用してよいが、線図では多くても 2 箇所まで**
  — 線図で多用すると強調点がぼけて逆効果。装置図は質感・立体感が出て有効。網掛けは
  縮小で潰れるので「大きすぎるかと思うほど」粗く選ぶ。（作図 §4.5.3）
- **肉眼で観察できない現象は可視化図面（3D 表示＋隠れ線処理＋軸の誇張拡大）にする** —
  詳密なほど写真的印象を与え脳が瞬時に理解する。電子顕微鏡等で示せない面性状などを
  PC で 3D 化・隠れ線処理し、必要なら Z 軸を誇張拡大する。（作図 §4.5.1）
- **生データ図でも枠の役割を確保する** — 補助スケールだけの表示は失敗しやすい。枠線が
  描けないときは実験条件などの説明文字を**四隅**に配して領域を確定する（振幅が大きい
  ときは振幅中心に白線を入れると変化が明確）。（作図 §4.5.1）
- **プロット点記号の内部を平滑曲線が貫通してはならない** — 白抜き記号の中を曲線が突き
  抜けると記号内部の模様が読めず識別性が落ちる。曲線をプロット点手前で止める（隠れ線
  処理）か、黒ベタ記号を使う。（作図 §4.3.2 ⑦ / 知的論文 §4.3.2 ⑦）
- **フローチャートは誤解されない範囲で最大限ブロックを集約する** — 1 万行のプログラムも
  誤解なく 12 ブロックに集約できる。複雑な手順は「図面フローチャート技法」（複数の図を
  矢印で連結して流れで示す）を使う。（作図 §4.2.1 / §4.5.4）

### F. 表のデザイン（罫線最小・強調・配置）

- **表は罫線を最小限にする** — 縦横すべてに罫線を引かず、両端の縦罫を付けず左右を
  開放する。最上段（見出し下）だけ太罫（0.2〜0.3 mm）、他の横罫は細く（0.05〜0.1 mm）。
  行ごとの罫線は引かない。（作図 §5.1 / スライドデザイン p.124）
- **表の共通単位は毎セルに書かず、表全体共通なら右上、列共通なら列見出しに入れる** —
  各セル併記は表内がごちゃごちゃする。（作図 §5.1 / スライドデザイン p.127）
- **表でも注目データを文字色・セル色・吹き出しで目立たせる** — 数値の羅列では「どこが
  要点（例：従来比 2 倍）」か伝わらない。表内でも目線誘導が要る（網掛けアクセントは
  1〜2 箇所に限定）。（作図 §5.3 / スライドデザイン p.124）
- **表番号・タイトルは表の上部中央に置く（図は下、表は上）** — 図と表で番号位置の慣習が
  逆。英文は `Table`（略さない）。（作図 §5.1）
- **表中の文字は本文より小さく、図中の文字と同ポイントにする。補足説明はさらに小さく** —
  本文と同サイズだと不自然。（作図 §5.1）
- **説明表は本文より簡潔に要約し、各項目の説明は 2 行までに収める** — 大分量の説明表は
  かえって理解しにくい。冗長な記述は本文へ回す。（作図 §5.2）
- **罫線の代わりに網掛け（ケイ線省略法）で表構造を示せる** — 見出し行を濃い網＋白抜き
  文字、データ行を薄い網にすると罫線なしでも明確で本文と区別できる（横項目 3 列以上には
  不向き）。装置形態が要る表にはイラスト/写真を組み込むとアピール性が上がる。（作図 §5.3）

### G. 作成プロセス（図を先に・紙芝居でストーリー）

- **本文を書き始める前に図・表を準備する（フリーハンドでよい）** — 理科系の文書では図表
  が最も大切な役割を演じることが多く、図表を先に作ると「何を書くか」が明確になる。立案
  段階では印刷用に仕上げる必要はない。（木下 §2.5.3 / §9.6）
- **執筆前に核データ中心の「Figure の紙芝居」を作り、図の並びでストーリーを設計する** —
  ① 核データの図を作る → ② それを中心にストーリーを考える → ③ 前後に足す図を置く。
  虫食い部分＝これから集める材料。図の順序がそのまま Results のストーリーになる。研究室
  ミーティングで紙芝居を見せ筋が通るか意見を求める。（なぜ論文 Q17/Q18）

### H. 写真・スライド図・図中ラベルの細則

- **写真は伝えたい箇所を中心にトリミング・拡大し、特定箇所は丸囲み・矢印で指す** — 全体
  写真のままでは見せたい部分が小さく埋もれる。（スライドデザイン p.130）
- **図中の凡例・ラベル・注記でも修飾語を被修飾語に直結し、二義を避ける** — 図ラベルは短い
  名詞句が多く係り先が二義になりやすい（「強い価格面以外の競争力」→「…の強い競争力」）。
  並列項目はテンでなくナカテン（A・B・C）で区切る。（本多『日本語の作文技術』）
- **図キャプションの引用・出典は原文に忠実にし、引用範囲をカギで厳密に区切る** — 出典を
  解釈混じりで書くと改竄になる。自分の言い換えはカギを外す。（本多『日本語の作文技術』）

---

## 🎤 木下の講演テクニック (『理科系の作文技術』10章)

### 時間配分

- **1/3 則 (comprehension gradient)** — 前 1/3: 全員分かる話、中 1/3: ほとんどが「分かった気がする」、後 1/3: 専門家向け
- **1/4 則 (section 時間)** — 10 分講演の場合: 序論 2.25 + 方法 2.25 + 結果 2.25 + 議論 2.25 min

### 原稿と memo

- **である体・短文** で原稿作成 → memo card 抽出 → 登壇は memo のみ
- 原稿を読み上げるな — 「ついていくのは非常な努力」
- memo card: **A6 (15×10.5 cm)、片面のみ**、スライド番号は**色枠で囲む** (薄暗い場で目立たせる)
- 机上に左から右へ順に並べる (Q&A で瞬時に参照)

### 話し方

- **丁寧すぎはNG**: 「報告させていただきます」「であります」 → 「します」「です」
- **ズバリ話法**: 「〜ということのようです」 → 「〜とわかりました」
- **黙る technique**: 注意が散った聴衆に再注目させるには、声を大きくせず**2-3 秒沈黙** → 「つまり...」で key line

### スライドのハンドリング

- **「戻る」指示は混乱の元** → 必要なら slide を**複製**して後方にも置く
- **Q&A backup 方程式 slide** を projection 係に事前預け
- 添字は**大きすぎるくらいで丁度よい** — プロジェクタ投影で縮小される

### 道具

- マイクは **必ず外して左手**に持つ (スライドに振り向いても途切れない)
- 前の発表者の本番時に道具 (clicker / pointer / mic) を確認

---

## 🎨 『研究発表のためのスライドデザイン』の取込 (v0.9.0 OCR 取込)

**出典**: 宮野公樹『研究発表のためのスライドデザイン』(講談社, 2010; OCR 取込 2026-04-23)。image-only PDF を OCR で取込、下記の原則集として整理。

### 第 1 部 「わかりやすい」スライド構成にするために

**S1. 作成前に全体構成を検討**: スライド作成を始める前に、**内容を精査し各スライドに掲載する情報を明確化**する。ホワイトボード + 付箋で要素を書き出し、論理展開をフローチャートで可視化してから着手する。

**S2. オリジナリティ/客観性 self-check**:
- オリジナリティ: 問題設定 / アプローチ / 結果のうち **どこに独創性** があるか明言できるか。類似研究との違いが言えるか。
- 客観的・論理的か: 背景説明が聴衆のレベルに合うか。**「何の問題をどう解くか」を一言で**説明できるか。
- データ完全性: 単位は適切か、再現性確認したか、飛躍していないか。事実紹介で終わらず考察まで到達しているか。

### 第 2 部 「わかりやすい」スライドを作成する技術

**S3. Kiss の法則 (Keep It Simple & Short)**: メッセージを短く単純に。

**S4. 「1 枚に伝えたいメッセージ 1 つ」**: 4 つのメッセージがあれば 4 枚に分ける。

**S5. コントラストの活用** (強調の 3 技法):
- 文章中のキーワードに色を付ける
- 重要な一文の書体 (フォント) を太くする
- 図中の重要箇所のみに色を付ける

**S6. 同じ種類の要素を縦方向に揃える** → スライド全体がスッキリ。

**S7. 行頭記号の原則**:
- **大きめの行頭記号を乱発すると ごちゃごちゃする**
- 項目の段階に合わせて行頭記号で上下関係を表す
- 行頭記号が目立ちすぎないようにする

**S8. グラフの選択**:
- 円グラフは細かい角度比較しにくい → **帯グラフ / 横棒グラフ** を優先
- 項目数が多い → 横棒グラフ (項目名が表示しやすい)
- **Excel のグラフを初期設定のまま使わない**
- 3D グラフを安易に使わない
- 軸目盛りが細か過ぎない / 凡例の位置がグラフから離れていない
- グラフの伝えたい箇所へ目線を自然と導く

**S9. 強調技法** (視線誘導):
- 注目させたい棒の色を変更して強調
- 中心から少しずらして強調
- セル / 文字の色を変更
- 矢の色を変えて目線を自然に導く

**S10. 矢印の使い方**:
- 三角形部分の長さに注意
- 直線が図と重なる場合は 点線 / 曲線 に
- 矢印を縁取りしない

**S11. 文字・書体の原則**:
- **文字サイズは 40 pt 標準** (body)
- 書体は適材適所、**全部同じにしない** (メリハリ)
- 汎用性の低い書体を使わない
- **丸ゴシック体を使いすぎない**
- 派手な文字装飾機能を使わない
- すべて同じ書体 = メリハリなし

**S12. 配色の原則** (3 色使い):
- **ベース色・メイン色・アクセント色を使い分ける** (計 3 系統が基本)
- 色数が多すぎない
- **スライド間で色使いに統一感**
- 無意味なグラデーションを使わない

**S13. 写真・イラスト・図形の原則**:
- イラストのテイストに統一感
- 伝えたいことと直結していないイラストを入れない
- **全スライドにロゴを入れない** (聴衆の目線を惑わせる)
- 角丸四角形の角の丸さを揃える
- テンプレ集のぎこちない図形を使わない
- 楕円を多用しない
- 背景写真に依存しない (主張が背景に見合うか)

**S14. 下線・装飾の原則**:
- **長文の文字列すべてに下線を引かない**
- 点の下線を引かない
- 直線ばかりが目立つスライドにしない
- 円弧など適度に使い分ける

### 第 3 部 スライド全体の構成を聴衆に伝える工夫

**S15. 進行位置の明示**: 目次項目を各章の冒頭に表示し、**現在位置を色で示す** (「今どこを話しているか」を聴衆に持続的に伝える)。

### Tool 対応 (v0.9.0 新規 / 既存)

| 原則 | 対応 tool |
|---|---|
| S11. 40pt body | `presentation_check_pptx_font_size(min_body_pt=20)` ← 宮野基準なら 40pt に上げる |
| S12. 色数制限 | **`presentation_check_color_count_per_slide`** (v0.9.0 新規、3-5 色推奨) |
| S14. 下線密度 | `presentation_count_underlines` (tex 向け、slide の下線は目視) |
| S7. 行頭記号過多 | `presentation_check_bullet_count_per_slide(max_bullets=5)` |
| S9. 強調過多 | `presentation_check_color_accessibility` + 色数 tool |

**⚠️ anti-Goodhart**: 上記 S1-S15 はすべて **hint**、机上の数値を 0 化することが目的ではない。発表の内容 / 聴衆 / 会場サイズで最適解は変わる。最終パスは **tool を閉じてスライドショーを実演**、声に出して 12 分 / 15 分を測って違和感を潰す。

---

## 📐 Workflow Phase 2.5: 5 段階 figure pipeline (作図力学)

storyboard → draft の間に figure 作成の正式プロセス:

1. **下図** (pencil sketch) — 縦横比・サイズを手書きで確定
2. **基本図** (software trace) — 軸・枠・marker のみ、データはまだ
3. **作図力学 tuning** — weight 階層を別 pass で適用
4. **推敲 pause** — 完成後しばらく見ないでから再チェック
5. **実印刷サイズプレビュー** — 最終縮小率で subscript 判読可能か確認

---

## 📚 追補: 6 冊からの新規ルール (2026 mined)

既存の宮野 (S1-S15) / 木下 (10章) / 作図力学 と重複しない、**新規または実質的に拡張された**規則のみを以下に集約する。出典は b01『研究発表のためのスライドデザイン』(宮野系) / b02『作図力学』/ b03『理科系の作文技術』(木下) / b04『なぜあなたの研究は…』/ b05 英文 academic writing / b06 知的プレゼン (OHP系) / b08 本多『日本語の作文技術』。

> **⚠️ すべて hint**: 机上の数値を 0 化するのが目的ではない。最終パスは tool を閉じて声に出して実演し、違和感を潰す (anti-Goodhart)。

---

### 追-A. タイポグラフィ / レイアウトの新規定量 (宮野 b01)

既存 S5-S14 を補う、**未収録の具体寸法・余白規則**:

- **行間 = 文字サイズの約 1.3 倍。** 0 倍は詰まって読めず、2.0 倍は離れて散漫。1.3 が最も読みやすい。あわせて行数を減らし 1 枚の文字量を抑える。(宮野 b01 p.92)
- **ハイライト (強調色) は画面面積の 10% 以下。** 色付け・太字・下線の合計が 10% を超えると、かえって何も目立たなくなる。超えたら強調対象を厳選。(宮野 b01 p.92 / p.79) — 既存「太字/下線は key number のみ」の**定量版**。
- **空白 (余白) を無理に埋めない。** 埋めようと不要なイラスト・背景写真を足すと逆効果。余白は読みやすさの一部。(宮野 b01 p.79, p.136)
- **改行位置は文字数を揃えるより文脈 (文節) で。** 文字数優先で単語の途中や中途半端な位置で切らない。(宮野 b01 p.102 図A)
- **写真上のテキストは文字枠に背景を付けない。** 枠が写真より目立ち情報を隠す。背景なし + 写真と十分コントラストのある文字色で読みを確保。(宮野 b01 p.132)
- **漢字二文字熟語は前後に半角スペース。** 「科学」等が他語と続くと境界が見えにくい。半角アキで語の切れ目を明確化。(宮野 b01 p.102 図B / 本多 b08 と同趣旨)

---

### 追-B. 文章を図に変換する 5 ステップ + 図の型カタログ (宮野 b01)

既存「Phase 2.5 作図力学 5 段階」は**印刷 pipeline** だが、これは**文章→図の論理変換** (別物):

- **イラストレーション 5 ステップ**: ①伝えたい内容を文章化 → ②キーワード抽出 → ③分類 → ④相互関係を表現 → ⑤可視化。いきなり作図せず、キーワードを抽出・分類・関係付けしてから可視化する。(宮野 b01 p.70-73)
- **内容の型に図の型を対応させる**: 流れ → プロセス図 / 論理展開 → 分岐・合体・要素分解 / 分布・拡散 → サテライト・マッピング / 複数の軸 → マトリクス (2軸)・バブルチャート (3軸)。(宮野 b01 p.75-78)
- **長い説明文は体言止めで簡潔化** + **不要な接続句 (「そこで本研究では」「次に」) を画面から削除**。繋ぎ言葉は口頭で言えば足り、スライド上では情報価値ゼロ。(宮野 b01 p.50-52)

---

### 追-C. 聴衆に合わせた背景設計 + 事前調査 (宮野 b01)

- **研究背景は聴衆の専門に応じて「どこから話すか」入口を変える。** 例: 知識ゼロ → 地球規模の問題から / 問題既知 → 高効率発電から / 他分野研究者 → 材料開発から / 専門家 → 当該テーマから。全員に同じ背景は冗長または不足。(宮野 b01 p.32)
- **発表前に会場・聴衆属性 (年齢/業種/専門/人数) を調べ、用語と背景の深さ・準備物を変える。** 機材 (投影等) も事前テスト。聴衆構成が分からなければ入口も用語選択も決められない。(宮野 b01 p.20)
- **背景・目的の self-check (既存 S2 の詳細版)**: 独創性 (先行研究との差はどこか) / 客観性・論理性 (目的と手法がズレていないか) / 正確性 (主張の根拠は十分か)。**「結局、何の問題をどう解くのか」を一言で言えるか**を作図前にクリアする。方法スライドは「その装置スペック・メーカ名を載せる必要があるか」を吟味して情報を絞る。(宮野 b01 p.38-39)

---

### 追-D. 媒体別スタイル — プレゼンでは可、論文では不可 (b02 / b04 / b05)

スライドと論文で**許容スタイルが異なる**点。混同しない:

- **表中の指差し (人差指) / 矢印による強い注目誘導は、プレゼンの投影 (スライド/OHP) でのみ可。** 掲載論文の表では「ふざけた表」と取られうる。(作図力学 b02 §5.3 表5.7) — 既存「視線誘導 (矢印/枠)」の媒体限定条件。
- **コロン後の名詞列挙 (`Processes include: oxidation, hydration, ...`) はプレゼンスライドで機能するが、論文では声に出して不自然なので避け**、`Several processes occur, including ...` の自然文に直す。(b05 §16.11)
- **学会発表は「タイトルだけ読めば内容が分かる」+「データ提示と同時に解釈を示す」**。発表は制限時間があり聴衆が自分で解釈する時間がないので、論文以上に解釈を踏まえた提示にする。タイトルをそのスライドの結論メッセージにする。(b04 Q31)
- **学会発表では Discussion (考察) を薄く (スライド 1-2 枚)。** 入れるなら過去文献との特に重要な比較か研究の limitation に絞る。Result に解釈が織り込まれるため独立 Discussion の比重は小さい。(b04 Q31)
- **導入 (Introduction) を全体時間に対し過大にしない。** 10 分発表に 8 分の導入は不要 (20 分インタビューで 10 分前置きがあったら離れる)。導入が長いほど中身が乏しいと見抜かれる。(b04 Q25 / b05 §14.4)

---

### 追-E. 強調副詞とパラレル構造 (英文発表 b05)

- **重要メッセージは短文で提示し、強調副詞を文頭に置く** — `importantly` / `interestingly` / `remarkably` で (i) 注意を引き (ii) 内容を強調し (iii) 理解を助ける。**ただし deck 全体で 1-2 回まで** (乱発で効果消失)。口頭プレゼンと同じ原理。(b05 §8.8)
- **箇条書き各項目は同じ文法構造 (パラレル) で揃える。** `to evaluate... / to make...` のように全項目を同一構造 (全て不定詞等) に統一。並列構造は複雑な情報の要素を整理して理解を助ける。3 項目以上や後で再言及するなら番号を使う。(b05 §8.5)

---

### 追-F. 図の口頭説明 + 結論の繰り返し (知的プレゼン b06)

OHP 世代の本だが、**図の語り方**として今も有効:

- **図の記号は「名称を前に付けて」呼ぶ。** 「縦軸の記号 σ1 は金属の残留応力であり…」と軸量記号→関係を説明し、以降も「金属の残留応力 σ1」と名称付きで呼ぶ。記号だけ言われても聴衆は意味を記憶していない。論文より説明分量が大幅に増える。(b06 §7.1.2)
- **図と図の連結部に橋渡しの解説を入れる** (唐突に図を羅列しない)。「いま提示した図で理解いただいたことを踏まえ、未解決のこの現象を次の図ではこう解析します」。(b06 §7.1.2) — 既存 outline coherence の図遷移版。
- **講演論文集の図を映すときは図番を読み上げ、同じ図番を投影面にも付ける。** 会場後方でスクリーンが見えない聴講者が手元論文集を参照でき、前の図に戻りたい人にも対応できる。(b06 §7.1.2)
- **結論は最大 3 項目に厳選し、重要項目は本文中 + 結論で計 3 回繰り返す。** 多くの結論を提示しても全部は理解されない。厳選した結論を完全に理解させる方が成功。(b06 §7.1.2) — 既存 take-home / repetition の具体回数。
- **箇条書き項目は一度に全部見せず小出しにする** (pptx は appear animation、OHP は紙で投影光を遮る)。隠された項目への予想・期待で知的緊張が維持され理解度が高まる。(b06 §7.3.1) — 既存 T4 の「段階提示」を原則化。
- **企業報告会では結論 (成果) を冒頭で具体数値とともに述べる。** 「結論から先に述べますと…性能 1.5 倍、価格 80%」。上司・幹部は企業活動への寄与に関心を持つ (成果オンリーの世界)。(b06 §7.1.2)
- **発表原稿の分量目安 m = (0.7〜0.8) × t** (m: 400 字詰め枚数, t: 分; 0.7 でゆっくり, 0.8 でやや早口)。15 分なら 10.5〜12 枚。初心者は盛り込みすぎが逆効果 — 「10 話して 1 分かるより 5 話して 4 分かる」。(b06 §7.1.2) — 既存「20min 4800-5000 字」の枚数版・原則。
- **回答不能の質問には正直に「現段階では…のため分かりません」と答え、喧嘩はしない。** 知ったかぶりは逆に弱点を突かれる。弱点を弱点として認識してもらうのは研究発展の好機。感情的応酬は理性を疑われ研究も知的でないと推断される。(b06 §7.2.3) — 既存 Q&A 復唱に追加すべき態度規範。
- **流暢さより「理解させたい熱意」を前面に。** 弁舌爽やかさより熱意のこもった訥弁の方が説得力を発揮する (考えながら話すことが言葉に力を与える)。(b06 §7.2.1)

---

### 追-G. 耳で聞く日本語の語順 — 一読で誤読されないスライド文 (本多 b08)

スライドの短文・見出しは読み返せないので、**一瞬の誤読が致命的**。本多『日本語の作文技術』の語順四原則を見出しにも適用:

- **見出し・キャプションは「どの語がどこにかかるか」が一読で取れる語順に。** 修飾語を被修飾語に直結し、長い修飾語・節を先に置く。二義になるなら語順を変えるかカギカッコで区切る (例: 「野蛮な文明の敵」→「文明の野蛮な敵」)。(本多 b08 第2-3章)
- **否定文は限定の「ハ」を補って肯定/否定を一義に。** 「完全に回復しなかった」(=全然回復していない?) は誤読を招く。「完全には…ない」「日本のようには…ない」と否定動詞とセットにする。数量の全部/一部を示す否定は特に注意。(本多 b08 第6章)
- **無色の接続助詞「ガ」で文を続けず、句点で切る。** 「ガ」は聴衆に逆接を予期させ思考を一瞬乱す。逆接でない「ガ」を見つけたら切り、一項目一思想にする。(本多 b08 第6章)
- **カナ連続を避け、漢字・カタカナ・半角アキで語のまとまりを視覚化。** 同種の字が続くと一字ずつ拾い読みになり、瞬時に読ませるスライドでは致命的。専門語を半分カナにしない (読みはカッコで添える)。(本多 b08 第5章) — 追-A の半角スペース規則の上位原則。

関連 tool: `presentation_check_misuse_japanese` / `presentation_check_notation_variants` / `presentation_check_kanji_ratio` (台本・スライド原稿の和文 lint)。

---

## 🏆 CEFC/Compumag oral playbook — 自研究室の実録から (2026-06)

過去の自研究室 CEFC/Compumag **口頭発表**（CEFC2020 Yano: MMM-MSC / Compumag2023
Tanimoto: CLN-TEAM28 / CEFC2024 Sugahara: FP-CLN）から抽出した、この会場で通る型。
各 talk の field note は `talk_feedback.py`。**ポスターは別物**（`radia-poster` 参照）。

### talk arc（計算電磁気・手法提案の口頭）
1. **Title** ＋ 原案者クレジット（例「CLN は Kameari 2016」）
2. **Background / 位置づけ** — 聴衆が**実際に使う固有ツール名**に対して置く
   （OPERA-3D TOSCA / CST / COMSOL / Radia / FEM / BEM）。"de-facto standard は X"
3. **Gap** — 従来法の何が問題かを**専用スライド＋具体例**で（"Why MMM is not good?" /
   "FEM は air mesh が要る" / "BEM は時間がかかる" / "Radia の parallelepiped mesh 誤差"）
4. **Formulation** — 定式化（記号定義＋物理意味）
5. **Validation（最重要・厚く）** — 次を複数:
   - **TEAM Workshop Problem**（例 28＝磁気浮上）＝この分野の信頼の通貨
   - **商用ソフト一致**（CST / OPERA-3D）／**解析解**（"no need for FEM to build basis"）
   - **収束 study**（mesh 細分 N=1..10 を reference に対してプロット）
   - 例題を**複数・難度順**（T字→C字→quadrupole ／ 1D→1D→2D）
6. **計算コスト比較** — 劇的な1数字を前面に: "OPERA-3D 72 min → 2 min" /
   "commercial 40 min → 5 s" / "204 GB → 8.6 GB"
7. **Discussion / Next step** → 8. **Summary**（1行 take-home）

### この会場の作法（実録の共通項）
- **固有ツール名で位置づける**（汎用の"従来法"でなく OPERA-3D/CST/Radia/FEM/BEM）。
- **検証は必ず1段持つ（強い順 ladder）= TEAM problem ▸ 商用ソフト一致(CST/OPERA-3D) ▸
  解析解 ▸ 実測 ▸ 二定式化の相互照合(A-φ vs T-Ω) ▸ 手法間比較**。実録(16本)では
  **下2段（二定式化照合・手法間比較）が最頻**で、TEAM は自研究室でも稀（差別化要素で
  あって baseline ではない）。無い段は「今後」と**正直に**。(→ § Corpus update)
- **収束・mesh 依存性の plot を厚く**（この聴衆は数値手法の収束を必ず見る）。
- **劇的なコスト数字を1つ**前面（speedup / memory）。
- **先行手法・原案者をクレジット**（PVL/POD/Kameari…）。喧嘩を売らず差分を明示。
- **ノートは SSML/間で pacing**（`<break time="1s"/>`、CEFC2020 実録）→ 時間厳守。
  `presentation_speaking_pace_estimate` と併用。
- backup スライド（formulation 詳細）は本編の後ろ。

> 使い方: 新 talk を作る前に **`presentation_qa_from_history`** ＋ 本 playbook と照合。
> 発表後は実 Q&A を `talk_feedback.py` に `recorded` で還流。

---

### Corpus update (+12 orals, 2026-06)

Mined 12 additional decks (CEFC2016/2022/2026, Compumag2017/2025, IGTE2023×3
/2024/2026) and cross-checked against the existing playbook (§"talk arc",
§"この会場の作法"). This update QUANTIFIES the validation spine over the broader
corpus, adds conventions the current playbook does not yet state, and pins
down the **venue-genre split** (a CEFC/Compumag *oral* vs an IGTE *research
seminar* are different animals).

#### 1. Validation strategy across the 12 (the headline finding)

The existing playbook calls TEAM-problem validation "事実上必須の通貨". The
broader corpus shows that is the *aspiration*, **not** what most of these
decks actually did. Counting the 8 decks that are genuine results talks
(excluding new02 figure-bank, new04/new07 non-technical seminars, new05
title-only-stub):

- **TEAM Workshop problem: 1/8** (new06_igte2023, TEAM #28). The currency is
  real but **rare** even in this lab's own output — treat "we have a TEAM
  number" as a *differentiator*, not a baseline.
- **Commercial-software cross-check (CST/OPERA/COMSOL): 0/12.** None of the
  new decks cross-check against a commercial code. (The 3 *original*-playbook
  decks did — so the lab's strongest validation precedents pre-date this
  batch.)
- **Analytic closed-form solution: 2/12** (new01_cefc2016 sphere→cylinder→
  ellipsoid→elliptic-cyl→periodic; new02_compumag2017 the figure twin of it).
- **Measurement / experiment as ground truth: 3/12** (new09 Toyo-Univ B-H
  loops; new11 F.W.Bell-8030 rig vs in-house model; new12 measured-loop
  reconstruction). Measurement is **as common as analytic** here.
- **Formulation-vs-formulation mutual cross-validation: 2/12** (new10
  A-φ vs T-Ω, 0.21%/0.11%; new03 CLN Method-1/2/3 vs the underlying
  2-D multi-slice FEM). **This is the single most reusable substitute when
  no TEAM/commercial number exists** — solve the same problem two
  independent ways and report the small mutual error.
- **Internal convergence / sensitivity study (mesh-size, mesh-type,
  formulation, hysteron-count): ≥4/12** (new06, new08, new09, new12) — this
  matches the existing playbook's "収束 study を厚く" and remains universal.
- **NO validation at all: 3/12** (new04, new07 self-intro/travelogue; new02
  figure-only) — see §3 genre warning.

**Refinement to the playbook's validation rule:** the priority ladder it
should teach is **TEAM ▸ commercial-agreement ▸ analytic closed-form ▸
measurement ▸ two-formulation mutual cross-validation ▸ method-vs-method**.
A new oral needs at least ONE rung at or above "two-formulation
cross-validation"; the corpus proves the bottom two rungs (cross-validation,
method-vs-method) are how the lab most often clears the bar in practice.

#### 2. NEW conventions / refinements (not in the current playbook)

- **Headline number is frequently ABSENT — and that's the #1 fixable gap.**
  Of the results talks, only new10 carries a quantified headline (0.21%/0.11%
  cross-validation error). new03 has a *dedicated* "comparison of
  computational time" slide but the actual speedup figure did not survive;
  new08 (CMA-ES) and new12 (energy-Stop) have NO dramatic number at all. The
  playbook's "劇的な1数字を前面" rule is **the most violated** — add it to the
  pre-talk checklist as a hard gate, not a nicety. A cost slide is not a cost
  *number*.
- **Analytic-validation LADDER ordered by increasing difficulty.**
  new01_cefc2016 is the template: isotropic sphere → uniaxial cylinder →
  biaxial ellipsoid → elliptic cylinder → periodic. One geometry per
  setup+result slide-pair, blank/formula divider between. The *ordering*
  (isotropic→anisotropic→periodic) signals to the reviewer you stress-tested
  the hard cases, not one easy one. New rule: **when validating a BC/open-
  boundary/integral method, use a difficulty-ordered ladder of closed-form
  geometries, one auditable rung per slide-pair.**
- **Numbered competing-method spine + fixed color legend.** new03 (CLN
  Method-1 green / Method-2 blue / Method-3 red, reused on every results
  slide) and new08 (Warm-Start CMA-ES vs GA vs plain CMA-ES) both make the
  talk a *controlled comparison*, not a single-method pitch. The fixed legend
  across all results slides lets the audience read accuracy AND cost at a
  glance. New rule: **if you have method variants, number them, color-lock
  them once, and reuse the legend verbatim — the structure does the
  persuading.**
- **Earn trust BEFORE claiming speed (mechanism-before-cost).** new03 shows
  the reduced CLN circuit reproduces the FEM B-field and eddy-current J
  *from the basis vectors* (S13/S15) before the climactic time slide, so it
  reads "fast AND faithful" not "fast but lossy". State the explicit DOF
  count (new03: M·Ns) so reviewers can size the problem. New rule:
  **reduced-order / acceleration talks must show field fidelity first and
  put problem-size DOF bookkeeping on its own slide.**
- **Recurring agenda slide as a progress divider.** new08 reuses one
  identical agenda slide (S2/S6/S12/S22) before each act; new03 uses a
  *duplicate* TOC (S7 + S16) to split a formulation half from a results half.
  This is a concrete realization of the existing "navigation" idea — promote
  it to an explicit convention for methodology-heavy talks.
- **Two-case generality check ("Case A, then Case B as additional
  verification").** new08 answers the reviewer's "does it only work on your
  one example?" *before it is asked* by framing Case B (iron loss / V-shape)
  as generality evidence. New rule for optimization/design talks:
  **always carry a second, structurally-different case explicitly labeled as
  generality verification.**
- **Problem → trap → fix triplet on one tight slide group.** new12 shows the
  non-monotone shape-function *trap* (S9) immediately answered by the
  ridge-penalty *fix* (S10). Showing the failure you hit and its cure is more
  convincing to an expert audience than a polished-only method.
- **Credit-the-data-source on the slide where data appears** (new09 stamps
  東洋大/Toyo-Univ on every B-H slide), distinct from credit-the-method
  (Egger-2025 in new12's title, ESRF/Radia "we forked it" in new10). Both
  forms appear; the data-credit form is new vs the current playbook's
  method-credit-only guidance.
- **For experimental EM talks: validate the apparatus before the physics,
  then close the loop measure→model→design.** new11 puts a dial-gauge
  parallelism slide and a named sensor (F.W.Bell-8030) + specs slide *before*
  any field result, and ends by feeding the validated model's gradient
  forward into a coil design. New rule: **name+spec the instrument, show the
  calibration result first, and end on the forward use, not the agreement
  plot.**
- **Remote/scripted-delivery hygiene.** new10 was delivered remotely (flight
  disruption) with verbatim SSML-paced notes and tool-primer slides
  (CoefficientFunction vs GridFunction) so a mixed audience could follow —
  reinforces the existing SSML note but adds: **front-load a 1–2 slide primer
  on the enabling abstraction**, and keep cost-vs-DoF + extra background as
  *titled backup slides* after Thank-you (new10 does exactly this).

#### 3. Venue differences (CEFC vs Compumag vs IGTE seminar) — IMPORTANT

The current playbook is written for the *oral* and assumes every deck is one.
The corpus shows **the filename venue tag is NOT a reliable indicator of
genre** — several "igte2023/2024" and "cefc2026" files are research seminars
or progress decks, not podium orals.

- **CEFC / Compumag oral (4-min digest register):** short, results-spine,
  validation + headline-number mandatory. new01/new03/new08/new10 fit this.
  This is the genre the existing playbook targets — keep applying it as-is.
- **IGTE "Research Seminar" (TU Graz) register:** exploratory, work-in-
  progress, broader and more narrative. new05 (Kato, coupled-IH CLN), new06
  (Fujita, CLN gauging) are *technical* seminars — they DO carry formulation
  + a benchmark (new06 names TEAM #28 on slide 3) but in a looser,
  one-message-per-slide spoken-outline form with dedicated citation slides up
  front and **backup/appendix slides after the Summary** (new06: Cubit→Gmsh,
  CLN-history). Harvest their *credit-early + name-the-benchmark-early*
  discipline.
- **IGTE/seminar NON-technical (the trap):** new04 (Keiko Sugahara
  self-intro), new07 (Sugahara "What I learned in Austria" travelogue) are
  tagged with conference names but have **no validation, no headline number,
  no positioning** — they are community/relationship talks. **Do NOT reuse
  these as oral templates.** Harvest only their narrative/credit conventions
  (bookend motif, named-collaborator slides, explicit closing call-to-action
  / "HOMEWORK" slide). Before treating any mined deck as an oral template,
  verify it has the validation/positioning/headline-number spine — these two
  prove a venue tag can lie.
- **Working/backup decks underpinning an oral:** new09 (ヒスフィッティング
  curve-fit pipeline), new12 (energy-Stop status), new13 (its sibling) are
  Japanese data-conditioning/status decks that *feed* a polished English
  oral. Convention: **keep the dense fit machinery (17th-order polynomial,
  rational envelope) in a separate backup deck; present only the
  measure→fit→generate storyline in the oral.**

#### 4. One-line additions to drop into the existing checklist

- Add headline-number as a **hard gate**: "a cost *slide* is not a cost
  *number* — if no quantified speedup/accuracy/iteration figure exists, the
  oral is not ready." (Most-violated rule in this corpus.)
- Add the validation **substitute ladder** explicitly (TEAM ▸ commercial ▸
  analytic ▸ measurement ▸ two-formulation cross-validation ▸ method-vs-
  method) so authors without a TEAM number know the next-best rungs the lab
  actually uses.
- Add a **genre check** to "使い方": confirm the reference/template deck is an
  *oral* (has the results spine), not a seminar/self-intro/figure-bank.

---

## 📡 Talk feedback loop — CEFC/Compumag field notes (2026-06)

**最強の長期 training は「実際の talk で何を聞かれたか」を貯めて還流すること。**
`meta/bug_patterns.py` の presentation 版として `talk_feedback.py` に会場別の
Q&A・座長/客の反応・教訓を蓄積する catalog を持つ。

### Before a talk（準備）
- ★ **`presentation_qa_from_history(topic=...)`** を必ず実行 → 過去の同種 talk で
  実際に出た（or 予測した）質問を一覧。各質問に **30 秒で答えられるか**を素読み確認。
- `presentation_talk_feedback_lookup(venue="CEFC")` で会場特有の作法・教訓を確認。
- 予測 Q&A を `status="anticipated"` で 1 エントリ FIELD_NOTES に追加しておく。

### After a talk（24 時間以内に）
- その `anticipated` エントリを **`status="recorded"`** に更新し、**実際に出た
  質問**・座長/客の反応・うまくいった点・改善点を埋める。これが最も濃い venue
  signal — 記憶が新しいうちに。教訓は `rules` に蒸留し、次回 `qa_from_history`
  から自動的に効く。

### CEM (CEFC/Compumag) で必ず来る質問（seed = TM-O1 より）
1. FEM と比べて精度/コストは?（この聴衆は全てを FEM と比較する → 早めに定量で）
2. ACA/H 行列近似で精度・反復回数は劣化しない?
3. N に対するメモリ/時間 scaling は? sub-cubic?（**過大主張しない**・正直な caveat）
4. TEAM 問題 or 実測で検証した?（**必ず先読み**・正直に・できれば次回までに TEAM 追加）
5. FMM-BEM (ngsolve.bem) でなく MMM+HACApK の理由は?
6. 3D か? 非線形鉄は扱えるか?（明示的に早く言う）
7. コードは公開? C/Fortran?

tool: `presentation_qa_from_history` / `presentation_talk_feedback_lookup` /
`presentation_talk_feedback_stats`（catalog: `presentation/talk_feedback.py`）。

---

## 📤 PDF → PPTX 変換: `pdf2ppt-pdfgear` skill (PDFgear 上級モード) を使う (2026-06-14)

**POLICY**: 既存の PDF スライド (beamer 出力等) を PPTX に変換する場合は、
**`pdf2ppt-pdfgear` skill** を使う。 これは PDFgear desktop app を
**UIAutomation 経由で自動駆動** し、 **上級モード (Advanced Mode) を ON** に
して editable な PPTX を出力する。 手動で GUI クリックしない。

**Canonical home** (実装ファイルの真の置き場所):
`packages/radia-mcp/skills/pdf2ppt-pdfgear/{SKILL.md, pdf2ppt_pdfgear.py}`

**Discovery path** (Claude が skill 発見する場所、 canonical home への symlink):
`C:/Users/Administrator/.claude/skills/pdf2ppt-pdfgear/` (Windows symbolic link、
S:\ NAS drive 越しでも動作 — junction は対応不可)

**Setup** (一度きり):
```powershell
# canonical home に実装が無ければコピー、 ある場合は不要
Remove-Item "C:\Users\Administrator\.claude\skills\pdf2ppt-pdfgear" -Force -ErrorAction SilentlyContinue
New-Item -ItemType SymbolicLink \
  -Path   "C:\Users\Administrator\.claude\skills\pdf2ppt-pdfgear" \
  -Target "S:\Radia\01_GitHub\packages\radia-mcp\skills\pdf2ppt-pdfgear"
```
junction (`-ItemType Junction`) は S:\ (\\192.168.11.100\work) で「再解析ポイント
バッファー無効」エラー。 SymbolicLink を使うこと (2026-06-14 切替済)。

**使い方** (discovery path 経由 — symlink で canonical home を解決):
```powershell
python "C:/Users/Administrator/.claude/skills/pdf2ppt-pdfgear/pdf2ppt_pdfgear.py" \
       "input.pdf" -o "S:/path/out.pptx"
```
canonical home を直接呼んでも OK:
```powershell
python "S:/Radia/01_GitHub/packages/radia-mcp/skills/pdf2ppt-pdfgear/pdf2ppt_pdfgear.py" \
       "input.pdf" -o "S:/path/out.pptx"
```
- `-o` なし: PDFgear のデフォルト出力先 (`%USERPROFILE%/OneDrive/PDFgear/<stem> conv.pptx`)
- 上級モードは **skill の default** (`--no-advanced` で OFF; 通常使わない)

**上級モード ON で何が起きるか** (実測):
- 14 ページ beamer deck → **315 text box + 0 raster image** の editable PPTX
- タイトル / 本文 / 数式 / 図形がすべて text shape + vector shape として保持
- 日本語・LaTeX 記号も editable (OCR ではなく PDF text stream を解析)

**上級モード OFF (`--no-advanced`) では**:
- 1 ページ = 1 巨大 PICTURE shape (PyMuPDF レンダリングと同等)
- 編集不能、 ファイル肥大 → **通常は使わない**

**重要なカベアト**:
- **クラウドアップロード**: 変換は `apiw.pdfgear.com` に PDF を POST する
  → **未公開研究データには絶対使わない**。 機密はローカル PyMuPDF → python-pptx パス
- **sign-in**: 初回は PDFgear GUI で手動 sign-in が必要 (skill は無人実行のため)
- **数式 fidelity**: 上級モード ON でも beamer 数式は **glyph 単位で再 layout** される
  ので、 editable だが見た目を厳密維持したい場合は **`.tex` 原稿が真の editable master**
  (編集→ recompile が正攻法)

**典型的な用途**:
- 共著者が PowerPoint でコメント書き込みたいとき
- 学会が PPTX 提出を要求しているが、 原稿は beamer のとき
- 過去 PDF の slide 1 枚を新 deck に取り込みたいとき

**他選択肢との比較**:
- **`pdf2ppt-pdfgear` skill** (本 policy): クラウド OK な公開済み資料、 editable 必要
- **PyMuPDF → python-pptx (offline)**: 機密データ、 image-only でよい、 deterministic
- **PowerPoint で最初から作る**: 継続的に PPTX 必要なら beamer を捨てる
- **LibreOffice / Adobe Acrobat の標準 PDF→PPT**: 画像化されるため不推奨
- **`radia-doc-convert.doc_convert_pptx_to_pdf` は逆方向** (PPTX → PDF) なので別物

**関連 tool / skill**:
- `pdf2ppt-pdfgear` skill (本体)
- skill 内 `SKILL.md` に reverse-engineered detail (`pdfconverter.exe app1 PDFToPPT
  <conv_docs.json>` 経由起動、 UIA TogglePattern で 上級モード ON、 InvokePattern で
  変換ボタン発火)
- 出力 PPTX を `presentation_check_*` の lint tool 群でチェック (font / bullet / 数式)

---

## 関連 MCP

- `grant-writing` - 申請書の作文技術 (科研費 / JSPS / KDDI)
- `paper-writing` - journal 論文の作文技術 (IMRAD + reviewer)

三者は共通診断 tool (overfull hbox / sentence length / weak expressions) を持ち、閾値を用途別に調整している。
