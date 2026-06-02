# Presentation Skill — 学会発表スライドを通すための作文技術

研究発表スライドの作文ガイド。国内 (IEEJ 系 / 応物 / 機械学会) + 国際 (IEEE conference / APS March Meeting / Compumag 等) + 社内セミナー共通の設計原則。

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

## 関連 MCP

- `grant-writing` - 申請書の作文技術 (科研費 / JSPS / KDDI)
- `paper-writing` - journal 論文の作文技術 (IMRAD + reviewer)

三者は共通診断 tool (overfull hbox / sentence length / weak expressions) を持ち、閾値を用途別に調整している。
