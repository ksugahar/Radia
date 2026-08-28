# Paper Writing Skill — Journal 論文を通すための作文技術

Journal paper 執筆ガイド。CAE-AI Lab の対象は **電気学会 (IEEJ) / IEEE / 加速器系学会**に限定し、IMRAD 構造・contribution 明示・reviewer 対応の実戦ノウハウを統合する。

---

## CAE-AI Lab 投稿先ポリシー

- 正規ターゲットは `ieej`、`ieee`、`accelerator` の3系統だけとする。
- 対象外学会への journal-fit 最適化や投稿推薦は行わない。投稿先が未定なら、まず3系統の中から具体的な venue を決める。
- 電気学会では、和文本文と英文題目・abstract・keywordを別々にレビューする。
- IEEEでは英文を投稿原稿の正本とする。和文の思考整理稿がある場合も、英文との平均点ではなく独立に診断する。
- 加速器系では JACoW / PRAB / PASJ 等の具体的な公式テンプレートを先に選ぶ。ページ数やkeywordを「加速器系共通」と決め打ちしない。
- 投稿年度の公式テンプレートを source of truth とし、ページ数、abstract、keyword、reference、figure規則をその都度確認する。
- `paper_writing_target_venue_policy(target_venue)` で分類し、対象外をrejectする。`paper_writing_em_submission_gate(target_venue=...)` にも同じ投稿先を渡す。
- 和文と英文の可読性は単位も失敗様式も異なるため、スコアを平均しない。適用される言語のうち悪い側を総合ゲートとする。

Radiaの成果は、電気学会で日本語の理屈を磨き、IEEEで国際的な電磁界論文として確立し、加速器系で電磁石・運転・ビーム調整への効果を示す循環を基本とする。

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
Figure caption 単独で図の対象が理解でき、本文側でも主張が回収されているか。
- 目標: caption は plotted quantity・条件・series を同定し、主張・数値・解釈は本文にも必ず書く
- 推奨: "Fig. 3. <主張の動詞句>: <セットアップ>, <キー数値>" としつつ、同じ key result を本文の図参照段落にも置く
- caption が長くなる場合は本文を優先し、比較ロジック・物理解釈・結論を本文へ移して caption を短くする

### T5. Citation density → ✅ `paper_writing_related_work_density` (v0.13.0)
Introduction 内の引用密度・自己引用比率・年度分布を診断。
- 目安: Intro に 20-40 件 (研究動向 20, 直接比較対象 5-10)
- 自己引用比率 < 20%、直近 5 年分 40% 以上が目安

### T6. 制約検証 → ✅ `paper_writing_validate_pdf_pages` / `paper_writing_validate_abstract_length` / `paper_writing_check_overfull_hbox` / `paper_writing_check_pdf_edge_overflow`
IEEJ / IEEE / 加速器系 venue 固有の制約検証:
- ページ数は投稿年度の公式テンプレートで確認し、過去年度や別venueの値を流用しない
- Abstract は IEEJ / IEEE / JACoW / PRAB / PASJ の選択後に公式上限を適用する
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

既知の用語ゆれを直すときは `paper_writing_normalize_terminology` または
`paper_writing_normalize_terminology_file` を使う。既定 rule は
`cube=>立方体` で、日本語文脈に接した Latin 語だけを置換するため、
英語 caption・bibliography・`\includegraphics{fig_cube...}` は保持される。
ファイルに適用する場合はまず `dry_run=True` で置換例を確認し、問題なければ
`dry_run=False` にする。

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

#### Bilingual digest page-limit policy (EN venue-limit strict, JA synced)

IEEE conference / IEEJ研究会 / JACoW加速器論文などを **英語版 + 日本語版の
twin** として保守する場合、最初に必ず投稿先・年度・document type ごとの
page limit を確認する。英語版は投稿物なので、確認した venue limit
(例: 1 page / 2 pages / 4 pages) を厳守する。過去の別学会・別年度の
「1-page digest」慣習を無確認で流用しない。

一方、**日本語版は英語版と内容・数値・式・図・引用を完全同期**させるが、
英語版の page limit は課さない。和文は共著者レビュー・思考整理用の読みやすさを
優先し、英語版の page-limit compression を和文に機械的に移植して読みにくく
しない。英語版を削ったら日本語版も同じ主張範囲へ同期するが、和文の文章量
そのものを英語版に合わせて削る必要はない。

編集時は EN/JA を同じ change set で更新し、英語版は page-count を検証、
日本語版は bilingual sync と組版エラーを検証する。EN が venue limit を超える場合は
「冗長部削除 → 小さな `\vspace{-...}` 調整 → 低優先文の削除 → 英文圧縮」の
順で対応し、英文圧縮で不明瞭になるなら文を削る。

#### Reviewer Q&A driven revision policy

Digest / paper を修正するときは、まず査読者が突っ込みそうな点を
**Q&A** として書き出す。各項目は (1) reviewer question, (2) answer,
(3) manuscript action, (4) space が足りない場合に削る候補、の順で整理する。
この段階では文章が一時的に増えてよい。Q&A で論点を潰してから本文へ反映し、
page limit 超過は過剰な英文圧縮ではなく、冗長部削除または優先度の低い内容の
削除で解決する。

本文へ反映するときは、読みにくい compressed prose を作らない。関係を説明する
ために文中で `=` を助詞代わりに使わず、等号は数式・表・明示的な定義式に限る。
意味が明確でない括弧補足や、括弧内に文を押し込む書き方も避ける。括弧は
短い例、単位、略語初出、図表・式参照、citation に限定し、説明が必要な内容は
本文の通常の文として書く。

#### Fixed-iteration reporting and fair solver comparison

線形・非線形反復を固定回数で打ち切る実験は有効だが、固定した事実を隠さない。
方法節に、更新系列と回数、内側反復回数、緩和係数、前処理、残差または観測物理量
による採用条件、過去最良状態を保存する場合の選択規則を書く。固定回数は
「収束回数」ではなく **固定計算予算** と呼び、別途、採用解が所要精度を満たす
ことを示す。

残差ゲートを満たす履歴から観測量変化が最小の状態を選ぶ方法は、反復安定性や
メッシュ間整合性を示す根拠にはなる。ただし、解析解、独立参照解、測定値の
いずれもない場合は絶対精度の証明とは呼ばない。「内部状態より観測量の反復感度が
低い」と「観測量が真値に一致する」を区別する。

反復回数を直接比較してよいのは、更新の定義、初期値、残差norm、tolerance、
停止条件、最大回数が同じ場合だけである。一方が残差収束、他方が固定予算または
観測量ベースの最良解選択なら、同じ表の反復回数列で優劣を示さない。物理量、精度、
実測時間、実測メモリを比較し、異なる反復プロトコルは方法節または補足に分けて書く。

全問題規模で外側反復回数を固定したスケーリング測定では、その目的を
「問題規模による非線形反復回数の変動を除き、行列構築・分解・行列ベクトル積の
規模依存を評価するため」と明記する。この測定の近似指数は固定予算内のworkloadの
スケーリングであり、非線形収束まで含む計算量の指数とは主張しない。

#### platex+dvipdfmx 原稿の図は PDF 一択 — PNG は無警告で白紙になる

出典: SA-26-078（八戸, 2026-07-17）。IEEJ 系 platex→dvipdfmx テンプレートで
matplotlib の **PNG を \\includegraphics すると、図が本文にレイアウトだけ確保
されて中身が白紙**になる事故が起きた（TeX Live 2026: graphicx/dvipdfmx.def は
PNG を `em: graph` special で埋め込むが、現行 dvipdfmx バイナリがこれを
"Unparsed material ... ignored" と黙って捨てる）。**コンパイルはエラーゼロで
通る**うえ、ラボ標準の compile.ps1 は dvipdfmx 出力を Out-Null に捨てるため
警告も見えない。共著者が「無駄な空白が多い」と指摘して初めて発覚した。

ルール:
1. **platex+dvipdfmx 原稿の図は PDF（matplotlib はベクトル PDF, pdf.fonttype=42）
   のみ使う。** PNG しかない図（スクリーンショット等）は PIL で PDF に包む。
2. **図を再生成したら .xbb を必ず更新**（extractbb を明示実行）。stale .xbb が
   残ると旧図の縦横比の箱に新図が押し込まれ、潰れ・過剰余白になる。
3. **投稿前に「図が PDF に埋まっているか」を機械検証する**: pymupdf で
   `get_page_images`（ラスタ）+ `get_page_xobjects`（ベクトル Form）を数え、
   \\includegraphics の数と突き合わせる。ページを画像レンダリングして図領域が
   白紙でないか目視する（図の**ソースファイル**を眺めるだけでは検出できない —
   PDF 内のレンダリングを見ること）。
4. dvipdfmx の stdout/stderr を握りつぶさない（少なくとも warning を grep）。

#### 誠実な定量報告 — 打ち切り値・in-sample 値・「データ結果」と「一般主張」

出典: SA-26-078（能動学習 B 入力ストップ同定, 2026-07-17）仕上げトレース。
数値そのものは正しくても、**何を意味する数値か**を偽ると reviewer 信頼を失う。

**打ち切り(censored)値・フォールバック値を「到達値」として報告しない**:
探索・反復・能動学習でしきい値に到達しなければ「試した最大予算」を返す実装は
多い（例: `n = next((i for i,v in enumerate(curve) if v<=thr), N_max)`）。この
`N_max` は「N_max で到達した」ではなく「N_max でも未到達→フォールバック」を
意味する打ち切り値である。これを「ランダムは N_max 回を要する」と書くと未到達を
到達と偽る。正しくは「N_max でも当該水準に届かない」。固定反復を「収束回数」と
呼ばない規律（上節）と同じ精神。実装の返り値がフォールバックか到達かを必ず
確認してから文章化し、可能なら「到達 vs 未到達」を非対称に書く（片方は達成、
片方は予算内未達）。これは弱主張化ではなく、むしろ正しく書くと主張が強くなる
ことが多い（相手が予算内で追いつけない、と言えるため）。

**in-sample の自己再現値を精度の見出しにしない**: 学習データ自身を再現する値
（例: 剥ぎ取りで構成した参照モデルが、その構成に使ったループを再現する誤差）は
近似的に自明であり、汎化性能ではない。見出し数値には out-of-sample / 未学習
入力の一致度を使う。in-sample 値を出す場合は「同定に用いた〜をよく再現する」と
質的に述べるに留め、tight な数値（0.x%）で generalization を示唆しない。
Sugahara は「0.6% はいいすぎ」と、この in-sample 値の見出し掲載を却下した。

**「このデータの結果」と「一般的主張」を配置で区別する**（[[feedback-result-claim-vs-general-claim]]）:
同じ数値でも、特定条件下の測定結果として書けば妥当だが、手法一般の性能として
書くと過大主張になる。数値の意味は置かれた文脈で決まる。abstract・結論・図
キャプションに数値を置くときは「この実験の結果」か「手法の一般性能」かを意識し、
後者に化ける書き方（"the method achieves X%"）を避け、"in this setup, X%" の
ように限定する。

**節をまたぐ同一数値は同一量を指すこと**: ある数値（例: "14"）が複数節に現れる
とき、それぞれが同じ量を指すか確認する。§A「能動が 14 で到達」、§B「ランダムが
14 でも未到達」のように同じ "14" が別物を指すと読者は混乱する。模型規模・判定
基準・測定対象が異なるなら明示して区別する（`title_abstract_conclusion_triangle`
の数値整合を節間にも適用する感覚）。

**比率は効果量を伝えるが、同時に情報を捨てる**: 「23分の1」「145倍」のような
比率は効果の大きさを直感的に示すのに有効であり、使用自体を禁止しない。一方、比率
だけでは分子・分母の絶対値、基準量、単位、評価条件が失われる。原則として元の2値と
共通条件を先に示し、比率は補助的な効果量として用いる。特に、異なるメッシュ分割・
試験群・集約規則で生じた最大値同士を割って一つの改善倍率にしない。また、基準値が
ほぼゼロのときは比率が過大に見えるため、絶対差または絶対値を主表示とする。抄録・
結論・図中の強調値を何でも比率にせず、読者が結果の規模を復元できる情報を残す。
同一条件で「10.922\%から0.468\%へ低下（約23分の1）」と書くことは許容されるが、
条件の異なる「0.159\%（$n=2$）と0.00110\%（$n=3$）の比は約1/145」は避ける。
`paper_writing_check_misleading_ratio_claims` は比率そのものではなく、元の2値がない
比率、条件混在、1文への比率の詰め込みを警告する。

**正則化・近似・高速化の役割を混同しない**: 物理量に罰則または制約を課す
Tikhonov 正則化と、低ランク近似・モード打切りによる数値的安定化は、結果として
いずれも高周波成分を抑えても目的が異なる。主たる正則化を物理量制約が担い、打切りを
高速化に使うなら、その順序を序論・定式化・結果で一貫して明示する。「打切り後も解が
変わらない」と主張する場合は、十分高いランクの**同じ正則化問題**を参照解とし、採用
ランクとの解差のノルム、判定許容値、採用ランクを報告する。特異値分布の一致だけでは、
正則化解の一致を検証したことにはならない。

**定式化を変えたら数値・図・表を再計算する**: 旧手法の出力ファイルを改名して新しい
定式化の結果として使わない。目的関数、制約、正則化項、行列分解、評価領域のいずれかを
変えた場合は、その原稿で使う実装から全数値と図表を再生成し、本文の式・コード・保存
データ・図キャプションが同じ問題を指すことを確認する。

**評価領域は成功領域だけに切り詰めない**: 空間分布や局所近似を評価する図は、関心領域
だけでなく、その外側で性能がどのように崩れるかも読者が判断できる範囲を示す。二次元・
三次元量を断面で示す場合は、代表する複数の位置でプロファイルを比較し、測定値を追加する
予定なら同じ座標・同じ軸上に重ねられる構成にする。

#### 用語・計算機名・製品名の規律（再現性に寄与する語だけを出す）

**用語は標準語で一貫させる**: 自作ジャーゴンを避け、標準語を統一して使う。計算結果を
保持して以後の要求に再利用する機構は「台帳(ledger)」のような造語ではなく **キャッシュ
/ メモ化 (cache / memoization)** と呼ぶ（実装コメントが memo と呼ぶなら論文も揃える）。
同じ機構を別名で呼び分けない——一方を「並進キャッシュ」他方を「ブロック台帳」とする
ような不整合を作らず、両方 cache/memoization に統一する。標準語がある概念に新語を
充てない。

**計算機名は一度だけ**: 使用した計算機（例: mdx）は方法節と謝辞で一度述べれば十分。
以後は「mdx計測」のように内部ホスト名を測定の枕詞にしない。読者に意味のない内部機器名
を繰り返さず、機能で「計算時間の評価 / 実測」と書く。表・図のキャプションにも機器名を
貼り重ねない。内部機器名それ自体は再現性を上げないので、プラットフォーム名は一度の
明記と謝辞に留める。

**製品名は再現性で判断する**: 製品名（商用ツール名を含む）を出すのは、その名を示すこと
で **再現性が上がるとき** だけである。汎用的で代替可能な作業——例えば球のメッシュ生成は
どのメッシャでもよい——では特定製品名を出さない。品質的にその製品でなければ再現できない
場合に限り、理由とともに名を挙げる。手法の中核が依存するフレームワーク（例: NGSolve の
要素・.vol）は再現性を上げるので明記してよい。

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
- ORCID / 著者 affiliation・郵便住所を、投稿時点の大学・研究機関の公式ページで確認

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

## 🚫 NG パターン 13 (journal / digest 特有)

1. **Abstract に数式 / 引用 / 専門略語** — display math (`\begin{equation}`, `\[...\]`),
   inline math (`$...$`), `\cite{}` / `[1]`, FEM/BEM/MCP/LLM などの
   domain-specific acronym はすべて abstract から除外。
   理由: (a) 検索エンジン (IEEE Xplore / Scopus / Google Scholar) は abstract
   をプレーンテキストとして index → math symbols が garble、(b) self-contained
   原則 (引用や略語表なしで主張が独立して読めるべき)。**例外なし** — 創設論文・先行
   研究・共著者論文を引用したくても abstract には名前で言及し `\cite{}` は
   本文に移す ([[project-esim-hollaus-coauthor-bib-2026-05-29]] 参照: 共著者
   への礼儀として abstract に cite したくなる衝動を抑える)。プレーン数値
   ("18%", "50 kHz") は OK、LaTeX math symbols (`$Z_s$`, `$O(h^3)$`,
   `$\mathcal{E}$`) は NG。略語は本文初出で `Finite Element Method (FEM)`
   のように展開する。`paper_writing_check_abstract_no_math_no_citation`
   で機械検証 (display + cite = fail, inline math + acronym = warning)。
2. **Contribution が discussion 調** ("We discuss..." で終わる) — 動詞で宣言
3. **Abstract に将来形 / 提案形** ("We will show") — 現在形 ("We show")
4. **fair comparison なし** — 比較対象の条件を揃える (表で明示)
5. **負の結果を隠す** — 正直に書く方が accept されやすい
6. **Related work が自己引用だらけ** — 競合他者を 5 件以上
7. **Figure caption だけに主張がある** — caption の数値・傾向・解釈は本文にも書く
8. **pragmatic な数値なし** — "significantly" だけ → "3.2x (Fig. 7)"
9. **cover letter で journal への fit が未言及** — "why this journal"
10. **response letter が箇条書きなし** — reviewer comment を line-by-line で
11. **supplementary material に本論の main result** — 本文に移す
12. **1-page digest に詳細を詰め込みすぎる** — abstract に個別誤差
    (0.04%, 0.001%, ...) を列挙しない。abstract は
    "sub-percent accuracy" のような包括表現にし、個別値は本文・図説明へ
    移す。Warburg 型要素 + Cauer Ladder Network (CLN) のように既知の
    組合せは「既知」と明示し、novelty は parameter-free Galerkin coupling,
    Schur complement, SIBC/HOIBC surface envelope などの差分に置く。
    HOIBC / Warburg などの専門語は
    本文初出で citation を置く。caption は plotted quantity の同定だけに
    近づけ、問題設定・解釈・比較は本文で説明する。caption に書いた数値・
    傾向・結論は本文にも必ず書き、caption が長い場合は本文を優先して
    caption を短くする。`rank-(1,1)` のような
    insider shorthand は "one bulk mode and one surface mode" のように
    物理的に言い換える。ただし "one bulk mode + one surface mode" は
    「円形導体全体が基底 2 個で表せる」と誤読されやすいので、有限次の
    bulk CLN basis と追加 surface envelope/block を明確に分けて書く。
    図・本文・caption のいずれかで、CLN 側の基底数 (`N_b=2`, `N=10`
    など) と表面インピーダンス側の次数 (`p_H=0`, leading/zeroth-order
    SIBC, `p_H>0` の HOIBC など) を必ず明示する。`2-rung CLN + SIBC0`
    が主張なら、`N_b` と dc 項の含有/非含有をこねずにそのまま書く。
    `uniform dc term is not included in this count` のような否定形の数え方は
    1-page digest では読者の負荷になる。`p_H=0` は
    high-order ではないので、現在のベンチマークを HOIBC 成果として
    売らず、「この例では leading SIBC term、高次項は曲面・3D などで
    必要になり得る拡張」と切り分ける。`N_b` は何を数えるかを必ず定義し、uniform dc term を
    含めるのか、zero-boundary bulk correction functions だけを数えるのかを
    書く。`N_b` と `N` を同じ図・節で使う場合は、`N_b` は混合モデルの
    体積補正基底数、`N` は体積のみ CLN 梯子の段数、というように役割を
    明示して表記ゆれに見せない。L-term / R-term は本文で
    L-terminated = inductive last rung, R-terminated = resistive last rung
    と定義する。CLN+表面インピーダンス結合の提案を示す digest で `N_b=1` だけを出すと、
    「DC + IBC だけではないか」と読まれやすいので、少なくとも bulk
    補正基底を 2 個使うか、`N_b=1` が十分な理由を本文で説明する。
    混合 Galerkin / Schur 補の中心式は `\[` ではなく番号付き `equation`
    環境に置き、本文で `\eqref{...}` する。`K_{bb}`, `K_{sb}` などの
    ブロックは `K_{bb}(s)` のように周波数依存を明示し、体積・表面・結合
    Galerkin ブロックのどれかを式の前後で定義する。Schur complement は
    algebraic elimination の結果であり、Dirichlet-to-Neumann / Steklov--Poincare
    map はその作用素としての解釈なので、同じ文でいきなり同一視しない。
    まず「bulk 変数を消去すると Schur 補が得られる」と書き、次に
    「それが surface Dirichlet data を Neumann flux に写す」と説明する。
    "high-frequency SIBC scaling" とだけ書くと曖昧なので、意図が
    admittance の `f^{-1/2}` / `s^{-1/2}` tail なら指数を明示する。
    CLN と Warburg/表面項の接続問題は「Galerkin 縮約の外」と狭く書かず、
    「有限段 CLN と接続する際に経験的な遷移周波数でモデルを閉じる」
    という汎用的な課題として書く。円形導体などのベンチマークは、半径だけ
    でなく導電率・透磁率・参照解まで本文に書き、再現可能にする。
    "wall band" / 「壁帯」は少なくとも和文では伝わりにくいので、
    "skin-effect transition region" / 「表皮効果の遷移領域」を使う。
    球・立方体・多面体などの secondary geometry は 1ページ digest では
    原則落とす。
13. **Human review の汎用知見を使い捨てる** — coauthor / human review で
    出た指摘のうち、特定原稿だけでなく他の paper / digest / slide にも
    再発しそうなものは、その場の修正だけで終わらせず、radia-mcp の
    paper-writing policy / skill / checker / tests に反映する。具体的には、
    (a) 文章ルールなら本ファイルと `_em_paper_style.py` に追加、(b) 機械検出
    できるなら `paper_writing_check_*` に warning として追加、(c) bad / clean
    test を最低 1 組追加する。原稿固有の数値・著者判断・未確定な約束は
    汎用ルール化しない。
    図については、本文の該当 section より前に main result figure を浮かせない、
    "Verification." を太字段落に埋めず必要なら `Numerical Example` などの
    番号付き章にする、図中の `$f_N$` などの基準線は caption または図中ラベルで
    何を示すか明示する、縦幅を潰して可読性を犠牲にしない、という指摘を
    digest 汎用ルールとして扱う。図は原則として `figure` 環境を使い、
    TeX ソース上も該当 section の直後に置く。`float` package の `[H]`、
    `\refstepcounter{figure}` + 手書き caption、`minipage` による非float図、
    `\clearpage` / `\newpage` による強制配置は最終手段であり、まず本文量・
    図幅・caption長・figureのソース位置を調整する。caption にしかない
    結果説明を残さず、本文の図参照段落に移してから caption を削る。
    full-paper 予告は、実際に検証計画がある場合だけ
    "will examine extension to three-dimensional conductors" 程度に留め、1-page
    digest の現在成果として 3D 結果を匂わせない。
    In research plans and proposal notes, do not title sections with the
    author's drafting activity ("plain explanation", "notes", "miscellaneous").
    Use section titles that state the theoretical or design role the reader
    gets from that section.  Do not present immature ideas as result-like
    integrated graphs; show them as an iteration flow or working hypothesis.
    Keep internal tool improvements, MCP implementation notes, and generation
    provenance out of public manuscripts; move them to development notes when
    they matter.  When claiming a coordinate transform or perturbation, state
    what remains invariant and what changes.  For differential-form arguments,
    say immediately after the equation whether the exterior derivative,
    conservation law, unknown, Hodge operator, metric, or coefficient changed;
    do not merely say "the nonlinearity is linearized."  For proposal-stage
    figures, prefer either (a) real calculation or measurement data that ground
    the proposal, or (b) an iteration flow showing how an unverified hypothesis
    will be tested.  Do not write defensively that a combined figure is "not
    ready"; state what evidence or procedure each shown figure contributes.
    Put unverified implementation ideas in a clearly marked inventory table or
    working-hypothesis list, not in a summary-result figure.
    座標形状を示す図は、視覚的な縦横比が物理形状を偽らないよう `axis equal` 相当を
    原則とし、最終 PDF に埋め込んだ状態で文字サイズ・縦横比・表示範囲を目視確認する。

---

## 📚 reference.bib lab style (2023 Compumag review 由来)

**ファイル名ポリシー (研究室標準, 2026-07-18)**: 書誌ファイルは
**`reference.bib`（単数）** に統一する。理由は「単数が正しい英語だから」ではなく、
paper-writing MCP の引用ツール (`cite_a_claim` の既定 `bib_path`、
`verify_citation`、`em_submission_gate` の `bib_policy`、`_em_paper_style`) が
この名前を **既定・強制** するため — ツールに合わせれば変更ゼロで全体が一貫する。
プロジェクト直下に **1つだけ** 置き、EN/JA 双方の `.tex` から同じ
`\bibliography{reference}` で参照する（単一ソース）。`references.bib`（複数形）や
`references-<suffix>.bib`（例: `references-es.bib`）は使わない。

**経緯**: 2023 年 Compumag → IEEE TMag 投稿 (`public-safe curated corpus`) で reference.bib の書き方に複数の review 指摘が入った。 以下を **default** style として全 paper / digest / poster で適用すること。

### Rule 0: `.bib` は毎回外部検索で裏を取る

`reference.bib` を確認・追加・修正するときは、ローカル
`.bib` や過去原稿を信じてそのまま流用しない。毎回、DOI / publisher page /
Crossref / 公式リポジトリ / 著者ページなど一次情報に近い source を検索し、
その文献が **現在の主張に対して適切な引用先か** と、metadata
(authors, title, journal/conference, year, volume, issue, pages/article number,
DOI) が正しいかを確認する。検索で裏を取れない repo-only citation は、
論文・preprint・公式 documentation が存在しないかを確認し、存在するなら
そちらを優先して引用する。

実装したアルゴリズムの出典には、そのアルゴリズムを実際に記述した査読付き論文を
引用する。隣接テーマのポスター、題名だけが近い会議要旨、著者名のない内部資料を、
実装根拠の代用にしない。内部資料は研究経緯の確認には使えても、公開論文の主引用には
しない。

外部確認後は、(1) bib entry の metadata、(2) 本文での引用位置、(3) reference
list の初出順、(4) abstract に引用が残っていないことを同時に確認する。

IEEE 系の近年論文では、通常のページ範囲ではなく article number / Early Access
metadata が使われることがある。`pages = {1--1}` や `PP(99):1-1` は通常の
ページ範囲として意味を持たない場合があるため、IEEE Xplore / Crossref で
volume, issue, article number, pages, online publication status を確認する。
article number が確定している場合は article number を優先し、まだ
`Volume PP, Issue 99, Page(s) 1-1` 型なら `note = {early access}` と DOI を
残し、`1--1` を実ページ範囲として見せない。

### Rule 1 (最優先): 著者の苗字は `{}` で囲んで case-mangling 保護

**Why**: IEEEtran.bst をはじめ多くの `.bst` ファイルは BibTeX 出力時に著者名の case 変換を行う (大文字保持しない、 全文小文字化、 surname を ALL CAPS 化など)。`{}` で囲んだ文字列は `.bst` が触らないので、 `McDonald` → `Mcdonald` や、 アクセント付き surname `Bíró` → `bíró` の歪みが起きない。

```bibtex
% ✅ GOOD: surname を {} で保護
author = {K. {Hollaus} and J. {Sch{\"o}berl}}
author = {Y. {Sato} and H. {Igarashi}}
author = {H. A. {van der Vorst}}
author = {Niels {K{\"o}ster} and Oszk{\'a}r {B{\'\i}r{\'o}}}

% ❌ BAD: 苗字 protect なし — .bst によって case が壊れる
author = {Yousef Saad}          % "saad" になる可能性
author = {A. McDonald}          % "A. Mcdonald" になる可能性
author = {Niels K{\"o}ster}     % accent が壊れる可能性
```

- **Lastname, Firstname 形式** (`Sato, Y.`) は comma の前が surname と BibTeX が認識するので、 比較的安全だが、 厳密な case 保証は `{}` 必須
- **複合姓** (`van der Vorst`, `de la Cruz`) は全体を `{}` で囲む
- **複合姓 + アクセント** (`Bíró`, `Köster`, `Sch{\"o}berl`) は `{B{\'\i}r{\'o}}` 等、 アクセント macros ごと brace で囲む

### Rule 2: 1 paper = 1 entry (重複 entry 禁止)

同じ論文に複数 cite key を作らない。 同じ DOI が複数の `@article{...}` に書かれていると bibtex は両方を reference list に出力し、 **同じ論文が連番 [3] と [14] に二重に出る**。

```bibtex
% ❌ BAD: Kameari 2018 が 2 entries
@ARTICLE{Kameari, author = {...}, doi = {10.1109/TMAG.2017.2743224}, ...}
@article{kameari2017cauer, author = {...}, doi = {10.1109/TMAG.2017.2743224}, ...}
% → references 欄に同じ論文が二重に
```

追加前に `grep "<doi>" reference.bib` で重複確認。

### Rule 2.5: Numbered references follow first citation order

For IEEE-style numbered references and manual `thebibliography` lists, order
the reference list by the first appearance of each `\cite{}` key in the text.
After moving related-work sections, tables, or citation paragraphs, re-check
the first-citation order against the bibliography.  Mixing author order, year
order, and drafting order breaks the meaning of in-text [N] labels and makes
the prior-work flow harder to follow.

### Rule 3: Cite key は `AuthorYearTopic` 形式に統一

```bibtex
% ✅ GOOD
@article{Kameari2018, ...}
@article{Sugahara2017team28, ...}
@article{Koester2021PGDCLN, ...}
@article{BiroKoester2022scalar, ...}

% ❌ BAD: zoo
@article{Karl, ...}              % 苗字のみ
@article{POD, ...}               % acronym
@article{Lanczos, ...}           % topic only
@article{CLN_open, ...}          % topic-only
@article{MOR_T.H&S.C, ...}       % `&` は BibTeX special char — parse error
@article{MOR_Y.STO&H.Igarashi, ...} % 同上
```

- 同一著者 + 同一年は **小文字 suffix** で区別: `Sugahara2017team28`, `Sugahara2017kelvin`
- 共著者の場合 first author を採用 (`Biro` first author なら `BiroKoester2022scalar`)
- `&`, `空白`, `日本語` は cite key に含めない

### Rule 4: Journal name を一貫した abbreviation に統一

```bibtex
% ❌ BAD: 同じ journal に 3 つの書き方
journal = {IEEE Transactions on Magn.}       % 中途半端
journal = {IEEE Trans. Magn.}                % full abbrev
journal = {IEEE Transactions on Magnetics}   % full

% ✅ GOOD: paper 単位で 1 形式に統一
journal = {IEEE Transactions on Magnetics}
```

- **IEEE 推奨は full name** (`IEEE Transactions on Magnetics`)。 IEEEtran.bst は full でも自動で短縮しない
- abbreviation 使うなら IEEE 公式 `IEEE Trans. Magn.` (period 含む)
- COMPEL → `COMPEL --- The international journal for computation and mathematics in electrical and electronic engineering` (full 標準)
- 旧 entries も新 entry も **paper 全体で 1 形式に揃える**

### Rule 5: 著者は `and` 区切り (comma 区切り禁止)

```bibtex
% ❌ BAD: comma 区切り — BibTeX は最初の "H. Karl" を "Karl, H." 解釈してしまう
author = {H. Karl, J. Fetzer, S. Kurz, G. Lehner, and W. M. Rucker}

% ✅ GOOD: and 区切り
author = {H. {Karl} and J. {Fetzer} and S. {Kurz} and G. {Lehner} and W. M. {Rucker}}
```

### Rule 6: 正しい entry type を選ぶ

```bibtex
% ❌ BAD: web リソースを @article で書いてある
@ARTICLE{TWP28_1,
  author = {...},
  title = {Description of TEAM workshop problem 28},
  howpublished = {\url{http://ics.ec.lyon.fr/team.html}}
}

% ✅ GOOD: @misc + howpublished + url
@misc{TWP28_1,
  author       = {H. {Karl} and J. {Fetzer} and S. {Kurz} and G. {Lehner} and W. M. {Rucker}},
  title        = {Description of {TEAM} Workshop Problem 28: An Electrodynamic Levitation Device},
  howpublished = {\url{http://ics.ec.lyon.fr/team.html}},
  year         = {2009},
  note         = {accessed YYYY-MM-DD}
}
```

| Entry type | 用途 |
|---|---|
| `@article` | journal / transaction 論文 |
| `@inproceedings` | 学会 proceedings 論文 |
| `@book` | 書籍 |
| `@incollection` | book chapter |
| `@techreport` | 技術報告書 |
| `@phdthesis` / `@mastersthesis` | 学位論文 |
| `@misc` + `howpublished` + `url` | web リソース、 standard 文書 |
| `@unpublished` | preprint / submitted / in preparation |

### Rule 7: Year は **published volume/issue の発行年**

```bibtex
% ❌ BAD: preprint や online-first の年を採用
@article{kameari2017cauer,
  year   = {2017},        % 実際は Vol 54 No 3 March 2018 publication
  volume = {54},
  number = {3}
}

% ✅ GOOD
@article{Kameari2018,
  year   = {2018},
  month  = mar,
  volume = {54},
  number = {3}
}
```

IEEE は早期公開で year ずれが頻発する。 **必ず published volume/issue の月の year** を採用。 `month` は BibTeX macro (`jan`, `feb`, ..., `dec`) を使う (`{Mar}`, `{March}` 等の文字列は inconsistent)。

### Rule 8: 必須フィールドの欠落禁止 + 空 brace 禁止

```bibtex
% ❌ BAD: 空 brace は欠落と同じ
@article{Foo,
  volume = {},     % 空 → 削除すべき
  number = {},
  ISSN   = {}
}

% ❌ BAD: 必須 field 欠落
@inproceedings{Bar,
  author = {...},
  title  = {...}
  % booktitle, year, pages 欠落 → references 欄に "in ?, ?? (n.d.)" と出る
}

% ✅ GOOD: 必須 field を完全に
@article{Foo2024bar,
  author  = {...},
  title   = {...},
  journal = {...},
  year    = {2024},
  volume  = {12},
  number  = {3},
  pages   = {1100304},
  doi     = {10.1109/...}
}
```

| Entry type | 必須 |
|---|---|
| `@article` | author, title, journal, year, volume, **pages or art. no.** |
| `@inproceedings` | author, title, booktitle, year, pages |
| `@book` | author, title, publisher, year, **address** (出版地) |
| `@misc` | author, title, year, howpublished (or url) |

### Rule 9: 文字エンコーディング — UTF-8、 mojibake 禁止

```bibtex
% ❌ BAD: cp932 から崩れた Japanese mojibake が残っている
@ARTICLE{Hiruma,
  journal = {�d�C�w�� �d���E��͂̍����x���Z�p�������ψ����}
}

% ✅ GOOD: 英訳 + (in Japanese) annotation
@article{Hiruma2017lanczos,
  author  = {S. {Hiruma} and H. {Igarashi}},
  title   = {On {L}anczos algorithm for non-self-adjoint matrices, Part 2},
  journal = {Technical Committee on Static Apparatus and Rotating Machinery, IEEJ},
  year    = {2017},
  month   = aug,
  note    = {in Japanese}
}
```

`.bib` は UTF-8 で保存。 cp932 から copy & paste した日本語は mojibake になりがち。 英訳 + `note = {in Japanese}` で対応。

### Rule 10: タイトルの大文字保護 (固有名詞 / acronym)

```bibtex
% ❌ BAD: BibTeX/IEEEtran.bst が小文字化する
title = {Cauer Ladder Network Representation of Eddy-Current Fields}
% → "Cauer ladder network representation of eddy-current fields" になる

% ✅ GOOD: 固有名詞 / acronym を {} 保護
title = {{Cauer} Ladder Network Representation of Eddy-Current Fields}
title = {Generating a {C}auer Ladder Network Representation of Eddy Current Fields}
title = {Application of {A}-input {C}auer Ladder Network Method to {MOR}}
title = {Cauer Ladder Network Representation of the {T}eam Problem 28}
```

`.bst` ファイルが title-case を強制する場合 (IEEEtran.bst デフォルト)、`{}` で囲まれた文字は **そのまま大文字保持**。 固有名詞 (Cauer, Padé, Lanczos, Kelvin, ...) と acronym (CLN, MOR, BEM, FEM, IEEE, TEAM, SIBC, ...) は全て `{}` で囲む。

### Rule 11: 投稿前 PDF の References 欄を音読確認

bibtex compile 後、 PDF の References ページを **連番で 1 つずつ音読** する。

- ✅ 著者名が正しい case で出ているか (`Sato` `Köster` など)
- ✅ Journal name が一貫しているか
- ✅ Page 番号が連続しているか (`1--4` / `1100304` 等)
- ✅ Year / volume / number が一致しているか
- ✅ Title の大文字保護が効いているか (`Cauer` が `cauer` になっていないか)
- ✅ 重複している論文がないか

**実例**: 2023 Compumag review (Reviewer 3) で `reference [14]` に typo 指摘が入った。 compile 後 PDF の音読を skip すると **bibtex で正常 parse できる typo** (例: `proceeings`, `Trnasctions`) は検出できない。

### 自動化 (radia-mcp 既存 tool)

- `bibliography_lint(bib_path)` — cite-key 規約 + 必須 field 欠落 + lab style 違反 を検出
- `bibliography_cite_validation(tex_path, bib_path)` — `\cite{key}` の key が bib に存在 / orphan entry 検出
- `bibliography_dedupe(bib_path)` — DOI / title fuzzy-match で重複 entry 検出
- `bibliography_normalize_journal_names(bib_path)` — journal name の inconsistency flag

**今後追加候補** (まだ未実装):
- `bibliography_check_surname_braces(bib_path)` — Rule 1 違反 (`{}` 無し surname) を検出
- `bibliography_check_title_acronyms(bib_path)` — Rule 10 違反 (acronym が `{}` 無し title) を検出
- `bibliography_check_authors_and_separator(bib_path)` — Rule 5 違反 (comma 区切り author) を検出

---

## 💡 実例集 (汎用 template)

### Case 1: 国内和文論文誌 (IEEJ / 電気学会 論文誌 クラス)
- 構造: 和文 8 ページ、図 6 枚、引用 25 件
- Hook: 1 段落で「既存手法の限界 → 提案」を定量的に示す
- Contribution list: 3 件を intro 末尾に箇条書き
- Reviewer 対応: 反論 figure を必ず 1 枚追加する余裕を持つ

### Case 2: 国際英文誌 (IEEE Transactions)
- 構造: 英文 6 ページ double-column、図 8 枚、引用 35 件
- Abstract: 200 words、problem→method→result→impact の 4 文構造
- Cover letter: "This work addresses a long-standing tradeoff in [domain]"

### Case 3: 加速器系 (JACoW / PRAB / PASJ)
- 具体的なvenueと投稿年度を先に確定し、公式テンプレートを使う
- 電磁石・電源・ビーム運転のどの問題を改善するかをIntroductionで明示する
- solver精度だけでなく、調整時間、再現性、運用手順への効果を示す

## ✅ 完成検証チェックリスト

### Abstract
- [ ] 4 要素 (問題/手法/結果/意義) が揃っている
- [ ] 単独で論文の主張が分かる
- [ ] 語数 (EN 200 / JP 400) 以内
- [ ] 弱気修飾語ゼロ
- [ ] **数式ゼロ** (display `\begin{equation}` / `\[...\]` / inline `$...$` すべて)
- [ ] **引用ゼロ** (`\cite{}` / `[N]` 全形式) — 例外なし、共著者・創設論文も本文に移す
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
- [ ] Caption に書いた数値・傾向・解釈が本文にも書かれている
- [ ] Methods の tense が過去形 (我々の行為)
- [ ] Discussion の hypothesis が現在形、自分の結果が過去形
- [ ] Bullet の文法が並列
- [ ] Phase 0 artifacts 5 点揃い済

### R. 使える診断ツール (paper_writing_*)

既存 + 新規で計画されている:
- `paper_writing_check_subject_verb_distance` — 主述の物理距離
- `paper_writing_check_paragraph_length` — 段落字数の範囲
- `paper_writing_bilingual_readability_check(text_or_tex_path)` — **和文と英文を別基準で査読者視点から評価**。和文は字数・読点・接続・述語・段落内の論証役割、英文は語数・従属節・名詞化・句読点連鎖を使う。数式・図表・keyword・TeX command は除外し、`japanese` / `english` に別の diagnostic score と source line を返す。両言語の score は平均せず、悪い側を `paper_writing_em_submission_gate` の判定に用いる。
- `paper_writing_check_abstract_background_ratio` — abstract 内 background 比率
- `paper_writing_check_abstract_no_math_no_citation` — **abstract に数式 (TeX math)、`\cite{}` / `[1]`、domain-specific acronym (FEM/BEM/MCP/LLM など) が混入していないかチェック** (IEEE / Elsevier / Springer / Nature / Science 共通の慣習: abstract は self-contained + 検索エンジン indexing 可能であるべき)。Display math (`\begin{equation}`/`\[...\]`) と citation は `status="fail"`、inline math (`$...$`) と acronym は `status="warning"`、両方なしなら `status="clean"`。
- `paper_writing_check_digest_human_review_triggers(tex_or_text)` — **1-page digest の human review trigger を検出**。abstract の個別 percent error 列挙、Warburg/CLN 既知事項を novelty として見せる書き方、経験的な接続周波数の問題を Galerkin 限定に見せる表現、HOIBC/Warburg 初出 citation 不足、説明過多 caption、`rank-(1,1)` 型の不透明 shorthand、「基底 2 個で円形導体を表せる」と誤読される minimal-basis 表現、重要な Schur/ブロック式の無番号・未ラベル・本文未引用、`K_{bb}`/`K_{sb}` などの `(s)` 依存や体積・表面・結合ブロック定義不足、`N_b` の定義不足、`N_b=1` が DC+IBC と誤読される構成、CLN 基底数 / HOIBC 次数の未記載、円形導体ベンチマークの導電率・参照解不足、"wall band" / 「壁帯」の不明瞭語、main result figure が Numerical Example / Verification より前に浮く配置、検証節が太字段落に埋まる構成、`$f_N$` など未説明の基準線、球・立方体・多面体まで盛る scope creep を `status="warning"` で返す。`paper_writing_em_submission_gate` からも自動実行。
- `paper_writing_check_undefined_acronyms(tex_path)` — **略語 (IH, MQS, FEM, BEM, ...) が初出時に full name と並記されているかチェック**。`Full Name (ACRONYM)` または `ACRONYM (Full Name)` パターン (初出の ±80 文字以内)、または Nomenclature / Acronyms / Abbreviations section に listed があれば OK。万人共通の略語 (PDF, USA, CPU, USB, ...) のみ whitelist、研究室 EM 専門用語 (FEM/BEM/MQS/IH) は意図的に whitelist 外 → 必ず spell out 必要。`extra_whitelist="ABC,XYZ"` で institutional 略語追加可。
- `paper_writing_check_citation_keys_exist(tex_path, bib_path)` — **`\cite{key}` の key が `.bib` の entry に存在するか静的チェック**。`status="fail"` = 引用キーが bib にない (compile 時 `[?]` で render される)、`status="warning"` = bib にあるが cite されていない entry あり (cleanup 推奨)、`status="clean"` = 1-to-1 一致。`\input{}` chain も自動 resolve (`auto_resolve_inputs=True` default)。bibtex compile 前の sanity check として使う。
- `paper_writing_check_ref_label_consistency(tex_path)` — **`\ref{}` / `\eqref{}` / `\autoref{}` / `\cref{}` / `\pageref{}` の key が `\label{}` に対応するか静的チェック**。`status="fail"` = dangling ref (PDF で `[??]` 表示)、`status="warning"` = orphan label (本文で言及していない figure/eq/table = digest で空間浪費)、`status="clean"` = 1-to-1 一致。digest review の典型指摘 "Fig. 3 is never referenced" を pre-compile で catch。
- `paper_writing_check_ieee_keywords(tex_path)` — **`\begin{IEEEkeywords}` (IEEEtran) の存在 + 個数 (3-7 推奨) + 各 keyword 長 (3-50 chars) チェック**。これはIEEEプロファイル専用である。IEEJの `jkeyword` / `ekeyword` は `paper_writing_em_submission_gate(target_venue="電気学会...")` が別に検査し、加速器系は選択した公式テンプレートの規則に従う。`status="missing"` = block 不在、`status="warning"` = 個数 or 長さ問題、`status="clean"` = OK。
- `paper_writing_check_pdf_unresolved_markers(pdf_path)` — **compile 後 PDF を pymupdf で text 抽出し、`[?]` / `[??]` rendered marker を検出**。上 2 tool (cite key / ref label) は pre-compile 静的 check、本 tool は post-compile の safety net (bibtex 再 run 忘れ等で漏れた未解決参照を catch)。各 marker の page 番号 + ±60 char context を返す。digest 提出直前の最終チェックに最適。
- `paper_writing_check_tense_consistency` — discussion の 3 部時制
- `paper_writing_check_figure_caption_showing` — caption が showing vs telling か、caption の主張を本文にも置く方針を返す
- `paper_writing_check_strong_adjective_budget` — 強調副詞の過剰
- **`paper_writing_check_word_repetition`** (v0.8.0) — 同一単語の近接障害 (中島・塚本)
- **`paper_writing_check_sentence_ending_variety`** (v0.8.0) — 文末表現の単調さ (中島・塚本)

### 和文・英文の readability は別物として扱う

同じ「長文」でも、和文では読点の連鎖、述語の遅延、定義・条件・結論を
一文へ入れる構造が負荷になる。英文では 30 words を超える文、従属節の
入れ子、nominalisation と technical noun の積み重ねが主な負荷になる。
このため、文字数と word 数を一つの平均へ混ぜない。

投稿前は `paper_writing_bilingual_readability_check(tex_path)` を実行し、
`HIGH` を先に修正する。修正順は、(1) 段落冒頭に reviewer takeaway、
(2) 理由、(3) 条件・定義、(4) 数値・記号である。レイアウト、参照、数式が
clean でも、この check が `fail` なら「技術的には正しいが読み手が論旨を
復元しなければならない」状態なので投稿を止める。

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

## 📚 追補: 作文技術 6 冊からの新規ルール (2026-06 re-OCR mined)

Wallwork / 佐藤『なぜあなたの研究は』/ 本多『日本語の作文技術』/ 中島・塚本『知的な科学・技術文章』/ 問題な日本語 / 作図力学 を再OCR全文から再学習。既存と重複しない 73 規則。
> ⚠️ **語長・同一語反復・態は EN/JA・section で文脈分岐** させており、既存の定量 lint 閾値（≤60字 / no-word-repetition / active-voice）を上書きしない。和文 60 字目安は既定、論理的に必要な長文は文節境界に述部+接続詞を入れる条件で可。英文はキーワード反復可・Method は受動可。

### 作文技術フォルダ 12 冊の取込境界 (2026-08 再監査)

`09_作文技術` 直下の OCR backup を除く 12 冊を再監査。

- **paper-writing に直接反映**: 『理科系の作文技術』、『まんがでわかる理科系の作文技術』、『知的な科学・技術文章の書き方』、『日本語の作文技術』、『中学生からの作文技術』、『問題な日本語』、Wallwork『論文の書き方・アクセプト術』、『なぜあなたは論文がかけないのか』。漫画版と入門版は別の規則数として水増しせず、原著の誤読防止と例の交差確認に使用。
- **figure/paper に反映**: 『作図力学』。図表選択、本文と図の重複排除、掲載寸法での可読性は `radia_mcp.figure` と共有。
- **presentation に反映**: 『研究発表のためのスライドデザイン』。論文本文へ「体言止め」を誤輸入せず、スライドの短文表現に限定。
- **grant-writing に反映**: 『申請書作成の理論』、『いちばんわかりやすい科研費申請書の教科書』。審査項目への応答など申請書固有の技法は paper に持ち込まず、読者の判断負荷を下げる原則のみ共有。

### 書き始める前に「目標規定文」を 1 文で置く

木下の目標規定文は、論文内にそのまま印刷する定型句ではなく、執筆中の取捨選択の判定基準。実験・解析の結論を確定した後、例えば「本稿は、[対象] に対し [手法] を用い、[定量的に何が分かったか] を主張するために書く」と 1 文に固定する。

- 目標規定文に寄与しない材料、図、枝葉の実装説明は削除または補遺へ。
- 執筆中に主張が変わったら、目標規定文だけを差し替えず、Title/Abstract/Conclusion と本文の材料配置を再監査。
- 診断は `paper_writing_title_abstract_conclusion_triangle` と、下記の `paper_writing_check_conclusion_first_use` を併用。

### 結論の「新しさ」は洞察であり、本文未出の証拠ではない

Wallwork 第 19 章は、結論を単なる要約とせず、結果の統合から得られる新しい視点・影響・将来への道筋を示すことを求める。ただし、その展望は論文中で示した結果から正しく導かれる必要がある。

- **可**: 本文で示した結果の重要度づけ、分野への含意、限界と対応する現実的な将来展望。
- **不可**: 結論で初めて現れるソルバ名、略語、変数、定量値、比較条件、引用。これらは Method/Result/Discussion のいずれかで導入・裏付け。
- `paper_writing_check_conclusion_first_use` で、結論で初出となる acronym-like な技術語、数式記号、数値、citation key を機械抽出。引用不要の一般略語と意図的な展望は whitelist/目視で除外。

## 📐 図表の使い分けと最小化 (作図の技術)

論文中の「図 vs 表 vs 文章」をどう選ぶか。`radia-figure` server は
作図そのものの品質 (font / no-in-figure-title / size) を見るが、
ここでは **そもそも表にすべきか・図にすべきか** という上流の判断を扱う。

### 結果は「まず図化できないか」を検討してから表にする

- **原則**: 解析結果は『文章では煩雑・図でも表現しにくい』ときに **限って**
  表にする。多くの結果は表より図の方が理解が容易。表は数字の箇条書きで
  定量把握はしやすいが一見では理解できず、読者に脳内再構築を強いる。
- **手順**: 結果表を作ったら必ず一度図化を検討する。図にできるなら図を
  優先し、細かい数値が要ると言われたら図のプロット点近傍に数値を添える。
  文章・表・図の 3 表現を比べて最も分かりやすいものを選ぶ。
- (作図力学 5.2)

### 表は最小限 (5〜6 ページで最大 1〜2 枚)

- 表の掲載頻度は図よりはるかに低いのが通例。表が所狭しと並ぶ論文は
  知的な書き方とはいえない。表は本当に表でしか提示できない情報に絞り、
  多数の表は図や文章に振り替えられないか見直す。
- (作図力学 5.序説)

### 線図の変化過程を文章で逐一なぞらない (図面-文章重複障害)

- **NG**: 一見して分かる線図の変化 (「A が増えると B が減り…」) を文章で
  重複させる → 「図を見れば分かる」と見識を疑われ、最も主張したい新知見が
  ぼける。
- **OK**: 読者が変化過程を了解している前提で、図で主張したい点 (例: 従来説を
  覆す箇所) だけを強調する。「図 3 より明白なように…である」と要点だけ書く。
- (作図力学 2.5.2 / 中島・塚本)

---

## 🖼️ 図の実装は radia_mcp.figure 正準で (2026-07-17)

論文の**全グラフ図**は `radia_mcp.figure` の `paper_figure()` +
`emit_paper_figure()` で作る (スタイル詳細は radia-figure server の
`paper_figure_recipe` / `paper_figure_quality_rules` が正典; ここには
**論文ワークフロー側の運用ルールと実地の落とし穴**だけを置く)。

### Rule 1: 1 枚指摘されたら全 figure を監査する

図の品質指摘 (単位の角括弧・凡例重なり・小さすぎ等) は**指摘された図だけの
問題であることはまずない** — 同じ生成コード/習慣で作った全図が同罪。指摘を
受けたら文書内の **figure を全数チェック**し、生成スクリプトごと正準化する
(2026-07-17 SA-26-078: 3 枚の指摘 → 実際は 6 枚全部が非正準で、生成スクリプト
の無い raster 図も 1 枚発掘された)。図ごとに**再生成スクリプトを必ず持たせる**
(生成コードの無い図は監査も再現もできない)。

### Rule 2: 作図幅 = 埋め込み幅 (クラス実測値でプロファイルを作る)

- 紙面 10 pt を守る前提は「作図した幅のまま 100% で埋め込む」こと。まず
  文書クラスの実寸を測る:
  `\typeout{COLW=\the\columnwidth TEXTW=\the\textwidth}` を挟んで 1 回
  platex → log から mm 換算 (IEEJ-tec.cls 実測: 82.17 / 174.35 mm。
  IEEJ Trans 用既成プロファイル 88/180 mm とは違う)。
- `emit_paper_figure()` は**プロファイル幅へ強制リサイズ**する
  (`paper_figure(rel_width=...)` は上書きされる)。実測幅の
  `dataclasses.replace(IEEJ_SINGLE_COLUMN, width_mm=82.17)` プロファイルを
  作り、`paper_figure()` と `emit_paper_figure()` の両方に**インスタンスを
  渡す** (`get_profile` は文字列以外にインスタンスも受ける)。
- tex 側は `width=\columnwidth` / `width=\textwidth` / 実寸 mm で 100% 埋め込み。

### Rule 3: 1 段組に多パネルを詰めない → figure* に格上げ

サブパネル幅が **~40 mm を切る配置は不可** (例: 84 mm 段に 3 ループ図 =
各 26 mm)。横並び多パネルは**両欄 figure\*** (実測 \textwidth) に格上げする。
その際、近接する単独図 (誤差 vs 振幅など) を **(d) パネルとして統合**すると
紙面消費が旧 2 図とほぼ相殺され、ページ限界内に収まる。

### Rule 4: auto_tighten の罠 — 全ラベル格子は warn + 明示マージン

`emit_paper_figure(on_fail='auto_tighten')` は効率目標 (0.72) に**届かない
場合ラベルが切れるまでマージンを削り続ける**。全パネルに軸ラベルの付く
2x2 / 1x4 格子は構造的に効率 0.5-0.65 で頭打ちなので、`on_fail='warn'` +
明示 `fig.subplots_adjust(...)` に切り替え、効率警告は受容する。
`aspect` は**行あたり** (総高さ = 幅 x aspect x nrows) な点にも注意。

### Rule 5: ゲートに素直に従う (実際に効いた例)

- 凡例重なり検出 → 凡例をやめ**曲線終端の直接ラベル** (`ax.text` /
  `label_curve_endpoints`) へ。
- 色ゲート → 独自色 (#E8000B / tomato / steelblue / viridis 曲線族) は全滅。
  モデル=vermillion #D55E00・比較=blue #0072B2・基底族=**グレースケール**・
  参照=黒、で大抵足りる。キャプションの色名も追随 (赤→橙)。
- 単位は丸括弧 (A/m)。`[A/m]` は IEEJ/IEEE 規約違反。
- **重ね描き比較（モデル vs 参照）の線種**: 点線や疎な破線は「抜け」が
  大きく曲線として読めない。**dash ≥ 2.5×gap の密な破線**
  (例 `ls=(0, (5, 2))`) で上に重ね、下の参照線（黒/灰の実線）が gap から
  見える密度にする (Sugahara 2026-07-17:「点線より破線のほうが抜けが
  少なくてよい」)。

### Rule 6: 検証は「レンダして見る」まで

図差し替え後: `extractbb` (stale .xbb は寸法崩れ) → 再コンパイル →
**pymupdf で該当ページを画像化して目視** (ラベル欠け・パネル間衝突・
キャプション不一致はこれでしか見つからない) → raster 画像 0 / Type3
フォント 0 を確認 (`pdf_font_embed_check` / get_fonts)。「図ファイルが
正しい」ことと「紙面で正しく出る」ことは別 (PDF 一択の節も参照)。

再描画の重い図 (能動学習ループ等) には **`--plot-only`** パス (保存済み
JSON/選択点から決定論的に再構成) を用意しておくと、図スタイル反復が
分単位で回る。

---

## 🎤 発表スライド (presentation — 2026-07-17 統合)

**presentation は独立サーバをやめ、この paper-writing スキル群に統合された**
(菅原判断: スライド本体は現状 AI が end-to-end で作れない。単独サーバに値する
のは「AI が作れる」ものであって、スライドはまだ違う)。

- **運用の分担**: スライドは**人が作る**。AI の役割は lint / 抽出 / 予算配分 /
  台本照合 — `presentation_*` ツール群 (密度・箇条書き・発表時間 1/3・1/4
  ルール・タイトル動詞・数式/図スライド適合・台本カバレッジ・PPTX テキスト
  抽出・TTS 埋め込み等) は `mcp-server-paper-writing` が配信する。
- **論文とスライドは同じ主張の別レンダリング**: 数値・主張・図は論文側 (この
  スキルの誠実な定量報告・figure 正準) と同期して監査する。論文で打ち切り値を
  到達値と書かないのと同じ規律をスライドの箇条書きにも適用する。
- スライド作文の詳細原則 (speakability・単一メッセージ・理系ミニマリズム・
  発表時間規則ほか) は `radia_mcp/presentation/skill.md` を参照 (統合後も
  実装ホームは `radia_mcp.presentation` モジュール)。

---

## 🔁 フレーミング移行時の figure/数値監査ポリシー (2026-06-14)

論文の理論フレーミングを更新したとき (例: Warburg-Schur → Mixed Galerkin)、
**本文 (タイトル・abstract・章題・キーワード) と figure/* と figure 内の数値
は不可分の三点セット**として監査する。これを怠ると、概念は新フレームだが
**図と数値は旧フレームのまま**という最も悪いタイプのドリフトが発生する。

### 実例 (2026-06-14 lab incident)

IGTE 2026 ダイジェスト:
- ✅ Title: "Mixed Galerkin Reduction" に更新済
- ✅ Abstract / §3 本文: Mixed Galerkin で記述
- ❌ Fig. 1: 旧 `circle_warburg.pdf` のまま (Warburg-with-$d$ で 17% wall band)
- ❌ Fig. 1 caption: "rank-(1,1) Mixed Galerkin **specialization** — equivalent to the historical Warburg-Randles cell — within **17%**"

実際の no-$d$ Mixed Galerkin は **0.064%** (270× 改善)。本文には正しく
0.04% と書いてあるが、図と図キャプションだけが旧 Warburg 時代の 17% を
報告し、自己矛盾していた。reviewer がこれを見つけたら「結果が再現してない」
と即 reject。

### 監査ルール

**Rule 1: 移行は三点セット**

| 三点セット | 監査項目 |
|---|---|
| 本文 (text) | フレーム名・abstract 数値・章題の用語 |
| 図 (figures/*) | ファイルが新フレーム由来か / 内部の数値が新フレームの実測か |
| 数値 (caption + body) | 図キャプションと本文の数値が一致しているか |

何か一つだけ更新するのではなく、**3 つを同じコミット**で更新する。
本文だけ更新して図を後回しにすると、ほぼ確実に忘れる (人間も AI も)。

**Rule 2: 同じファイル名で意味が変わる場合は rename**

旧フレームの `figures/circle_warburg.pdf` を新フレームの内容で上書きしないこと。
Warburg-with-$d$ の図と Mixed Galerkin no-$d$ の図は**別ファイルとして共存**
させ、本文の `\includegraphics{...}` を差し替える:

- ✅ 旧: `figures/circle_warburg.pdf` → 残す (歴史的参照、git blame で追跡可)
- ✅ 新: `figures/cylinder_mixed_galerkin.pdf` → 新規追加
- `.tex` の `\includegraphics{figures/cylinder_mixed_galerkin.pdf}` で差し替え

**なぜファイル名を変えるか**: 同じファイル名が「Warburg ブランチでは Warburg
結果、Mixed Galerkin ブランチでは Mixed Galerkin 結果」と branch 依存で意味が
変わると、ブランチ切り替えや古いキャッシュ済 PDF を開いたときの取り違えが
起きる。**ファイル名に framework を埋め込む**ことで物理的に区別する
(`circle_warburg.pdf` vs `cylinder_mixed_galerkin.pdf`)。

**Rule 3: 移行コミットの差分 self-review**

フレーム移行コミットを作るときは `git diff` を以下の順に必ず確認:

1. `*.tex` の本文 (フレーム名・章題・abstract): 新フレームになっているか
2. `*.tex` の `\includegraphics{...}` パス: 旧ファイル名を指していないか
3. `figures/` の追加/変更: 新ファイルが追加されたか / 旧ファイルが意図的に
   残っているか
4. `*.tex` の数値リテラル (`17\%`, `0.064\%` 等): 図と本文で一致しているか
5. キャプション内の用語: 「旧フレーム名による特殊化」の記述が必要なら
   歴史的引用として残し、それ以外は新フレーム用語に置換

5 点全部が揃っていなければ、コミットを分割するか保留する。

### Cross-link

- 関連: `figure_style_guide` (figure 品質の純粋な技術ルール)
- 関連: `cross_lint` (本文と citation の整合性チェック)
- (本ポリシー自体には自動 lint なし: 上記 Rule 1-3 は人間 (or AI) が
  目視で行う運用ルール。将来 `paper_writing_check_framework_migration`
  tool 化を検討)

---

## 🔢 数値・単位・記号の組版規則 (英文 paper)

reviewer がスペル・組版の不統一を見ると「提示が素人」と判断し、
データそのものへの不信につながる。投稿前に機械的に揃える。

### 数値と単位の間隔: 略語単位はスペース、記号は詰める

- **単位の略語 (語の短縮)** は数値との間にスペース 1 つ: `34 kg`, `10 mm`,
  `2,400 mL`。これは元が単語 (`2,400 milliliter`) なので詰めると奇異。
- **単位の記号・%** は詰める: `20%`, `4°`。等号・不等号・± と数字も詰める
  慣例: `P<0.001`, `P=0.02`。
- **ただし数学・物理系で数式として書く場合は等号・不等号の前後にスペース 1 つ**:
  `x = y + 1`, `P < 0.001`。
- 名詞修飾の位置ではハイフンで連結: `a 5.2-kg infant`, `a 40-mm-diameter tumor`。
- (佐藤 Q27 / 木下 — SI 表記)

### 市販の機器・薬品は発売元を明記する

- 再現性のため購入元を特定できる必要がある (慣例)。機器/試薬名の直後に
  括弧で: 米国製は (メーカー名, 市名, 州名)、米国外は (メーカー名, 市名, 国名)。
- (佐藤 Q27)

### 量記号はイタリック・添字は下付き

- 量記号を立体にすると他記号と区別しにくい。変化量はイタリック + 下付き添字。
- 分数単位の分母は () で括る: 比熱は `[J/(kg·K)]`。`[J/kg·K]` は
  「kg で割ったか K で割ったか」が二義的。
- 英文ハイフネーションは辞書の音節境界で改行 (`con-clu-sion`)。
- (中島・塚本 基本ルール20 / 4.3.2)

---

## ✍️ 英文センテンス構造 (Wallwork 補完)

skill.md の Given-New (§D)・showing-not-telling (§H)・hedging を超えた、
語順・句読点レベルの実戦規則。日本語話者が直訳で崩しやすい箇所。

### 基本語順 主語+動詞+目的語 を 4 要素離さない

- 英語ネイティブはこの語順を厳密に期待し、崩れると思考の流れが分断される。
  `last week` / `for the second time` などの副詞句を主語・動詞・目的語の間に
  割り込ませない。
- (Wallwork §2.2, §2.8-2.9)

### 主語と動詞の間に 8〜10 語以上を挿入しない

- 主語と動詞の間に挿入された情報は読者に軽視され、動詞到達時には主語を
  忘れている。`The result, after the calculation has been made, can be used` →
  `After the calculation has been made, the result can be used`。
- 主語の導入を `It is...` で遅らせない: `It is probable that this is` →
  `This may be`。助動詞 (may/might/could/should) や副詞 (surprisingly) に
  書き換える。
- (Wallwork §2.4-2.7, §5.12)

### 否定語・目的を文頭近くに置く

- **否定語 (no/not/none)** はできるだけ文頭近くに: `...were not available` →
  `No data were available`。肯定的に始まって途中で否定に転じる文は読者を裏切る。
- **理由 (根拠) より目的を先に** 述べる: 読者はまず目的を知り次に手段を知りたい。
- (Wallwork §2.11-2.13)

### 形容詞・名詞の語順

- 形容詞は名詞の前に置く (または関係代名詞節で後置)。`a paper particularly
  interesting` → `This paper is particularly interesting` / `a paper that is
  particularly...`。`The main document contribution` → `The main contribution
  of the document`。
- **名詞を数珠つなぎにして形容詞を作らない**: `art state technology` /
  `mass destruction weapons` は誤り。3 語以上の名詞連結を避け前置詞・動詞で
  関連を示す (`weapons of mass destruction`)。考案した連結は Google Scholar で
  実在を確認 (ヒット 10 万件未満なら使わず前置詞・動詞へ)。装置・手順名
  (`an Oxford SATW EDX detector`) は例外的に連結可。
- (Wallwork §2.15-2.17, §12.8)

### キーワードは反復する — 類義語で置換しない

- 長文を分割するときキーワード反復は必然となり、それが可読性を高める正しい
  テクニック。`survey` や `English` などのキーワードを間を空けず繰り返す。
- **類義語に置換すると微妙な意味差を読者に憶測させ混乱を招く**。
  it/that/this/former/latter/which 等の代名詞より、指示対象のキーワードそのものを
  反復する (中島・塚本「同一単語の近接障害」とは逆方向の指針 — 英文では明確さ
  優先でキーワードを反復、和文では幼稚さ回避で類語化。文脈で使い分ける)。
- (Wallwork §4.4, §6.4, §20.8)

### あいまいさの最大原因は代名詞: it/they を先行詞より先に出さない

- 後方指示 (前文の名詞を指す) は文頭の代名詞 OK だが、前方指示 (後続の名詞を
  指す) は読者を待たせる。`Although it is a stable material, the composition of
  beeswax...` → `Although beeswax is a stable material, its composition...`。
- 共著者がリストの順序を変えると former/latter が別物を指す事故が起きる →
  キーワードを反復する。
- (Wallwork §2.10, §6.3, §20.8)

### ピリオドを最大限活用し、コンマ多用・セミコロンを避ける

- コンマの多用は無計画な手抜きライティングの兆候で、読者は思考の流れを何度も
  修正させられる。列挙以外のコンマ多用文は短文に分割する。
- 現代英語で情報追加のセミコロンはほぼ不要 (ピリオドで読者は無意識に小休止)。
  セミコロンは項目をグループ化して関連を示すとき (`Spanish, Italian and
  Romanian; German and Dutch; ...`) のみ。
- 説明的な括弧は流れを阻害し文を長くする。括弧は短い具体例リスト
  (`e.g. Spain, France, Germany`) に限定し、説明文や非リスト実例は入れない。
- (Wallwork §4.12-4.15, §6.9)

### and/which/接続詞/動詞-ing/in order to で文が伸びたら分割する

- これらの構造を多用すると文が際限なく伸びあいまいさが増す。`and` の前で
  ピリオドを打ち、必要なら次文を `Also/In addition` で始める。長い `which` は
  `this` (this fact/this method) に置換。
- (Wallwork §4.8-4.11)

### 制限用法と非制限用法で which/that を使い分ける

- コンマで挟む `who/which` 節は付加情報 (非制限)、`that` 節は限定情報
  (制限、コンマなし) で意味が変わる。動詞-ing形や which/that の先行詞が直前の
  正しい名詞になっているか確認する。
- (Wallwork §6.10-6.11)

---

## ✂️ 英文の冗長削減 (Wallwork §5)

skill.md の `check_prose_density` / `suggest_concept_drops` は **悪い圧縮** を
検出するが、ここでは **そもそも語数を減らす正攻法** を扱う。語数を減らすほど
ケアレスミス (前置詞・綴り) が減り要点が明確になる。

### 名詞より動詞を使う (補助動詞+名詞 → 1 語の動詞)

- `make a comparison` より `compare`。補助動詞の選択ミスもなく流れがスムーズ。
  `X was used in the calculation of Y` → `X was used to calculate Y`。
  `make an analysis`→analyze, `carry out a test`→test, `reach a conclusion`→conclude,
  `show an improvement`→improve, `undergoes a rapid rise`→rises rapidly。
- (Wallwork §5.13-5.14)

### 冗長な「一般語+具体語」と無意味な抽象語を削る

- `small in size`→small, `The process of registration`→registration,
  `twice a year in June and December`→`in June and December`,
  `two countries (Italy and France)`→`Italy and France`。
- `activity, case, character, condition, factor, instance, operation, phase,
  phenomenon, problem, procedure, process, situation, step, task` は伝達に
  貢献せず記憶に残らない → まず削除を試み、無理なら具体語に置換。ただし
  `freedom/love/fear` 等の明確な抽象語は残す。
- (Wallwork §5.3-5.5)

### 著者目線のメタコメントをしない

- 読者は `we observe that` / `we find it interesting to note` / `As we can see` /
  `in the rest of the paper we...` の連続を好まない。
  `As in the previous case we observe that there are three...` → `There are three...`。
  `It is now time to turn our attention...` → `The rest of the paper focuses on...`。
- 見出し直後・Results/Conclusion 冒頭の不要な導入語句 (`The salient results are
  summarized in the following.` / `In conclusion, we can say that...`) を削除し
  直接本論に入る。
- 評価形容詞 (`novel/interesting/appropriate`) は理由を直後に具体的に説明できる
  ときだけ使う。同義形容詞の連続 (`absolutely necessary`, `completely
  different`) は冗長。
- 「新規性を置く」「新規性を〜に置く」のように執筆者の配置作業を述べない。
  先行研究が含まない物理項、アルゴリズム、検証範囲を具体的に述べ、その差がもたらす
  効果を直接書く。
- (Wallwork §5.10-5.11, §5.15, §8.13)

### 原稿は基本的に短く (同内容なら 20 ページより 15 ページ)

- 原稿が長いほど主張があいまいになり主旨理解が困難になる。査読者は厚い原稿に
  否定的反応を示す (短い論文が引用されにくいという証拠もない)。40 ページが
  本当に必要な情報量か、推敲を避けた結果でないか自問する。25% 短縮を求められたら
  重要内容を保ちつつ削る (質はほぼ下がらず、むしろ改善する)。
- (Wallwork §5.20, §20.4)

---

## 🎬 トーン調整とヘッジングの精密制御 (Wallwork §10)

skill.md の hedging 禁止・強調副詞 budget を **段階調整** のレベルまで掘り下げる。
英米トップ誌では 100% の確信や尊大な態度は嫌われる一方、ヘッジしすぎると
日本人読者には謙虚でも欧米では「自信の無さ」と取られる。

### 主観的・未確認の仮説にはヘッジング、客観的事実には不要

- `prove/demonstrate/is` の代わりに `would seem to indicate / may be / could /
  appears to`。`we believe / as far as we know / to the best of our knowledge` を
  主張の前に置く。賛否が分かれそうな成果に **限り** 適用。
- `Table 2 shows that X had higher values` のような客観的事実にヘッジは不要。
- (Wallwork §10.2, §10.4, §18.9, §19.7)

### 確からしさの語を使い分ける (確率の目安)

- `must/cannot` = 100%、`may/might/could` = 中程度、`possibly/conceivably` =
  50-70%。`seem/appear/tend` を強い動詞の前に重ねる (`seems to show`)。
  `somewhat/relatively/quite/to a certain extent` で和らげる。`significantly` は
  「統計的に偶然でない」の意味でのみ使う。
- (Wallwork §10.4-10.7)

### ヘッジの強度を 1 文の中で混在させない

- `It is clear that yellow may be...` は `may` が `clear` を打ち消し信頼を欠く。
  4 つの婉曲表現を重ねると自信皆無の印象。断定すべきは断定し
  (`It is clear that yellow is preferable`)、和らげるべき箇所だけ `may` を使う。
- (Wallwork §10.11)

### 他研究を批判するときは建設的・肯定的に扱う

- 批判対象の提案者が査読者かもしれず、あからさまな批判は反感を招く。H1 を
  否定する前に H1 の良さを代弁する (「H1 提案以降に新データが出た」「H1 は
  小サンプルだった」「H1 著者自身が限界を認めていた」)。`although/however/
  moreover` の使い過ぎに注意 (否定的トーンが強まる)。
- (Wallwork §10.9-10.10)

### 自分の研究を `our` で明示する (Who+Did+What)

- 受動態は動作主を隠すため、著者自身の研究か他者の研究か読者が判別できなくなる
  (考察で頻発)。非人称スタイルでも `our` という語は自分の研究を示すために
  積極的に使う (`in our survey`)。we 禁止の誌では時制で区別 (自分=過去形、
  定説=現在形)。
- (Wallwork §7.2, §7.4, §7.6, §18.6)

---

## 🏷️ タイトル・要旨・各セクション (Wallwork 補完)

skill.md にタイトル・abstract の一般則はあるが、Wallwork の具体的な
書き換え規則・時制ルールは未収録。

### タイトルの 6 条件と前置詞・動詞-ing

- 効果的なタイトルは執筆スキルの半分。研究のユニークさを 3〜5 語で表現し、
  名詞の羅列を避け、動詞を最低 1 つ・前置詞を 2〜3 使う。
- 5 語以上なら前置詞 (by/for/from/in/of) で語句間関連を明確化:
  `Depression measuring inventory` → `Inventory for measuring depression`。
  `<名詞①+of+名詞②>` では名詞①に the。冠詞判断は Google Scholar で類似
  タイトル検索。
- 抽象名詞より動詞-ing形: `The Specification and Evaluation of...` →
  `Specifying and Evaluating...`。`novel/innovative` は何が新規か示せず誰も
  検索しないので避け、`low-cost/high-performance/pain-free` 等の具体形容詞を使う。
- 無駄な前置き (`A study of` / `An investigation into`) を削る。疑問形
  (`Will women always live longer than men?`) や 2 分割型 (コロンで前半に
  問いかけ・後半に専門説明) も候補 (特に学会抄録で効果的)。
- (Wallwork §12.1-12.7, §12.12-12.15)

### スペルチェッカーを過信しない (タイトル・キーワード)

- `hearth/form/recorder/through` 等の実在語の誤用はスペルチェッカーで
  検出されない。`company was funded`(正: founded)。著者は内容を熟知するため
  誤記を見落とす。スペルミス 2 つで査読者は完璧な校正まで受理を見送りうる。
- キーワードの綴り/句読点ミスは検索エンジンに検出されない → 第三者に確認させる。
- (Wallwork §12.18, §20.11)

### 要旨で「述べてはならない」もの

- 既知の常識的背景・無根拠の主張・キーワード定義・**数式**・**他論文への引用**・
  過度な数詞 (`many/several/few` は具体数値に) / 主観形容詞
  (`innovative/interesting/fundamental` は理由提示か削除)。
- 第一センテンスに陳腐な定型句 (`This paper deals with` / `The aim of this paper
  is to` / `This article explores`) を使わない (要旨は論文の広告)。
- 背景情報は全体の 1/4 以内に抑え、`In the last few years` / `Recently there has
  been` で始めない。
- (注: skill.md `check_abstract_no_math_no_citation` / Wallwork §13.15-13.22)

### 要旨・本文の時制を使い分ける

- 確立した現在の知見 = 現在形、過去から現在に至る背景 = 現在完了形
  (`In the last few years there has been` / `To date there has not been`)、
  研究成果 = 過去形 (`we found/we showed`)。要旨を躍動感ある現在形で統一する
  手法もある。
- (Wallwork §13.9, §13.30, §17.5)

### 結論の第一センテンスから定型前置きを削る

- `This paper describes` で始め `This paper has described` で終わるのは語数の
  無駄でメイントピックの開始が遅れる。最初の 5〜8 語を削除し、
  `In this study it is concluded that compression plays...` → `Compression
  plays...`。直接的すぎるならヘッジ (`reveals`→`seems to reveal`) で和らげる。
- 明確な結論を導けない場合も、発見できなかったこと自体が重要なので簡潔に
  述べ、`Although it is too early... two patterns seem to be emerging` のように
  将来展望を加えてネガティブで終わらせない。
- (Wallwork §19.7-19.8)

### 方法の手順は 1 文 1 手順にしない・箇条書きは自然な文に

- 単一手順=単一文を繰り返すと全文が名詞始まりで単調。連続する 2 手順を 1 文に
  統合する (ただし主語と動詞を分断する挿入は避け、3 つ以上の独立情報は別文へ)。
- 箇条書きは可読性を上げるが、声に出して不自然な羅列 (`Processes include:
  oxidation, hydration, ...`) は避け `Several processes occur, including A, B, C` の
  ように自然な文にする。時系列手順は番号、それ以外はブレット。
- (Wallwork §16.8-16.13)

### インフォーマルな英語を投稿前に排除する

- `doesn't`→does not, `kids`→children, `a lot/big/tiny/nice` を避ける,
  `so/till/like`→thus/until/such as, `check out/work out`→examine,
  `going to`→will/現在形, `you` の使用を避ける。
- (Wallwork §20.10)

### 投稿前の一貫性チェック

- 図表番号・キャプションと本文解説の一致、用語の統一、数式後のコンマ有無、
  `i.e./e.g.` 後のコンマ、米英綴りの統一をダブルチェック。一貫性の欠如
  (`Figure1` と `fig5a` の混在等) があると査読者は「提示が素人」と判断しうる。
- (Wallwork §20.3, §20.9)

---

## 📊 佐藤『なぜ』補完 — Result/Discussion の分量と棲み分け

skill.md は佐藤の北極星・Phase 0・Discussion A+B→C・Reviewer A-D を収録済。
未収録の **分量規律と細部** を追加。

### Result の各 Figure は 3 段で書く

- 各 Figure に対応して **①何を (なぜ) 調べたか → ②Figure の解説 → ③結果の
  まとめ (任意)**。①は特に簡潔に (1 文 25〜30 語まで、条件の細部は後を読めば
  分かるので冒頭で詳述しない)。③は結果が単純なら省略可。
- (佐藤 Q28)

### Result 本文と Figure legend で同じ文を重複させない

- 本文・legend がそれぞれ単独で読めて、両方読んでもくどくないようにする
  (3 条件で点検)。legend には矢印・スケール・略語の説明を入れ、本文の
  「上のパネルは〜下のパネルは〜」は legend にあるなら削る。本文語数が厳しければ
  Result を薄くし legend で詳述する裏ワザもある (legend は通常語数制限外)。
- (佐藤 Q28)

### Result と Discussion を文体で棲み分ける

- 断言できること (`A was B`) は Result、`〜と考える/推測する/示唆する` が必要な
  ものは Discussion。ただし確からしさが高く次の実験へ確実につながる強い結果なら
  Result 内でも `This result suggested...` でまとめてよい (固定ルールではなく
  確からしさの程度問題)。結論を否定するデータは絶対に隠さない。
- (佐藤 Q29-Q30)

### Discussion の分量は規定総語数の 1/4 (多くて 1/3)

- 総語 3,000〜4,000 なら約 1,000 語。1 段落 100〜150 語、第 1 段落と limitation と
  締めを除くと中心は 600〜700 語 = 3〜6 ポイント。これから段落数・文数が自動的に
  決まる (好きなだけ書けない)。データの出た順でなく結論への橋渡しに最適な順で
  組み替える。
- (佐藤 Q33)

### limitation は列挙でなく建設的コメントを付ける

- ただ挙げるだけだと取ってつけた印象。各 limitation に **①妥当性の主張** (限界は
  あるが○○の理由で妥当) か **②将来展望** (克服のため今後△△が必要) を付け、
  `Further investigation is necessary.` で終わらせず具体的に何が必要かを書く。
  研究目的外の事柄 (例: 副作用) を無理に limitation にしない。
- (佐藤 Q34)

### 知的作業 (Abstract/タイトル/Reviewer 返答) は朝のまとまった時間に

- Abstract・タイトル・Introduction 構成・Discussion 構成・Reviewer 返答は思考力を
  要する知的作業で、Methods/Result の単純作業とは別物。疲れた夜より朝に書くと
  明快になる。
- (佐藤 Q11, Q36)

---

## 🇯🇵 修飾の順序と句読点 — 本多『日本語の作文技術』

和文 paper の **わかりにくさの最大原因は修飾語と被修飾語の距離** および
**テン (読点) の打ち方の無原則さ**。中島・塚本 (近接障害・文末変化) を補完する
本多勝一の体系。まず自由に書き、読み返して「おかしい」と感じた箇所にだけ
これらの原則を当てはめる (暗記してから書くのではなく自己点検の道具)。

### かかる言葉と受ける言葉はできるだけ直結する (入れ子をはずす)

- 修飾語と被修飾語が離れると、文法的に正しくても一読では係り受けが取れず
  読み返しが必要になり、最悪は別の語にかかって正反対の意味になる。本書が
  最重要とする二大技術の一つ。
- 一文を書いたら各修飾語が「どの語にかかるか」を確認し、被修飾語の直前に
  置けるよう **語の位置だけ** を動かす。多重入れ子 (「私は小林が中村が鈴木が
  死んだ現場にいたと証言したのかと思った」) を作らない。係る側が出たら受ける側を
  必ず近くに置く (底抜け文章 = 受ける語の欠落を防ぐ)。
- (本多『日本語の作文技術』第二章)

### 修飾の語順 四原則 (優先順位 1>2>3>4)

1. **節 (述語を含む句) を先に、句 (述語を含まない句) を後に**: 句を先に置くと
   修飾が節の中の先頭名詞だけにかかる。「白い横線の引かれた紙」だと横線が
   白い意味になる → 「横線の引かれた厚手の白い紙」(節→句)。
2. **長い修飾語ほど前に、短い修飾語ほど後に**: 文字数の降順。短い「私は」
   「太郎さんが」は後ろに回す。原則 1 とほぼ同等の比重。(「主語・述語を近くせよ」
   式の俗説とは別物 — 主述は棚上げしてよい)
3. **大状況・重要内容ほど前に**: 長短が同程度なら最も大きな状況・重要な要素を
   先頭に (「初夏の雨がもえる若葉に豊かな潤いを与えた」)。
4. **親和度 (なじみ) の強い語は、直列なら近づけ、並列なら引き離す**: 隣接語と
   親和度が強すぎると誤読が起きる (「みどりがもえる夕日」で「みどりがもえる」と
   読む)。並列的にかかる親和ペアは離す。述語も相手を選ぶ (「雨が潤いを与えた」は
   翻訳調で不自然、「照り映えた」)。
- (本多 第三章)

### 形容詞は最後の名詞だけにかかる (語句全体にはかからない)

- 「美しい水車小屋の娘」より「水車小屋の美しい娘」。形容語を被修飾名詞の直前に
  置く。「積極的任務の遂行」より「任務の積極的遂行」。「西独製品が持つ強い
  価格面以外の競争力」は「強い」を「競争力」に直結。
- 否定述語にかかる修飾語は特に注意 (係り先取り違えが致命的): 「東京湾は瀬戸内の
  ようには…分離されていない」と限定の「ハ」を補う。
- (本多 第二章)

### テンの第一原則: 長い修飾語が二つ以上あるとき、その境界に打つ

- 句読点は字と同等以上に重要。同じ述語にかかる長い修飾語が複数あるとき、その
  境目に **だけ** テンを打つ。短い修飾語の境界には打たない (打つと重要なテンと
  区別がつかなくなる)。重文の境界も同じ原則に吸収。
- (本多 第四章)

### テンの第二原則: 原則的語順 (長い順) が逆順になっている場合に打つ

- 短い題目語「○○ハ」を冒頭に置くなど逆順にした文は、テンがないと係りが乱れる。
  「渡辺刑事は、血まみれになって逃げ出した賊を追いかけた」。逆順は筆者がその語を
  強調したい主観の表れ。倒置文もこの原則に吸収。
- (本多 第四章)

### 構文上必要なテンは上記二大原則だけ — それ以外は可能な限り打たない

- むやみにテンを打つと一つ一つの意味が弱くなり、本来構文上重要なテンの役割を
  侵害する (ゴチック多用が強調を無効化するのと同じ)。「自由なテン」(思想の最小
  単位を区切る・接続詞後の強調 `しかし、`) は強調したい所にだけ意識的に打つ。
- 文が終わったら必ずマル (句点) を打つ (字を抜く以上に重大)。連体形が終止形と
  同形のとき (`…ではない母は`) の後ろには絶対にテンを打たない (マルと誤読され
  文が切れ意味が正反対になる)。並列・同格はナカテン (・) で区切りテンの役割を
  侵害させない (`A 級戦犯容疑者・岸信介`)。
- (本多 第四章)

### 係助詞「ハ」は題目を示す (主語ではなく述語が日本語の大黒柱)

- 「ハ」はガ/ヲ/ノを兼務して文の題目を提示する (象は鼻が長い)。題目の「○○ハ」は
  文頭に置く義務はなく、長い文では述語の近くに置いて直結すると読みやすい
  (「突然現われた裸の少年を見て男たちはたいへん驚いた」)。強調したいときだけ
  題目を冒頭に出し逆順テンを打つ。
- (本多 第六章)

### 否定の動詞は限定の「ハ」とセットにする

- 「日本のようにガラス片がない」だと日本にも無いと読める。「日本のようには…ない」
  と限定の「ハ」を否定動詞とセットにする。「完全に回復しなかった」は「完全には
  回復しなかった」。一文 (一節) 中で対照の「ハ」は **三つ以上使わない** (二つまで。
  何と対照しているのか不明になる)。
- (本多 第六章)

### 無色透明の接続助詞「ガ」(逆接でない「そして」程度) を警戒し文を切る

- 逆接でない「ガ」は一切の関係・無関係をつなげてしまう便利語ゆえ曖昧な文を
  量産し、読者は「ガ」で逆接を予期し思考を一瞬乱される。「…というのですが」
  「…と思いますが」を見つけたらそこで文を切る (「…といわれています。」)。
  (中島・塚本「『…が、』は逆接だけに使用」と同趣旨。)
- (本多 第六章 / 中島・塚本 1.3.1)

### 「マデ」と「マデニ」を厳密に区別する (一字で論理が正反対)

- 「来週まで掃除せよ」= 一週間掃除し続ける、「来週までに」= 期限内に一度すれば
  よい。継続動作なら「マデ」、最終期限なら「マデニ」(「○日までにご投函下さい」)。
- (本多 第六章)

### 並列の助詞 (と/や/も/か) は最初の単語に付ける

- 日本語の並列助詞は前の要素と続けて発音される (英語の and とは逆)。「クジラや
  ウシ・ウマ…アザラシは哺乳類だ」(×「クジラ、ウシ…そしてアザラシ」)。「も」「か」は
  全体の最後にも付ける。
- (本多 第六章)

### 改行 (段落) は思想のまとまりの境界にだけ置く

- 段落は思想表現の単位。関節でない所で曲げれば骨が折れる (思想が引き裂かれる)。
  一行段落も一章一段落もありうる。「長くなったから」で改行しない。改行すべきか
  自分で分からないなら論理的に書けていない兆候。
- (本多 第七章)

### 漢字とカナの混合は「わかち書き」の効果を担う

- 同じ形の字が続くと拾い読みになり読みにくい。カナが長く続く所はまず漢字、次に
  傍点・カタカナ、それでもダメなら半角アキ。**漢語の熟語を半分だけカナにしない**
  (「こん虫」「両せい類」は禁。書けない語は言い換え `書翰`→`手紙` かルビ)。
  送り仮名は一人の書き手の中で統一する (混在禁止、誤読のおそれは多めに送る)。
- (本多 第五章)

---

## 🇯🇵 和文の誤用・表記チェック (lint 観点) — 『問題な日本語』系

reviewer・校閲委員が一読して引っかかる和文の誤りを機械検出する観点。
**過剰指摘 (false positive) を避ける** ため、定着した用法は誤りフラグしない
ルールも併記する。これらは将来 `paper_writing_check_misuse_japanese` 系の
lint ルールとして展開する候補。

### 仮名遣い・送り仮名の確実な誤り (機械検出可)

- **「やむをえない/ざるをえない」** を「やむおえない/ざるおえない」と書くのは誤り
  (語源「止む事を得ず」の格助詞「を」)。lint: `(やむ|ざる|余儀)お(え|へ)`。
- **「こんにちは/こんばんは」** を「こんにちわ/こんばんわ」と書くのは誤り
  (副助詞「は」)。
- **「雰囲気 (ふんいき)」** を「ふいんき」と書く/読むのは誤り (音位転倒、変換不可)。
- 動詞「言う」のかな表記は **「ゆう」ではなく「いう」** (現代仮名遣いの明文規定)。
  lint: `(そう|どう|こう)ゆう`, `ゆって`。
- **オ列長音は「う」で書く** が、歴史的に「ほ・を」だった語は「お」: 「通り=とおり」
  (×とうり、変換不可)、「お待ちどおさま」、一方「ぞうっと」(×ぞおっと)。
- **連濁「ぢ・づ」**: 二語連合は「づ」(鼻血=はなぢ、人妻=ひとづま)、現代語意識で
  一語なら本則「ず」(稲妻=いなずま、融通=ゆうずう)。
- (『問題な日本語』各項)

### ら抜き・レ入れ・サ入れ (活用の誤り)

- **ら抜き言葉** (可能) は書き言葉では避ける: 上一段・下一段・カ変の「見れる/
  来れる/食べれる」→「見られる/来られる/食べられる」。ただし五段由来の可能動詞
  (蹴れる・しゃべれる・すべれる) はラ抜きではない (**誤検出注意**)。
- **レ入れ** 「見れれる/食べれれる」は誤用。
- **使役のサ入れ**: 五段・サ変は「せる」が付く。「読まさせていただく」→「読ませて
  いただく」。五段未然形 (ア段)+「させて」を検出。
- 様態「そうだ/すぎる」: 助動詞「ない」は「知らなそうだ/読まなすぎる」が一般的
  (「知らなさそうだ」は俗用)。形容詞「無い」は「なさそうだ」(さ要)。接尾語「げ」も
  助動詞「ない」には付かない (「つまらなげ」は誤り、形容詞「ない」由来の
  「頼りなげ」は可)。
- (『問題な日本語』各項)

### 敬語の誤り

- **二重敬語** (尊敬動詞+れる) は一般に誤用: 「おっしゃられる/いらっしゃられる/
  召し上がられる/お見えになられる」(ただし「れる」が受身・可能なら誤用でない)。
- 「お求めやすい/お使いやすい」は誤り → 「お求めになりやすい」(「お…になる」文型)。
- 「…が来ていただく」は不適 → 「…に来ていただく」(授受の出どころは「に」)。
- 「おる (おります)」は謙譲語、尊敬には「おられる」(「お待ちになっております」は
  不適 → 「…ておられます」)。「あげる」は本来「やる」の謙譲語、目上には「差し
  上げる」(「猫に餌をあげる」は美化語化で許容)。
- (『問題な日本語』各項)

### 語彙の混同・混交

- **混同を避ける**: 「おざなり」(その場限りでいい加減・一応はする) ≠「なおざり」
  (無視して放置・しない)。「耳ざわりのよい」は誤用 (本来「耳障り」= 不快) →
  「聞き心地のよい」。「なにげに」は誤用 →「なにげなく」。
- **混交 (コンタミネーション)**: 「合いの手を打つ」→「合いの手を入れる」、「明るみに
  なる」→「明るみに出る」、四字熟語「侃々諤々/喧々囂々」を混ぜない、「胸先三寸」は
  誤り→「胸三寸」。
- **同音異義の漢字使い分け**: 開放 (開け放つ)/解放 (束縛を解く)、改定 (定め直す)/
  改訂 (書物訂正)、絞める (首)/締める (きつく)、空ける (からに)/開ける (ひらく)、
  受 (受ける)/授 (さずける) →「授賞式」が正 (×受賞式)。
- **年齢に「個」を使わない**: 「二個上」→「二つ上/二歳上」。
- (『問題な日本語』各項)

### 形式名詞・「こと/とき」・「的」

- **形式名詞「こと」はかな**で書く (「言うこと」「…たことがある」)。「事」は事件・
  事態など実質名詞のときだけ漢字。**経験は「…たことがある」**で表す (「…たときが
  ある」を経験の意で使うのは誤り。「とき」は時間の存在)。
- 接尾語「**的**」の和語・口語的乱用を避ける: 「わたし的/気持ち的/暮らし的」は
  不自然 (「的」は本来漢語名詞に付く、漢語+的 = 金銭的は可)。
- (『問題な日本語』各項)

### 口語・話し言葉を改まった文で避ける

- 文頭接続詞の **「なので」** は文章語として未定着 →「そのため/したがって」。
  形容詞+「です」(「ないです」) は口頭語では可だが文章語ではやや落ち着かない →
  「ありません」。「すごいおいしい」は連用形「**すごく**」に。話題転換・断定回避の
  「っていうか」「みたいな」は改まった文では避ける。不要な「(…の) ほう」のぼかし
  (金額など正確さが要る場面) を避ける。接客の「…になります」(変化がない場面) →
  「…でございます」、過去形「よろしかったでしょうか」(現在の確認) →「よろしい
  でしょうか」。
- **「全然」** は否定と呼応させる (改まった文では「全然いい」より否定呼応か
  「まったく/きわめて」)。
- (『問題な日本語』各項)

### 外来語表記

- 日本語音に従う: 「メール」が標準 (「メイル」も誤りではない)。長音「ー」は原則
  残す (短い語で末尾「ー」を省くと不明瞭: メーカ/ユーザ/シャッタ →「ー」付き推奨)。
  「ィ」で半長音を表すのは違反 (「フロッピィ/ファジィ」→「フロッピー/ファジー」)。
  原音への近さより日本人に読みやすい慣例・単純形を優先 (ヴィエトナム×、ベトナム○)。
- (『問題な日本語』各項)

### 過剰指摘を抑制する (定着した用法は誤りフラグしない)

- 「とんでもありません/とんでもございません」は今や許容。重言だが定着した
  「犯罪を犯す/歌を歌う/選挙戦を戦う」(動作の結果生じるものをヲでとる正用) は
  追放しない。「汚名挽回/名誉挽回」は許容 (「挽回」に巻き返しの意)。移動・経過・
  離脱点等を表すヲは自動詞 (「街道を行く/山を登る/学校を卒業する」は他動詞誤用
  ではない)。
- **これらを機械的に誤りフラグしないことが lint の品質**。
- (『問題な日本語』各項)

---

## 📑 和文論文の体裁規律 — 中島・塚本『知的な科学・技術文章の書き方』(補完)

skill.md は §1.3.5 (近接障害)・§1.3.2 (文末変化)・§1.3.1 (接続詞)・§6.4 (校閲対応)
を収録済。未収録の **題目・結論・緒論・記述様式** の規律を追加。

### 結論は箇条書きで 3〜4 項目以上 (平均 3.6、学会賞論文 4.6)

- 非箇条書きの結論は「水増し/結論把握困難」と判定されやすい (285 編中 83% が
  箇条書き)。項目が 2 つだと 1 つが既知と判定された瞬間に掲載否になりうる。
- 「得られたおもな結論は以下のとおりである」+ 箇条書き。重要な結論ほど上位に
  置く (本文出現順でなく。第 1 項目に平凡な結論を置くと後続の重要結論を読まれない)。
  水増しでなく実質を再検討して漏れた成果を探す。
- **残留課題 (今後の課題) は原則避ける** (「さまざまな問題が残されている」は
  完結性・独創性の欠如を疑わせ、論敵型読者には手の内を明かす愚行)。書くなら
  「未解決の欠陥」でなく「次の発展目標」に限る。卒論・修論・紀要は引き継ぎ目的で可。
- (中島・塚本 2.6.2)

### 緒論と結論の同一文章重複を排除する (重複障害)

- 読者は緒論を読んだ後すぐ最終ページの結論をめくる。緒論末と結論冒頭が同文だと
  「二度手間を取らせる」と評価が下がる。第 1 稿後に緒論-結論を並べ通読し、重複
  箇所をマークして結論側を短縮 (「本研究で得られた結論は以下のとおり」一文化が
  最も簡単)。和文摘要も緒論と同文で済ませない (内容を選別し再構成)。
- (中島・塚本 2.6.2, 2.11.1)

### 英文摘要は和文の直訳にせず検索キーワードを必ず含める

- データベースのキーワードは主に英文摘要から抽出される。キーワード選定が
  粗雑だと類似テーマの検索から漏れ成果が埋もれる。専門領域を的確に網羅する
  キーワードを先に選定し、それらを含む英文摘要を構成する。
- (中島・塚本 2.11.2)

### 実験方法は「独創性抽出法」で書く

- 細かい手順を冒頭からダラダラ列記すると緒論・結論しか読まない読者を引き込め
  ない。独創的箇所・具体名称・既手法との差を冒頭に: 執筆事項を列挙→既報告手法を
  除外→独創箇所を抽出→装置/手法に具体名称 (例「瞬間検出装置」) を付け「本装置の
  特徴は…にある」と冒頭で要約・優位性主張。細部手順はその後。
- 実験条件は本文に流し込まず **表 (または箇条書き) で別出し** する (「その他の
  実験条件は表 1 のとおりである」)。実験方法を過去形だけで書かない (追試で再現
  可能だから過去形の必然性はない、全文過去形は単調で迫力に欠ける)。
- (中島・塚本 2.4.2, 1.3.5)

### 他者研究の批判は中立に紹介したうえで控えめに (ケナシ制御法)

- 内容紹介なしに「致命的欠陥」と断じると独断的中傷で最低の優位性しか主張できず、
  第三者から逆に評価を下げられる。引用文献の内容を厳正中立に紹介→「…のため…に
  無理があった」程度の控えめ指摘に。実名引用時は非難でなく「初めて綿密に測定した」
  等の称賛を添える。(Wallwork §10.9-10.10 と同趣旨。)
- (中島・塚本 2.4.2, 2.6.2)

### 図番・表番は文末でなく文頭に置く

- 読者は本文と図を交互に参照する。3 行以上の説明文で図番が文末だと、読み終える
  までどの図の話か分からず参照のタイミングを逃す。「図 6 は、…を示したもので
  ある。」と図番を冒頭に。短文なら文末・末尾の (図 3) 併記も可。
- 図を本文で参照するときは、図番号だけを添えて済ませない。本文側で「何を示す図か」
  「どの構成要素を見るべきか」「次の議論でどの性質を使うか」を少なくとも一つ明示する。
  特に概念図・回路図・模式図では、caption に任せきりにせず、本文の最初の参照文で
  図の役割を説明してから式変形や議論に入る。
- (中島・塚本 1.3.4)

### 和文の構文障害 (読み返して直す)

- **修飾語と被修飾語の物理的距離を短縮** (本多と同趣旨)。**親亀・子亀・孫亀**
  (名詞を動詞で前後から挟む三段重ね) を避ける (「出力電圧に比例する回転数に変換
  するインバータ」→ 動詞を移動し分割)。**() の多用を避ける** (() 内に文章・特に
  複数文を入れない)。**連続する「の」は 2 回まで**。**1 文中に同じ助詞「は」「が」を
  複数回使わない** (名詞を記号に置換して検算)。複雑な分類・関係は **図表化** する
  (専門用語を記号に置換すると論理矛盾が露見)。
- 長文は悪文ではない (短文化は字数増・文脈の整合崩れを招く)。文の長短より、
  文節境界に述部+接続詞を入れ意識に落とし込めているかが基準。慣用句の受け答え
  (「…理由は」→「…ためである」) を順当にする。
- (中島・塚本 1.3.6, 1.3.3, 1.3.5)

### 科学・技術文章の基本表記ルール

- 横書き・「である」口語文章体。主語が執筆者なら省略 (**「私」は主語にできない**、
  文献紹介では「著者」「著者ら」のみ可)。句読点は「,.」または「,。」、テン「、」は
  横書きで使わない、並列名詞は中点「・」(`科学・技術専門家`)。形式名詞・補助動詞・
  接続詞・指定の副詞は平がな、常用漢字外・当て字は使わない。漢数字は固有専門用語・
  慣用語・数詞のみ (「12 ビット」「1996 年 7 月 13 日」はアラビア数字)。
- 本文中で等号を助詞代わりに使わない。`回路の高さ=電池の起電力` のような書き方は
  「定義」「同一視」「対応」のどれかが曖昧になる。 prose では
  「回路の高さは電池の起電力に対応する」「縦の長さが電池の起電力に等しい」のように、
  関係を言葉で明示する。等号は数式・表・定義式に限定する。
- (中島・塚本 基本ルール)

---

## 🛡️ 研究倫理と引用 — 中島・塚本『知的な…』(補完)

skill.md の Reproducibility/Open Science ツール (`reproducibility_open_science_
check`) と相補的な、引用・倫理の実務規律。

### 引用範囲は厳密に明記し、他者の文章はそのまま転記しない (無断引用の禁止)

- 参考文献番号を付けても、実際の引用範囲 (式・図・文章) を超えて自説のように
  見せるのは無断引用・盗作。文章の丸写しは引用文献明記でも盗作と断定されうる。
  転載した式・図・文章はすべて参考文献明記、他者が新定義した術語は「」で引用、
  文章は可能なら自分の言葉で再構築。(Wallwork §11 paraphrase と同趣旨。)
- (中島・塚本 3.1.3)

### 参考文献は公表出版物のみ・否定的引用は著者名を出さない

- 「投稿中」「校閲終了掲載待ち」「将来の後続論文」「私信」は引用しない (未公表・
  入手困難、時系列矛盾)。私信由来のアイデアは本文/脚注に「〜氏との私信で想起した」と
  付記。否定評価で実名を出すのは失礼 (番号のみ)、肯定的引用・式や概念の転載では
  著作権回避のため実名明記。
- (中島・塚本 2.10)

### 脚注・付録は極力付けない・二重投稿禁止・データねつ造禁止

- 脚注 = 持て余し論文、付録 = 水増し論文の印象 (校閲過程の欠陥補修と疑われる)。
  脚注内容は本文に挿入。長い付録は 2 論文に分ける。
- **二重投稿禁止** (1 論文 1 学会主義)。追試論文は原論文を必ず引用 (優先権侵害
  回避)。類似論文を必ず調査・引用する (独創性の防衛、未引用は無断引用と疑われ
  掲載否)。(skill.md NG#6 self-citation, §K self-plagiarism と相補。)
- **データのねつ造禁止** (欠損点の創作だけでなく **外れ値の削除/移動も**)。平滑/
  理論曲線から外れた点もバラツキ自体が意味を持つ。外れ点は削除せず、分散が実験
  手法由来か真の現象かを明確化する。**不作為の行為** (欠点・弱点を隠す、追試
  不能なほどノウハウを伏せた虫食い論文) も避ける。
- (中島・塚本 2.8, 3.1, 3.2, 3.5.2)

### 数式の誤りは別経路で検算する

- 数式が誤ると以降の記述すべてが非論理と判断され掲載否。数式は別の導出方法で
  確認し、理論解析と実験解析・前報との矛盾を残さない。
- (中島・塚本 6.2.2)

### 回答文 (校閲照会への返信) は冷静に・全面妥協も全面拒絶もしない

- 校閲委員は前任委員の記録と回答文を比較する裁判官。全面妥協は「主体性なし」、
  全面拒絶は「頑固」と判断される。(a) 全面妥協しない (b) 指摘された独創性・
  信頼性は **新図/記述で完全補強** (c) 理屈で武装し冷静に。根幹でない正当な指摘・
  枝葉末節の照会は子供じみた反論を避け素直に修正。(skill.md §E/§F Reviewer 対応の
  和文版。)
- (中島・塚本 6.4)

### 投稿前に M 日 (7〜10 日) 寝かせ、第 4〜5 稿まで推敲する

- 完璧と思っても抜かりがある。最後まで残る誤りは執筆者が正しいと思い込んだ箇所で、
  他人/時間を置いた自分にしか発見できない。第 1 稿で論旨確定 (てにをはは後回し) →
  接続詞/文末/禁止事項/弱点補強で推敲 → 数日寝かせて再チェック → 指導教官の添削を
  必ず受ける。(skill.md Phase 0・Wallwork §20.2 と同趣旨。)
- (中島・塚本 2.3.3)

---

## 関連 MCP

- `grant-writing` - 申請書 (科研費 / JSPS / KDDI / パワーアカデミー)
- `presentation` - 学会発表スライド (IEEJ SA / IEEE conference / 社内セミナー)

三者は共通の診断 tool (overfull hbox / sentence length / weak expressions) を持ちつつ、それぞれの文脈で閾値と推奨を変えている。
