# grant-writing

Use this skill for grant proposal drafts, recommendation letters, and final
pre-submission checks.  It is intentionally parallel to the other document
servers:

- `paper-writing`: journal / digest manuscripts
- `figure`: publication-quality figures
- `grant-writing`: proposal logic, prose, feasibility, budget, social impact
- `presentation`: slides and oral-delivery artifacts

## Core Rules

1. Start from the reviewer's question: what public or technical problem is
   solved, why now, and why this applicant can finish it.
2. Convert "we will consider" into verifiable deliverables: a dataset, demo,
   prototype, public repository, report, measurement, workshop, or presentation.
3. Keep social implementation concrete: user, field, PoC, evaluation method,
   and handoff path.
4. Budget text must explain why each cost is necessary for the proposal's
   verification loop, not for generic productivity.  Do not make the requested
   budget artificially small: apply close to the program ceiling when the plan
   genuinely needs it, and make that look natural by showing careful itemized
   calculations.
5. Recommendation letters should say only what the recommender can plausibly
   attest: applicant ability, institutional fit, support, feasibility, and
   significance.
6. Do not turn absence in a small or bounded literature corpus into evidence
   that a whole field has failed to adopt a method or that an academic gap
   exists. Use that survey for background or case selection; state the research
   gap independently as an unresolved scientific condition, theory,
   comparison, or testable hypothesis.
7. Keep named software out of the title, summary, academic question, aims,
   novelty, and impact unless the software itself is the research object. Frame
   these sections at the category level, such as OSS, high-order finite-element
   OSS, or an OSS coupling platform. Use names in methods, preliminary results,
   collaboration evidence, rights, and costs when they make the plan verifiable.
8. Assume reviewers know the academic field, not every OSS/AI acronym, foreign
   institution alias, laboratory shorthand, or named benchmark. Use the
   established Japanese term first, explain a necessary acronym by its role,
   and use a benchmark for verification rather than as engineering significance.
9. A platform, repository, protocol, or AI interface is an enabling method, not
   the academic destination. Complete the chain from a field-specific object
   and measured quantity to a conditional law, changed decision, and
   falsification gate.
10. When inventing an interval, score, index, or composite criterion, separate
    calibration from validation. Define its components, freeze its formula and
    threshold before inspecting held-out data, and state what claim is rejected
    or limited when the gate fails.
11. Internal operation, a public repository, user counts, and external links do
    not establish inter-organizational feasibility. A preliminary pilot should
    name the artifact and organizations, show a bounded independent action and
    observed result, and state what remains unproven.

## Literature Evidence and the Academic Gap

限定した技術報告、論文群、会議録をキーワード検索し、ある手法名を確認
できなかったとしても、それだけでは「国内で普及していない」「実装障壁が
ある」「未解決の研究課題である」とは言えない。免責文を添えても、その直後
に同じ一般化を行えば論理上の飛躍は残る。

文献調査は、背景補助、予備調査、比較候補や実証対象の選定に用いる。科研費
の中心的な学術的空白は、既往研究が解いていない条件、理論、比較可能性、
再現・反証可能性、または検証可能な仮説として独立に記述する。普及実態その
ものを研究上の根拠にする場合には、検索式、対象範囲、選定基準、代替語、
集計方法を再現可能に示す。

Bad: 「調査した5報告でXを確認できなかった。したがって国内ではXが普及せず、
研究障壁となっている。」

Better: 「個々の手法は高度化したが、異なる定式化を接続したときに設計判断の
再現性を保証する条件は明らかでない。文献調査は比較対象の選定に用いる。」

## Abstraction Level for Named Software

科研費のタイトル、概要、研究背景、学術的問い、研究目的、独創性、波及効果
では、個別ソフト名を研究概念そのものにしない。たとえば `NGSolve`、
`ONELAB`、`openCFS` と書く代わりに、その箇所で必要な抽象度に応じて
「OSS」「高次有限要素OSS」「マルチフィジックスOSS」「既存OSS連成基盤」
と書く。審査者が評価すべきなのは特定製品の採否ではなく、一般化可能な学術的
条件、方法、成果だからである。

固有名詞を全面禁止するものではない。研究方法の再現、予備成果、共同実績、
権利・ライセンス、予算積算では、実装を特定する名称が必要になる。その場合は
先に一般カテゴリと研究上の役割を述べ、固有名詞は「用いる一実装」または
「遂行可能性の根拠」として置く。解法名、定式化名、標準ベンチマーク名は科学的
条件を特定するために必要であり、単なるソフト製品名と同一に扱わない。

Bad: 「本研究の問いは、NGSolveとONELABを接続できるかである。」

Better: 「本研究の問いは、内部形式の異なるOSS解析モジュールを結合しても、
設計判断を再現可能に保つ条件は何かである。」

## Reviewer Vocabulary and Benchmark Role

科研費の審査者には当該学術分野の専門性を期待してよいが、ソフトウェア開発や
AI実装の語彙まで共有しているとは仮定しない。外国大学は初出から日本語の正式名
で記す。`OSS` は「オープンソースソフトウェア（OSS）」、`LLM` は「大規模言語
モデル」とし、`MCP` は英語名を展開するだけでなく「AIが利用する知識・実行
インターフェース」など研究上の役割を日本語で説明する。

分野内でも共有範囲が狭い略称や関数空間記号を申請書の顔にしない。たとえば
`MMM` は「磁気モーメント法」、`H(curl)` は概要・目的・意義では「辺要素」と
書き、必要なら方法節で数学的な適合性を補足する。厳密さを失うのではなく、
一般名、役割、詳細記号の順に情報を置く。

`$\eta=0$` のような短い条件式も、分野外の審査者には意味を伝えない。記号を
最初の見出しや目的文に置かず、(1)起きている現象、(2)実現したい工学的効果、
(3)その変化を測る量の名称、(4)記号と条件式、の順に導入する。たとえば
「運動量の異なる粒子を出口で同じ位置へ戻す」と述べ、その運動量差に対する
軌道位置の変化率を分散関数 `$\eta$` と定義してから、出口分散
`$\eta_{\rm out}=0$` と書く。概要・図・成果名では、可能な限り工学的効果を使う。

標準ベンチマークでは解析経路の整合確認、比較条件の校正、再現性確認の基準を
得る。その成果を肯定形で述べた上で、主成果は制約付き実設計課題における設計量、
選択則、意思決定の改善で示す。固有名は方法または予備結果で一度定義し、実設計へ
進む二段構成を審査者に見せる。

Bad: 「本研究の独創性は、TEAM Problem 28を高精度に解くことである。」

Better: 「単純化した公開基準問題は解析経路の初期検証に用いる。工学的有用性は、
制約付き機器設計で設計判断が保存・改善される条件により検証する。」

## Evidence, Equations, and Persuasion Hierarchy

慎重さを示すために、提示した成果を直後の免責文で打ち消さない。予備成果は
「何を確認できたか」を肯定形で限定し、次に主実証で検証する設計量、判断、適用範囲
へ接続する。限界を隠す必要はないが、成果段落の中心を「何を示さないか」にしない。

Bad: 「公開基準問題で一致した。この一致は工学的有用性の根拠にしない。」

Better: 「公開基準問題で誤差を4%に収め、結合前後の整合確認基準を得た。主実証では
制約付き機器設計に移し、設計順位が保存される条件を検証する。」

数式は厳密さの入口であって、説明の代替ではない。審査者に添字を解読させず、次の順で
置く。

1. 何を判断するための式かを平文で述べる。
2. 各記号を物理的意味と計算操作で定義する。
3. 数式を示す。
4. 式の大小、区間、閾値が変える設計判断を直後に述べる。

表示数式だけでなく、`$x=0$`、`$p>2$` のような短いインライン条件にも同じ順序を
適用する。記号を先に出して後から意味を補う構成は避ける。

本文には成立させる主張、必達実証、判定法を残す。否定、除外、失敗条件、代替策は一つの
リスク段落へ集約し、条件付き追加検証は将来展開へ移す。手法名・ソフト名・略語が一文に
並ぶときは、先に手法カテゴリ、比較軸、変わる判断を書き、固有名は方法節で補う。

申請本文に内部メモの識別子を残さない。`実証A/B`、`L1--L4`、`A・Bで`のように
別の凡例を必要とする略号は、磁気浮上設計実証、解法選択実証、単独検証、比較検証等の
科学的内容で呼ぶ。`（年度末: ...）`、`TODO`、`要確認`、`暫定案`は、到達点を述べる
完結した文へ直すか本文から除く。年度計画では「何を終えるか」を言い切り、作業台帳や
予算管理表の語調を持ち込まない。

## From Enabling Technology to Domain Knowledge

科研費でOSS、GitHub、API、MCP、AI、研究基盤を用いる場合、構築・公開・利用者増
だけを成果にしない。概要と研究目的には、対象物または現象、測る設計量・物理量、
そこから得る条件付き知識、変わる設計判断、反証条件をこの順に置く。道具の説明は
その後に置き、問いを実行可能にする役割へ限定する。

Bad: 「MCP付きOSS基盤を構築し、研究室間連携を促進する。」

Better: 「候補間の性能差に応じて低費用解析で順位を確定できる領域と高忠実度解析
へ進む境界を同定する。MCP付きOSS基盤は異なる実装でこの選択則を反証する手段と
する。」

## Version Control and AI Execution Are Different Layers

機械学習等をOSS共同研究の先例にする場合、「GitHubを導入すれば発展する」と短絡
しない。版付き実装、モデル、データ、共通ベンチマーク、変更レビューが、他者による
再現・比較・共同改良を累積可能にした機構を述べる。その上で、対象分野に共通基盤が
ない理由と、先例をそのまま移せない障壁を示す。

異種数値解析ではGitHubとAI向け実行層の役割を分ける。GitHub/Gitは、版、差分、
issue、試験、科学レビュー、採否、来歴を管理する。MCP等のAI向け知識・実行
インターフェースは、能力、支配方程式、仮定、物理量・単位、検証法、固有実行手順を
自己記述し、異なるコードを発見・意味解釈・実行する入口を与える。リポジトリがある
だけで相手資産を実行できない分野固有の理由を書かなければ、MCPの必然性は伝わらない。

各機関がコード、知的財産、内部API、保守責任を維持し、外部には検証可能な能力と
実行窓口だけを提示する連邦型も許容する。成果は共通ソルバや利用者数でなく、異種資産
の結合で初めて得た分野知、別機関による再現・反証、共同改良版と論文で判定する。

## Derived Metrics: Calibration, Holdout, and Failure

申請者独自の「判定区間」「評価指標」「スコア」「指数」を中心方法にする場合、名称
だけでは再現可能な方法にならない。少なくとも、算出式またはアルゴリズム、観測可能
な構成成分、校正ケースの範囲と件数、校正に使わない保留データ、式・係数・閾値を
凍結する時点、合否閾値、不合格時の扱いを記す。校正データと検証データを同じにして
結果に合うまで係数を動かす計画は、検証ではない。

厳密な理論上界でない経験的指標を用いてもよい。その場合は適用範囲を明示し、保留
データで外れた条件を再調整で消さず、反例、適用境界、または高忠実度経路を必須と
する条件として残す。合格した場合に何を判断できるかも、順位、採否、停止条件等で
明記する。

## Cross-Organization Preliminary Evidence

研究室内で学生と教員が利用した実績は内部実行可能性を、公開リポジトリや公式リンク
は公開到達性を示す。しかし、研究室間の共同研究が成立する証拠にはまだ一段足りない。
予備実証では、どの機関のコード、モデル、データ、試験、変更を、別のどの機関が、
どの限定課題で再実行・変更・レビューし、残差、誤差、採否、失敗等の何を観測したか
を書く。さらに、独立した科学レビュー、設計量の判断、別課題への移転等、未達の範囲
を明記して本研究の出発点にする。

共同研究歴や共著論文は体制の根拠になるが、中心実証の予備結果とは区別する。単なる
利用者数、スター数、リンク、リポジトリ数を、科学的な大学間実証の代用にしない。

## Budget Policy

予算はほぼ上限いっぱいで申請してよい。重要なのは、上限近くであることを
遠慮して隠すのではなく、研究計画、PoC、計算資源、AI agent 運用、評価、
発表・社会実装に必要な経費として精査済みに見せることである。

予算欄では、単価 x 数量 x 月数/回数、年度配分、見積根拠、検証ループとの
対応を具体的に書く。読み手が「この計画なら上限いっぱいになっても不思議
ではない」と判断できる粒度まで積算する。

見積額は、公式料金表または機関の見積書へ遡れるようにする。各外部サービスに
ついて、提供者、料金表URL、料金年度・改定日、参照日、税込/税抜、最低購入単位、
有効期限、数量・期間、通貨と為替換算、端数処理を記録する。申請年度の料金が未公表
なら、現行料金による暫定積算であることと再確認時期を明記する。公式単価から計算
した額、価格変動への予備幅、研究上の使用上限を混同しない。サブスクリプションと
従量課金APIのように別契約となる費用は分ける。

総額は費目別・年度別の双方から再計算し、研究種目の直接経費上限、研究期間、間接
経費を含めるか否かと照合する。旅費は開催地、人数、泊数、航空運賃、登録費、宿泊費、
日当・現地交通へ分解する。将来の会議は、開催地・日程が公式公表済みか、開催地のみ
公表か、未公表かを区別し、直近大会または学内旅費規程を根拠に暫定積算する。

最大費目が何かを明記し、それが研究の中心的な実験・検証行為と一致することを
説明する。旅費が最大なら、一般的な学会参加ではなく、誰が何を再実行・変更・
レビューし、どの成果物を持ち帰るかへ分解する。AI、計算資源、サーバ、CI、
クラウドを「その他」にまとめる場合も、サービス、単価、月数、実行量へ分け、
サーバ運営費とCI・保存費等の二重計上を避ける。

AI・クラウド・HPCを併用する計画では、AI推論、常時稼働の低費用回帰、論文前の
短期集中計算を別費目として示し、各層に単位、期間、成果物を割り当てる。

GPU等の機種名を予算化するときは、アクセラレータ名とホストCPUのアーキテクチャを
公式仕様で区別する。製品名から別のCPUアーキテクチャを推測して研究条件を増やさない。
高速化の成功そのものより、異機種間で保存すべき残差、物理量、設計判断と、不一致を
適用境界として扱う手順を示す。

## Useful Tools

- `grant_writing_usage()`
- `grant_writing_health_report(text_or_path, program="generic")`
- `grant_writing_section_presence(text, program="generic")`
- `grant_writing_kddi_digital_check(text)`
- `grant_writing_kddi_power_electronics_focus_check(text)`
- `grant_writing_kaken_oss_platform_check(text)`
- `grant_writing_internal_evidence_to_external_scale_check(text)`
- `grant_writing_domain_outcome_chain_check(text)`
- `grant_writing_derived_metric_validation_check(text)`
- `grant_writing_cross_organization_pilot_check(text)`
- `grant_writing_named_software_abstraction_check(text)`
- `grant_writing_reviewer_vocabulary_check(text)`
- `grant_writing_persuasion_quality_check(text)`
- `grant_writing_literature_gap_evidence_check(text)`
- `grant_writing_collaborative_integration_risk_check(text)`
- `grant_writing_budget_alignment_check(text)`
- `grant_writing_analyze_sentences(text)`
- `grant_writing_count_weak_expressions(text)`
- `grant_writing_lint_bedrock(text)`
- `grant_writing_recommendation_letter_template(program="kddi_digital")`

For KDDI Foundation Digital Innovation / social implementation proposals,
use `program="kddi_digital"` so the report checks social issue, digital use,
PoC, schedule, budget, feasibility, and implementation outcomes.

For the power-electronics-board CAE-AI proposal, also use
`grant_writing_kddi_power_electronics_focus_check(text)`.  It checks that the
main subject is the power-electronics-board circuit / electromagnetic /
thermal CAE-AI environment, while companies below 1000 employees are kept as
the first users and implementation field.  It also warns when commercial CAE
is framed as an adversarial replacement instead of a powerful but hard-to-access
tool category for which AI/MCP provides an entry point.

When a proposal claims an existing pilot, operation, or preliminary result,
also use `grant_writing_internal_evidence_to_external_scale_check(text)`.
Internal success establishes feasibility, but does not by itself establish
external validity. The optional check asks who used the internal result, what
transferable unit leaves the original setting, which route and external actor
receive it, and how independent success is verified. It remains not applicable
for proposals that do not claim prior internal evidence.

For proposals centered on a platform, repository, protocol, OSS, API, or AI
interface, use `grant_writing_domain_outcome_chain_check(text)`. It requires the
proposal to terminate in a field-specific object and measurable quantity, a
conditional knowledge product, a changed decision, and a falsification gate.
It also checks that the enabling technology is explicitly subordinate to the
academic question.

When the proposal introduces a named interval, score, index, or composite
criterion, use `grant_writing_derived_metric_validation_check(text)`. It checks
the operational definition, observable components, bounded calibration set,
pre-test freeze, held-out validation, acceptance threshold, and consequence of
failure. Tuning and validating on the same cases does not establish validity.

When preliminary evidence is used to justify research across organizations,
use `grant_writing_cross_organization_pilot_check(text)`. It looks locally
around the claimed pilot for a cross-boundary actor, transferred artifact,
bounded task, observed outcome, independent action, and stated remaining gap.
Repository publication, links, and internal user counts alone remain L0/L1
evidence rather than a scientific inter-organizational pilot.

When a proposal uses keyword searches, non-detection, or counts from a bounded
set of reports or papers, use
`grant_writing_literature_gap_evidence_check(text)`. It warns when those
observations are promoted into field-wide adoption claims, causal barriers, or
the main academic gap. A scope disclaimer does not cure a conclusion that
still makes the same inference. Keep corpus observations as supporting
evidence or case-selection logic unless the proposal includes a reproducible
review design appropriate to the adoption claim.

For proposals that integrate software, models, solvers, data, or organizations,
also use `grant_writing_collaborative_integration_risk_check(text)`. Keep the
academic question distinct from the named protocol or product. Count the
provider's initial description, validation, and maintenance effort as well as
the user's setup effort, then evaluate total cost against reuse. Position the
proposal against existing standards and frameworks and reuse them where they
already solve the problem. Separate core experiments from optional transfers.
Treat unchanged rankings, failed coupling, and counterexamples as results when
their conditions are identified. For every collaborator, connect prior work,
available assets, assigned responsibility, and readiness. If process logs or
students are involved, analyze tasks and artifacts rather than individual
productivity and obtain an ethics determination. Confirm asset ownership and
maintainers, and provide a public benchmark or reference implementation when a
named private asset is unavailable.

For KAKENHI proposals that establish an OSS research platform in the AI era,
use `program="kaken_oss"` or call
`grant_writing_kaken_oss_platform_check(text)`. It checks that existing
technical reports are treated as knowledge sources, JP-MARs remains the
governing platform, domestic and overseas preliminary collaborations are
named, and executable outputs are portable across enterprise Windows/Linux,
mdx, HPC, GPU generations, and future architectures. Hardware acquisition
must remain subordinate to reproducible specifications, CI, containers,
AI-assisted environment setup, and long-term maintenance. The checker also
requires a why-now argument about AI accelerating duplicated lab-siloed work,
an upstream-first policy that surveys and reuses existing OSS before creating
new code, and scientific quality gates: test provenance, expected values and
tolerances, documented limits, CI, and independent re-execution. A no-warranty
license never substitutes for validation. Radia should appear only as
preliminary evidence that specialized software and methods can be integrated;
JP-MARs remains the governing platform.

For KAKENHI framing, also use
`grant_writing_named_software_abstraction_check(text)`. The integrated
`program="kaken_oss"` health report runs it automatically. It warns when a
named software implementation becomes the concept in the background,
question, aims, novelty, or impact, while allowing names in methods,
preliminary evidence, collaboration records, rights, and budget evidence.

The same KAKENHI health report also runs
`grant_writing_reviewer_vocabulary_check(text)`. It checks first-use Japanese
explanations for OSS/AI terms, Japanese institution names, readable domain
terminology, and whether a named benchmark is confined to verification or
calibration instead of carrying the proposal's engineering significance.

The integrated health report also runs
`grant_writing_persuasion_quality_check(text)`. It flags evidence immediately
cancelled by a disclaimer, displayed equations without a prose purpose or
nearby symbol definitions, missing post-equation interpretation, defensive
paragraphs, optional branches in the core plan, and acronym piles. The target
is a positive claim followed by its bounded verification, not confidence
created by deleting caveats.
