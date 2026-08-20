# grant-writing

Use this skill for grant proposal drafts, recommendation letters, and final
pre-submission checks.  It is intentionally parallel to the other document
servers:

- `paper-writing`: journal / digest manuscripts
- `figure`: publication-quality figures
- `grant-writing`: proposal logic, prose, feasibility, budget, social impact
- `presentation`: slides and oral-delivery artifacts

## Scope: Not a KAKENHI Skill

この技能群は**あらゆる競争的資金の申請書**を対象とする。科研費、民間財団、
JSPS各事業、公益法人、受託研究、いずれにも共通する欠陥がある。今日までに
規則の根拠として使った資料も、科研費だけではない。

| 制度 | 得られた規則 |
|---|---|
| パワーアカデミー（民間財団） | 空虚な動詞の置換、用語統一、未記入プレースホルダ、採択申請書の型 |
| JSPS外国人研究者招へい | 様式不適合による差し戻し（審査前の失敗） |
| 科研費 基盤(B) | 審査項目ごとの実測評点、問いの独自性、国際性 |

本文の節は次の二種に分かれる。**制度名の付かない節は制度非依存**であり、
どの申請書にも適用する。制度名を冠した節（`KAKENHI ...` 等）は、その制度
固有の様式・語彙・運用に限る。他制度へ流用するときは、対応する欄が存在
するかを先に確かめる。

新しい規則を足すときは、まずどちらかを決める。**一つの制度でしか成立
しない規則を、制度名を付けずに置かない。**

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
12. Write for the actual reviewing conditions: three KAKENHI criteria
    (academic importance, method validity, feasibility/environment), up to
    ~100 proposals read per reviewer in about a month, and possible
    monochrome printing. A draft that even a non-specialist reads smoothly,
    with monochrome-safe figures, identifiable publications, an explicit
    human-rights/legal box, and a complete funding-overlap box, wins.
13. State the central question once, and restate it with the same decisive
    nouns. A summary promising a 「境界」 and a body promising a 「条件」 read
    as two questions even when they mean one. Fix the wording, not the reader.
14. A lint score is not a quality score. These checks see mechanical defects;
    they cannot see whether the argument holds. Read the failures, ignore the
    number, and never edit a draft to satisfy a keyword list.

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

## One Question, Stated Once

申請書は中心の問いを概要と本文の二か所に書く。このとき、同じ問いを別の語彙で
言い直すと、審査者には二つの問いに見える。実測例（2027-2030年度 基盤C、
2026-08-19）では、概要が「順位が変わる**境界**を**定量化**し」と書き、本文1.1が
「優劣を覆さない**条件**を**記述・検証**できるか」と書いていた。意味は同じだが、
答えの形を表す名詞（境界／条件）も、行う操作の名詞（定量化／記述）も食い違って
いた。しかも本文側には、概要では曖昧だとして除いたはずの指示語「その差」が残って
いた。

規則は三つである。

1. **答えの形を表す名詞を一語に固定する。** 条件、境界、範囲、領域、基準、指標、
   選択則のうち何を与えるのかを決め、言い直しでも同じ語を使う。
2. **行う操作の名詞も一語に固定する。** 定量化するのか、記述するのか、検証するの
   かを決める。三つ並べると、どれが必達か読めない。
3. **概要と本文で役割を分ける。** 概要は研究全体の要約、本文は問いを導く論証と
   定義である。同じ文を二度置かない。本文の問いは、その節が固有に導入した定義
   （例: 解析モジュール）から導く形にする。

`grant_writing_central_claim_consistency_check` がこの食い違いを検出する。語彙
カバレッジ型の検査はこの欠陥を見つけられない。必要な語はすべて文書のどこかに
存在しており、欠陥は語同士が食い違っている点にあるからである。

Bad: 概要「順位が変わる境界を定量化する」／本文「優劣を覆さない条件を記述・検証
する」

Better: 概要「順位を確定できる条件をどのように定量化するか」／本文「こうして定義
した解析モジュールについて、順位を確定できる条件をどのように定量化できるか」

## Name the Operation, Not the Aspiration

「統合する」「連携する」「活用する」は、結果を約束して操作を書かない。審査者に
残る問いは常に同じで、**どうやって**である。

2026年度パワーアカデミー申請書を経験ある研究代表者が書き直した実測（2026-08-18、
改稿前後の語数）が、この編集が最優先で行われることを示している。

| 語 | 改稿前 | 改稿後 |
|---|---:|---:|
| 統合し | 3 | **0** |
| 双方向に連成 | 0 | **5** |
| 連成 | 2 | 8 |
| 一体的に | 1 | 5 |

置き換えの型は次のとおりである。

Bad: 「研究者三者が有する非線形磁気モデリング、物理ベース等価回路抽出、高周波損失
測定の技術を統合し、回路シミュレータ上で利用可能なモデルを構築する。」

Better: 「両モデルを巻線電流と誘起電圧を介して双方向に連成し、PoL変換回路の損失と
動作波形を回路シミュレータ上で統一的に解析可能な統合モデルを構築する。」

改稿後も「統合モデル」という名詞は残る（4→8件と増えてすらいる）。**問題は名詞では
なく動詞**である。何と何を、どの物理量を介して、どちら向きに渡すのかを書けば、
同じ「統合」でも審査者は再現できる。`grant_writing_vague_claim_verb_check` は、
同一文に操作を表す語（連成、接続、射影、介して、双方向、電圧・電流など）が無い
まま使われた動詞を検出する。

同じ改稿では、用語も一つに統一されている（点負荷変換回路 3→0、PoL変換回路 3→4）。
初出で「点負荷（Point of Load）変換回路」と定義した後は、以後すべて PoL に揃えた。
定義した略語と元の語を混在させない。

## What These Checks Can and Cannot Judge

この診断群は**機械的欠陥**に強く、**論証の欠陥**に弱い。2026-08-19に実申請書で
両者を突き合わせた結果は次のとおりである。

| 欠陥の種類 | 検査 | 人手の通読 |
|---|---|---|
| 逆茂木文、二重ハ、単調な文末、90字超 | 検出（12件の二重ハまで数える） | 読めてしまい見逃す |
| 色だけで区別した図、根拠なき「該当なし」、特定不能な業績 | 検出 | 見落としやすい |
| 中心の問いが二重・語彙不一致 | 専用検査を足すまで沈黙（10.0/10） | 検出 |
| 指示語の指示対象が不明 | 弱い | 検出 |
| 節構成の不整合（並列でないものを並列に置く） | 検出できない | 検出 |

したがって、**高得点は「機械的欠陥がない」ことしか意味しない**。学術的重要性、
方法の妥当性、遂行能力という三つの審査基準そのものは、語の有無では判定できない。

この区別は出力構造そのものになっている。`grant_writing_health_report` は
**findings**（位置を指せる欠陥、重要度付き）と **questions**（著者が判断する
問い、重要度も点数もなし）を分けて返し、`defect_counts` が正直な要約である。
`defect_score` は findings のみから計算し、**機械的欠陥の密度**だけを表す。
研究の良し悪しではないので、これを上げるための編集は無駄である。

検査を追加するときは、まずどちらに属するか決める。判定基準は上の一文である。
2026-08-20 の実例では、様式説明文の検出規則を作ったが、経験ある研究代表者が
6件を削除して13件を意図的に残しており、文面だけでは区別できなかった。
そこで欠陥報告から外し、件数を返す問いへ降格した。**区別できない規則は
検出器にしない。**

運用上の含意は二つある。

1. **点数を上げるための編集をしない。** 語彙リストに合わせて語を足す行為は、
   審査者から見れば無意味な水増しである。検査が挙げた個々の失敗だけを読む。
2. **節に無い根拠を要求されたら、それは検査の適用ミスである。** 研究目的節に
   予算の積算は存在しない。適用対象でない検査は `applicable: False` を返すよう
   になっており、点数にも入らない。それでも要求が出るなら、渡した文書と
   `program` の組合せを疑う。

## Which Evidence May Become a Rule

検査規則の根拠として使える資料と、使えない資料がある。区別を誤ると、雑音を
規則として固定してしまう。

**使える。**

- **経験ある研究者による改稿の差分。** 採否に依存せず、何を直すべきと考えたかが
  そのまま残る。本文書の「操作を書く」「用語を統一する」規則はこれに由来する。
- **公募説明会や事務部門の指摘。** 「人権欄は指摘が非常に多い」のような、実際の
  審査運用に基づく情報。推測で作った軸より確実に効く。
- **審査コメント。** 入手できれば最良である。

**使えない。**

- **無作為抽出や予算枠による不採択。** 内容の評価ではないため、そこから規則を
  作ると雑音を学習する。実例として、ある事業の2回連続不採択は無作為抽出による
  ものであり、申請書の欠陥を示していなかった（2026-08-20 確認）。
- **採否そのもの。** 採択された申請書にも欠陥はあり、不採択の申請書にも優れた
  記述はある。採否は分野、競争率、審査員構成に強く依存する。

したがって、採否を規則の根拠にする前に、**その採否が内容審査の結果かどうかを
確かめる**。確かめられないなら、その事例は規則化しない。改稿差分は採否を待たずに
使えるため、教師データとして最も扱いやすい。

## The Shape of an Adopted Proposal

> 出典は民間財団（パワーアカデミー）の採択2件。制度非依存の型として扱えるが、
> 欄の名称は制度ごとに読み替える。

採択された申請書を分析すると、書き方に再現可能な型がある。以下は
パワーアカデミー2021（個人型・100万円）と2022（チーム型・250万円、
粒子線がん治療用加速器電磁石）の採択申請書から抽出したものである。

**業績欄は情報密度で決まる。** 2022年の申請者欄は実質2ページを占め、
査読付論文10本前後を正式書誌（著者、題名、誌名、巻号頁、年）で列挙し、
受賞・採択歴を年月付きで、外部資金は助成名・テーマ・**金額まで**書いている。
「主要業績あり」と要約せず、審査者が確認できる形で並べる。

**研究背景は5段で進む。**

1. 社会的文脈の数字（「電力消費量の46%」のような実数）
2. 技術的ボトルネックの特定
3. 既往研究の限界を**明示する**
4. 本研究のアプローチ
5. 目的を1文で締める（「〜の確立を目的とする」）

**研究内容は番号付きサブテーマにする。** ①〜⑤に分け、各末尾に担当者を
書き、各サブテーマから図を参照する。担当者は役職名ではなく能力で紹介する。

Bad: 「近畿大学 菅原准教授と共同で行う。」

Better: 「CLN法の提案者の一人であり、CLN法の拡張に関する研究実績を持つ
近畿大学 菅原准教授の協力の元、申請者が行う。」

**産学連携は企業名を挙げる。** 採択例は「東芝エネルギーシステムズ、
日立製作所、住友重機械工業、バリアンメディカルシステムズ、メビオン
メディカルシステムズ、IBA」と具体名を列挙している。「関連企業」では
展開の実在性が伝わらない。

**設備備品は原則0である。** 上記2件とも設備費0で、消耗品と旅費を同程度に
配分し、チーム型では「その他」が最大費目（250万円中150万円）だった。
外部委託、実験補助、測定サービスに充てている。上限に近い額を申請する
場合は、この最大費目の内訳で説得する。

## What Reviewers Actually Marked Down (実測)

> 制度固有: 科研費の審査結果開示。評点表と項目文言は科研費のもの。
> ただし「問いの明確さと独自性」「国際性」「遂行能力」を見る点は他制度でも共通する。

2025年度 基盤研究(B)一般、小区分21010電力工学関連の不採択案件について、
審査結果開示が示した評点である。同一小区分の採択率は69件中18件で26.1%、
5名の審査委員が2段階書面審査を行った。

| 評定要素 | 本課題 | 採択課題平均 | 差 |
|---|---:|---:|---:|
| ①研究課題の学術的重要性 | **1.60** | 2.83 | **-1.23** |
| ②研究方法の妥当性 | 2.60 | 3.08 | -0.48 |
| ③研究遂行能力及び研究環境の適切性 | 2.40 | 3.33 | -0.93 |
| B. 国際性 | **1.60** | 2.70 | **-1.10** |

（評点基準: 4=優れている、3=良好である、2=やや不十分、1=不十分）

「やや不十分」または「不十分」を付けた審査委員の数は、次の項目に集中した。

| 指摘された項目 | 5名中 |
|---|---:|
| 研究課題の核心をなす学術的「**問い**」は明確であり、学術的独自性や創造性が認められるか | **3** |
| 学術的に見て、推進すべき重要な研究課題であるか | 2 |
| 研究方法等は具体的かつ適切であるか。研究経費は研究計画と整合性がとれているか | 2 |
| これまでの研究活動等から見て、研究計画に対する十分な遂行能力を有しているか | 2 |
| 着想に至る経緯や、国内外の研究動向と研究の位置づけは明確であるか | 1 |
| より広い学術、科学技術あるいは社会などへの波及効果が期待できるか | 1 |
| 研究環境（施設・設備・研究資料等）は整っているか | 1 |
| 研究目的を達成するための準備状況は適切であるか | 0 |

**研究経費について「問題がある」と評定した審査委員は一人もいなかった。**

ここから読める運用上の結論は三つある。

1. **最大の失点は「学術的『問い』の明確さと独自性」である。** 5名中3名が
   不十分と評価し、①の平均点を採択水準から1.23下げた。この技能群が
   `grant_writing_central_claim_consistency_check` と Core Rule 13 で
   問いの一貫性を扱うのは、この項目に直接対応するためである。
2. **予算の書き方は失点していない。** 誰も問題視していない。予算検査を
   全文書に適用して HIGH を出す設計が誤りだった（2026-08-20 に適用性
   ゲートで是正）ことを、この所見が裏付けている。
3. **国際性が採択水準から1.10低い。** 令和9年度から基盤(B)は国際性評価の
   高い課題へ重点配分が拡充されるため、基盤(B)に出すならここは配分に直結
   する失点である。

この所見は一件の実測であり、分野・年度・審査員に依存する。規則として
使えるのは、審査項目の文言そのもの（上表の左列）と、どの項目が実際に
低評価を受けたかという事実である。

## International Standing Must Be Shown, Not Declared

国際性を評価する制度は科研費に限らない。JSPS の国際事業、民間財団の国際
枠、受託研究の国際連携要件など、いずれも同じことを見る。そして
「国際的に展開する」と書くだけでは、どの制度でも点にならない。

示せるものは四つある。

1. **相手先の名前。** 「海外の研究者と連携する」ではなく、機関名と個人名を
   書く。審査者が確認できない連携は意図の表明にとどまる。
2. **既にある国際的成果。** 共著論文、国際会議発表、国際レビュー、国際
   ベンチマークへの参加。これから作るものだけを並べない。
3. **往来の双方向性。** 何を渡し、何が返るのか。招へい・派遣・共同実装・
   共著のいずれかで具体化する。
4. **自国発の価値。** 国内で発展した手法があるなら明示する。

四つ目は特に注意が要る。「世界水準に追いつく」という枠組みは、価値が一方向に
流れることを自ら認めている。国際性の評価項目が問うているのは逆で、けん引・
貢献・独自価値の創出である。

Bad: 「世界水準の解析技術に追いつくことを目標とし、将来的に国際展開を目指す。」

Better: 「日本発のCauer縮約・磁気モーメント法と欧州発の高次要素・実効表面
インピーダンス法を相互検証し、双方へ還流する。グラーツ工科大学での共同開発と
IGTE共著がその起点である。」

`grant_writing_international_standing_check` が四点を検査する。国際性に
言及がない申請書には適用しない。機関名を挙げていれば「国際」の語がなくても
適用する。

実績と予定は分けて書く。採択済み・掲載済みのものは状態を明記し、これから
行うものは予定と分かる形にする。両者を混ぜると、予定が実績に読める書き方に
なり、審査者が確認したときに信頼を失う。形成途上のネットワークは、途上である
と正直に書いた方が「本研究で何を進めるのか」に自然につながる。

助成期間の開始前に確定している往来（採択前に決まっている招へい等）は、年度
計画ではなく準備状況に書く。期間外の予定を年度計画へ入れると、計画の実行性が
かえって疑われる。準備状況に置けば、本研究が立ち上がる根拠として働く。



なお2025年度の実測（下記の審査結果）では、この軸が1.60対2.70と最も差の
大きい項目の一つだった。実体があっても書き方で失点する。

## Why This Partner: Answer It in Both Directions

国際連携を書くとき、双方の資金提供者が同じことを問う。**なぜその相手なのか。
近場では代えられないのか。**

2026年の実例がこの問いの形を示した。オーストリアの共同研究者が自国の資金
提供者から「なぜ日本か、欧州にはCERNがある」と問われた。しかし彼が必要と
していたのは施設ではない。自分の表面インピーダンス法を、**独立に発展した
別系統の手法**と突き合わせることであり、その系統は日本で発展し一次文献は
日本語だった。施設と手法系統を取り違えた反論には、その取り違えを指摘するのが
答えになる。

答えは常に三つの要素からなる。

1. **こちら側だけが持つ資産。** 手法名、ライブラリ名、原著者を挙げる。
   「連携する」では、なぜその相手かにも、なぜこちらかにも答えていない。
2. **相手がそれを求めている証拠。** 招請、共同研究の申し出、相手が自国で
   取り上げた事実。こちらの意欲ではなく、**相手の需要**が必然性を示す。
3. **近場で代替できない理由。** 独立に発展した別系統との相互検証である、など。

`grant_writing_collaboration_irreplaceability_check` がこの三点を検査する。
実測では、資産と代替不能性は書けていても**相手の需要だけが抜ける**申請書が
ある。自分の計画を書くうちに、相手から来た動きを書き落とすためである。

Bad: 「ウィーン工科大学と国際共同研究を進め、相互に交流する。」

Better: 「日本発の階層行列ライブラリを用いた積分方程式解法について、ミラノ
工科大学より議論の招請を受けている。独立に発展した別系統との相互検証は、
国内の近隣機関では代替できない。」

なお相手の需要を書くときも、実績と予定の区別は保つ。関心の表明は関心の表明
として書き、合意した共同研究として書かない。

## Writing the Question So Originality Is Recognisable

上表で最も減点された項目は「研究課題の核心をなす学術的『問い』は明確であり、
学術的独自性や創造性が認められるか」である。前半（明確さ）は
`grant_writing_central_claim_consistency_check`、後半（独自性）は
`grant_writing_question_originality_check` が扱う。

後半には三つが揃い、かつ**接続している**必要がある。

1. 問いが述べられていること
2. 何が新しいかを名指しする語（独自、独創、新規、初めて 等）
3. その主張が立つ**既往研究の限界**

3が無い独自性の主張は、根拠のない断言である。逆に、限界だけ書いて問いに
繋がっていなければ、それは他人の申請書の背景である。

限界の書き方は二通りあり、どちらも正しい。

一文で対比する:

> 既往研究は個々の手法の高速化を進めてきた**が**、その差を設計量へ伝播させる
> 方法は確立していない。

二文を接続詞で繋ぐ:

> 既往研究では、個々の手法の高精度化・高速化が進められてきた。**一方**、
> 手法間の差が設計量へ及ぼす影響を定量化する枠組みは体系化されていない。

検査は両形を受理する。一文形だけを要求すると、正しい日本語を不合格にする。

## Measured Against Real Submissions (負の結果を含む)

2026-08-20に、同一研究者の提出済み計画調書10件（採択3・不採択7）で検査を
検証した。結果は二つに分かれた。

**語彙が実物と合っていなかった。** 既往研究の限界を、私は「確立していない」
「明らかでない」という形で想定していた。実物はそう書かない。採択申請書は
能力の否定で書く。

> 従来提案されたコイル形状最適化は、経路情報を**考慮できず**、最適化空間が
> **制限される**。／従来提案されている手法を**そのまま適用できない**。

当初の語彙は10件中0件に一致した。採択申請書にも一致しなかった。「従来」を
単独語として登録していなかったことも同時に判明した（「従来手法」「従来の」
しか持っていなかった）。実物に当てて初めて分かる種類の欠陥である。

**そして、この検査は採否を予測しない。** 語彙を修正した後の実測値は次のとおり。

| | 平均対比文数 | 平均score |
|---|---:|---:|
| 採択 | 2.0 | 8.5 |
| 不採択 | 0.8 | 9.2 |

不採択の方が高い。松尾先生の採択申請書は対比文0で減点される。したがって
`no_gap_against_prior_work` は HIGH ではなく MEDIUM とした。**採択申請書に
重大欠陥を出す規則は、重大欠陥の重みを持てない。**

この検査が測るのは「既往研究との対比を書いたか」であって、「通るか」では
ない。審査結果が示すとおり審査者はこの観点で減点するが、**観点が重要で
あることと、私の代理指標がそれを検出できることは別**である。

なお10件中4件は適用外になった。2018年改革前の様式は「中心の問い」という
語彙を用いないためで、現行様式を対象とする検査としては妥当な挙動である。

## Reviewers Disagree, and the Disclosure Measures It

審査結果開示の「審査委員の数」欄は、単なる件数ではなく**審査員間の不一致の
実測値**である。2025年度の一件では、8項目のうち最も多く減点された項目でも
5名中3名にとどまり、**全員一致の項目は一つもなかった**。

| 項目 | 減点した審査委員 |
|---|---:|
| 学術的「問い」の明確さ・独自性 | 3 / 5 |
| 推進すべき重要な研究課題か | 2 / 5 |
| 研究方法は具体的かつ適切か | 2 / 5 |
| 十分な遂行能力を有しているか | 2 / 5 |
| 着想の経緯・国内外の位置づけ | 1 / 5 |
| 波及効果 | 1 / 5 |
| 研究環境 | 1 / 5 |
| 準備状況 | 0 / 5 |

同一の文書が、同一の項目で「十分」とも「不十分」とも評価されている。

読み方の規則は二つある。

1. **少数意見を欠陥と読まない。** 1/5 の指摘は、その審査員の関心や専門の
   ずれで説明がつく範囲にある。改稿の優先順位を決めるなら 3/5 から着手する。
2. **予算が無傷でも落ちる。** この案件は研究経費について「問題がある」と
   評定した審査委員が一人もいなかったが、不採択だった。逆に、金額記述に
   難のある申請書が採択された実例もある（下記の実測を参照）。予算の巧拙は
   採否を決めない。

ここから、申請書の書き方についての含意が出る。文書は一度評価されるのでは
なく、**複数回独立に評価され、その平均で決まる**。したがって、ある審査員を
強く感心させることより、**どの審査員からも減点されにくくすること**が効く。
機械的欠陥をすべて取り除く作業が意味を持つのは、この一点においてである。
それは採択を保証しないが、**揺れの下側を持ち上げる**。

## No Check Here Predicts Adoption (二度の実測)

2026-08-20に、同一研究者の提出済み申請書10件（採択3・不採択7）で二つの検査を
独立に検証した。どちらも採否を分離しなかった。

| 検査 | 採択の平均 | 不採択の平均 |
|---|---:|---:|
| 問いの独自性 | 8.5 | 9.2 |
| 予算記述 | 7.5 | 7.4 |

一つ目は不採択の方が高く、二つ目はほぼ同一だった。さらに個別に見ると、
採択された申請書が予算記述で10件中最悪（金額重複5件、2.5点）である一方、
別の採択申請書2件は10.0だった。**採択された文書どうしが割れている。**

予算記述の規則は、職業として文書を編集する人物が、後に採択された申請書へ
与えた助言に由来する。良い編集規則であることと、採否を決めることは別で
ある。

この二度の実測が支持するのは、本技能群の設計そのものである。検査は
**位置を指せる欠陥**を報告し、`defect_score` は機械的欠陥の密度だけを表す。
研究の良し悪しでも、採択見込みでもない。点数を上げるための編集は無駄で
あり、上げても通らない。

したがって新しい検査を足すときは、採否との相関を期待しない。期待すべきは
「著者がその指摘を見て、反論せずに直せるか」だけである。それが満たされて
いれば、採否を予測しなくても検査として正しい。

## Reference Proposals Have a Shelf Life

他者の採択申請書を型の参考にするときは、**その申請書がいつの様式か**を先に
確認する。科研費は平成30年度から審査システム改革を実施しており、評定要素も
調書様式も変わっている。手元の採択参考書類7件のうち5件は改革前（2012–2014年）
であり、現行様式の型としてはそのまま使えない。

改革後の採択例は、学術的「問い」という語を明示的に用い、文献引用が桁違いに
多い（1件で39件）という特徴を示した。改革前の例にはこの傾向がない。参考に
するなら、まず年度で分け、改革後のものを優先する。

## Where the Reviewer's Own Words Are

不採択の場合、**科研費電子申請システムで審査結果（所見）を確認できる**。
これは審査者が申請書について書いた唯一の一次資料であり、推測で作った
どの規則よりも価値がある。不採択のたびに必ず取得し、指摘を規則へ落とす。

ただし「Which Evidence May Become a Rule」の区別は保つ。無作為抽出や
予算枠による不採択には所見が付かないか、内容評価を含まない。所見の
有無と内容を確認してから規則化する。

## Program Changes That Affect What You Write (令和9年度公募)

> 制度固有: 科研費のみ。他制度の申請には適用しない。

文部科学省の制度説明（2026年8月）から、申請の書き分けに直接効くもの。

- **基盤研究(B)は国際性の評価が高い課題へ重点配分**が拡充された。基盤Bに
  出すなら、国際共同研究の実体（共著、渡航、相手機関の役割）を書く価値が
  資源配分に直結する。
- **挑戦的研究(萌芽)と基盤研究(C)の重複応募制限が39歳以下で緩和**。対象年齢
  なら併願が使える。萌芽は原則満額支給、基盤系は約7割充足。
- **挑戦的研究(開拓)は総合審査から2段階書面審査へ**変更。書面だけで伝わる
  構成の重要性が増す。開拓では萌芽からの発展性が確認される。
- **学術変革領域研究(B)の年齢要件が45歳以下から49歳以下へ**引き上げ。
- **審査支援AIエージェントが国立情報学研究所と開発中**。生成AIで応募が増え
  審査を圧迫する懸念が背景にある。人が読んで伝わることに加え、機械可読な
  構造（見出し、番号、図参照、正式書誌）を保つ意味が増している。
- 学際性への対応として、複数審査区分を選択できる仕組みが検討されている。

## KAKENHI Review Realities (in-house call briefing)

> 制度固有: 科研費。ただし「白黒印刷」「審査員の読む件数」「業績の特定可能性」は
> 紙で審査する制度に共通する。

科研費の審査基準は3つである: (1)研究課題の学術的重要性、(2)研究方法の
妥当性、(3)研究遂行能力及び研究環境の適切性。各セクションがどの基準で
読まれるかを意識して書く。

審査委員は約1ヶ月の審査期間に、多い場合で100件程度の計画調書を読む。
研究支援部門が毎年全調書を確認してきた実感として、専門外の読者でも
読みやすい調書が圧倒的に採択されやすい。効果的にアンダーライン・太字・
ゴシック体・図表を使い、一読で主張が追える構成にする。

カラーの図や写真は、審査時に白黒印刷される種目がある。色の違いだけで
系列を区別した折れ線グラフは、白黒では何も伝わらない。線種・マーカー・
直接ラベル・濃淡で区別し、白黒でも成立する図を作る(figureサーバの
モノクロ安全則と同じ)。

審査ではresearchmapが研究者番号で参照される。業績を羅列する専用欄は
なくなったが、実行可能性の根拠として調書に業績を書くことはできる。
その場合は業績を特定できる十分な情報(著者、誌名、年など)を添える。
応募前にresearchmapの更新と研究者番号の登録を確認する。

「人権の保護及び法令等の遵守への対応」欄は、例年審査委員からの指摘が
非常に多い。アンケート調査・動物実験・個人情報等を扱うなら、講じる
対策・措置(倫理審査、同意、匿名化等)を具体的に書く。該当がない場合も
「該当なし」と明記した上で、そう判断した根拠を一文添える。

## KAKENHI Program Strategy (充足率と重複応募)

> 制度固有: 科研費のみ。充足率・重複制限は制度ごとに異なる。

基盤研究等の多くの種目では、採択されても申請額の約7割程度に減額されて
内定する(充足率)。上限近くまで申請する方針は維持しつつ、減額後も検証
ループが成立する経費の優先順位を設計しておく。挑戦的研究は原則満額支給
(充足率ほぼ100%)だが採択率が非常に低いため、基盤研究との重複応募を
併用する。39歳以下は挑戦的研究(萌芽)と基盤研究の重複制限が緩和されて
いる。萌芽の最終年度前年度には開拓へ応募できるが、採択時は萌芽の残額を
返還する。開拓の審査は合議審査から2段階書面審査へ移行しており、書面
だけで伝わる調書の重要性がさらに増した。補助金種目は単年度会計(繰越
原則不可)、基金種目は研究期間トータルで執行できるため、年度をまたぐ
計算資源・PoC経費は基金種目で組みやすい。

## KAKENHI Compliance Layer (DMP / OA / 共用 / インテグリティ)

> 制度固有: 科研費。DMP・即時OA・設備共用は他制度では要求されないことが多い。

- 研究インテグリティ: e-Radで所属機関への誓約状況を登録していないと
  応募できない。研究倫理教育の受講修了も応募要件である。
- 応募・受入状況欄: 国内の競争的研究費だけでなく、国外資金、民間財団の
  助成金、企業からの受託研究費・共同研究費も全て記載する。2件目以降は
  本応募課題との相違点と応募する理由を書き、所属組織・役職
  (例: ○○大学教授)を添える。代表課題は分担者を含む金額、分担課題は
  自身の研究経費のみを書く。
- エフォートは、研究専従時間の割合ではなく、教育活動等を含む年間の
  全仕事時間を100%としたときの本研究への配分比率で書く。
- 研究データマネジメント: 原則全種目でDMP(研究データマネジメント
  プラン)の作成が求められ、公開した研究データの情報を実施状況報告書・
  実績報告書で報告する。
- 即時オープンアクセス: 学術雑誌掲載後、即時に機関リポジトリ等へ掲載
  する義務がある。成果公開計画にリポジトリ公開を織り込む。
- 設備共用: 直接経費で購入した研究設備・機器のうち条件を満たすものは、
  検索システムへの登録等により機関内外への共用が求められる。
- 経費細目: 設備備品費は単価10万円以上の物品。英文校閲は人件費・謝金の
  明細に計上する(「その他」ではない)。

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

科研費の基盤系種目では、採択されても申請額の約7割程度に減額されて内定する
場合が多い(充足率)。ほぼ上限で申請する方針は維持しつつ、減額された場合に
どの経費から削るか、削っても検証ループが成立する優先順位を設計時に用意して
おく。挑戦的研究は原則満額支給であり、満額前提の計画を組める。

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
- `grant_writing_kaken_review_format_check(text)`
- `grant_writing_central_claim_consistency_check(text)`
- `grant_writing_vague_claim_verb_check(text)`
- `grant_writing_budget_narrative_check(text)`
- `grant_writing_template_residue_check(text)`
- `grant_writing_question_originality_check(text)`
- `grant_writing_international_standing_check(text)`
- `grant_writing_collaboration_irreplaceability_check(text)`
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

Whenever a draft states its question in more than one place -- which every
KAKENHI proposal does, in the summary and again in the body -- run
`grant_writing_central_claim_consistency_check(text)`. It locates the claim
statements, compares their technical nouns, and reports HIGH when two
statements share a topic but promise different answer-shape nouns
(境界 versus 条件), MEDIUM when the operation nouns diverge (定量化 versus
記述), and LOW when the two are near-verbatim, which wastes the summary. It
is not applicable to a fragment carrying fewer than two claim statements.
The integrated health report runs it for every program.

Budget guidance is judged only where budget content exists.
`grant_writing_budget_alignment_check` reports `applicable: False` for a
research-plan or feasibility section, and the health report leaves it out of
the average, because such a section carries no itemization by design. Where
it does apply, a resource keyword counts only when a money token sits in the
same sentence: 評価 and AI appear throughout ordinary methods prose without
anything being costed.

For any KAKENHI draft (and most other Japanese proposals), run
`grant_writing_kaken_review_format_check(text)`. It encodes the in-house
KAKENHI call briefing: color-only figure discrimination (some categories are
reviewed as monochrome prints), missing safeguards when surveys, animal
experiments, or personal data appear, a bare 「該当なし」 without a stated
rationale in the human-rights/legal box (the box reviewers flag most often),
publication mentions that cannot be identified in the researchmap era, and an
incomplete funding-overlap box (相違点・応募理由・所属組織役職). For full
drafts it also checks coverage of the three review criteria and the presence
of emphasis and figure references, because reviewers read up to ~100
proposals in about a month. The result carries `briefing_notes` with the
review-reality reminders. The integrated health report runs this check for
every program.
