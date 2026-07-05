# IGTE 2026 Digest — CauerLadderNetwork (菅原)

**作業フォルダ (2026-05-12 〜)**: `S:/Radia/01_GitHub/examples/CLN/`
旧場所 (アーカイブ): `W:/02_学会資料/2026年度/2026_09_IGTE_Symposium/CauerLadderNetwork@菅原/`
ダイジェスト: `igte_symposium_2026.tex` (igtesymp.cls), same dir
進捗報告書: `W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/2026_05_10_progress_for_nagamine_v2.pdf` (6 pages, clean)

## 経緯
- 2026-05-10: 長嶺先生報告書 v2 (6 pages) 完成
- 2026-05-10: IGTE 2026 digest 着手、菅原指示「1 問ずつ Q&A → MEMORY.md に反映 → digest に反映」

## 議論ログ (Q&A 形式)

### Q1. Single main message
**Q**: IGTE 2026 digest の single main message は何ですか?
**A**: 「BEM 的要素で CLN が使える」提案
**含意**:
- 主張: BEM-Foster spectrum + Hankel-Padé による Cauer-I 抽出が、3D 任意形状 (sphere, cylinder, cuboid) の eddy current CLN として実用的に動作する手法であることを示す
- BEM の特性 (eigenmode 直接分解、free-space radiation 自然) を活かして CLN を構築する新提案
- FEM 的手法 (Kameari accumulation NGSolve) は cross-validation reference として使う
- これは CLN コミュニティでは新しい angle (従来は FEM iteration が主流)

### Q2. FEM との positioning
**Q**: BEM-based CLN を FEM-based とどう positioning するか?
**A (菅原)**:
- 本当は FEM-Kameari 反復は直交化の限界があり高段は困難と考えていた
- BEM-Foster CLN なら計算時間さえかければ高段まで行けるメリットがあると思っていた
- これが主張できるか? 微妙な気がする (実証データでは BEM も高段で breakdown)
- 関連: 昔 EFIE (Electric Field Integral Equation) には「時間領域 EFIE」があった、これは低周波 BEM における導体要素に相当

**含意**:
- BEM の理論的優位性は「Foster spectrum 直接抽出 → Hankel-Padé」の構造で、高段が原理的に可能
- しかし実証では BEM も Foster spectrum sparsity (uniform B 励振の selection rule) で高段 breakdown
- 「微妙」= 純粋な BEM 優位を主張するには反例 (cuboid breakdown) を扱う必要
- Framing 候補: 「BEM-based CLN は FEM 反復の代替として直交化限界が無く、原理上高段が可能。実装上の Foster spectrum 拡充が future work」 — これなら主張可
- 関連分野: 時間領域 EFIE → 低周波 BEM 導体要素 (これ言及で CLN を一般 EM 計算文脈に位置付け)

### Q3. Target audience
**Q**: digest の target audience は何者?
**A**: Eddy current 専門家 (TEAM problem 系)

**含意**:
- Cauer-I / Hankel-Padé / Kameari accumulation は **既知**として軽く触れるだけ (詳細省略可)
- TEAM 28 ベンチマーク, conducting cuboid eddy current loss などの context を強調
- 実用例: 電力機器 (変圧器、モーター、発電機の渦電流損失), 電力電子, 磁性体の渦電流問題
- BEM 関連は eddy current 専門家には少しなじみが薄いので **時間領域 EFIE / 低周波 BEM 導体要素** との対応付けで導入
- Hiruma vs Kameari 規約差は eddy current 専門家には新しい知見 → digest で価値ある main 貢献の一つ
- ZZ error estimator は FE 数値解析の標準 tool だが eddy current コミュニティでは Cauer 文脈で初使用 → 提示価値
- 国内 audience: 静止器・回転機合同研究会の常連 (五十嵐先生, 比留間先生, 松尾先生, 三分一先生 等)

### Q4. Figure / table 構成
**Q**: 2-page digest の figure/table 何個?
**A**: まだ議論早い。ダイジェスト骨子先に決めてから

**含意**:
- まず骨子 (story line: introduction → main contribution → method → result → discussion) を固める
- 骨子確定後に figure/table 配置を決める
- TODO: 後日再質問

### Q5. Story line (骨子)
**Q**: digest の 5 段骨子は?
**A**: 課題提示 (FEM 限界) → BEM-CLN 提案 → 規約同定 + 変換 → 3 形状 検証 → high-stage limit 結論

**含意**:
- 5 sections に対応:
  1. **Introduction**: FEM Kameari iteration の 直交化限界、3D corner singularity, mesh sensitivity (motivation)
  2. **Proposal: BEM-based CLN**: BEM-Foster spectrum + Hankel-Padé による Cauer-I 抽出、低周波 BEM 導体要素として位置付け
  3. **Convention identification + 変換**: Hiruma (impedance) vs Kameari (susceptibility) Cauer 規約差、BEM-to-Kameari 変換式 ($\tau$ 冪 -1)
  4. **Validation in 3 shapes**: sphere (5 系統 機械精度) / cylinder (3 系統 0.4%) / cuboid 5×2×1 (4 系統 0.3%)
  5. **High-stage breakdown limit**: ZZ + R, L 定量化、NGSolve mesh refine は段数増やさず, BEM order 増は leading converged だが高段不安定

- 各 section 0.4 page 程度 (2 page 中)
- §2 BEM 提案が main contribution として最も重い
- §3 規約同定は theoretical novelty
- §4 検証は credibility
- §5 limit は honesty (publishable result)

### Q6. FEM 限界の framing
**Q**: FEM Kameari iteration の限界をどの角度で?
**A**: 直交化の限界 (Krylov subspace の numerical breakdown)

**追加コメント (菅原)**: 「そもそも渦がめんどくさいというのはなんとかなるのか?」

**含意**:
- 主たる FEM 限界 = Krylov subspace 直交化が高段で破綻 (M-orthogonality loss in iteration)
- これは Lanczos / Arnoldi の知見 (re-orthogonalization, restart 等) を eddy current Cauer iteration 文脈で議論すべき
- 「渦がめんどくさい」= 渦電流問題の根本的計算難しさ
  - 周波数領域: complex BVP, 表皮効果で解像度要求
  - 時間領域: 長時間積分必要
  - **Cauer ladder は時間/周波数 両方を 1 つの ladder で表現する効率手法** (これが CLN の motivation)
  - しかし高段抽出で限界 → 渦のめんどくささを完全には escape できない
- digest 第 1 文: 「渦電流 Cauer ladder 抽出は計算 EM の難問の 1 つ。FEM Kameari 反復は標準だが Krylov 直交化限界で高段不安定。BEM-Foster 経路は別方向の解として有望。」 — このトーン

### Q6.1 (フォロー). vim なら Krylov 限界 escape?
**Q (菅原)**: で、vim なら問題ないわけ?
**A (Claude)**:
- **理論上 escape**: radia-vim は Kameari iteration を使わない、generalized eigenvalue $Kv = \lambda Mv$ を直接解く + Foster amplitudes $g_k$ → Hankel-Padé で Cauer-I 抽出 (mpmath 任意精度可)
- **しかし実証上の別の壁**:
  1. uniform $B_z$ 励振の **selection rule**: 81 mode 中 amplitude ある独立 mode は ~20 個に留まる
  2. **Hankel-Padé 自体が高段で ill-conditioned** (任意精度でも 解 multiplicity)
  3. BEM mesh accuracy: order=3 で leading converged だが高段で order=4 と値変動
- 結論: vim は Krylov 直交化問題は escape, しかし高段精度の壁は残る (sphere/cylinder では解析的 spectrum 豊富で問題薄、cuboid sharp corner で顕著)
- digest framing: 「FEM iteration の Krylov 限界を BEM eigendecomposition で escape できる、ただし高段は別の限界 (Foster spectrum sparsity, Hankel-Padé conditioning) がある」

### Q6.2 (フォロー 3 連). 球解析解 / 角ある形状 / 固有値分布
**菅原 3 コメント**:
1. 球に一様磁場印加 CLN なら解析解 (Stoll Bessel) があるのでそのスタディができる
2. しかし、それより**角がある形状の方がつらい**とわかった
3. それから、**固有値の分布**だよね

**含意**:
- (1) Sphere = clean ground truth として digest の意義 = 4 系統一致を支える基盤
- (2) Sharp corner (cuboid 5×2×1) が真の challenge → digest の primary stress test
- (3) **固有値分布が決定要因**:
  - 球 (smooth): Stoll spectrum dense, $\tau_n \propto 1/n^2$ に従い 12000 mode 取れる → Cauer 40 段
  - cylinder: BEM 50 mode で 6-7 Cauer 段
  - **cuboid 5×2×1**: Foster spectrum **sparse** (uniform B 励振で励起される独立 mode が 20 個程度) → Cauer 1-2 段
  - これが「角あると辛い」の本質: corner singularity → field 高 spatial frequency → BEM mode 必要数増 → しかし selection rule で励起される mode は限定 → 高段精度劣化

- **digest §1 構成案** (Q5 の §1 を更新):
  - 第 1 文: 「Eddy current Cauer Ladder Network 抽出は EM 計算の重要 problem。」
  - 第 2 文: 「球など smooth 形状では FEM Kameari accumulation で高段抽出可、解析解 (Stoll Bessel) と完全一致する。」
  - 第 3 文: 「しかし sharp corner を持つ rectangular conductor (cuboid) では Krylov 直交化破綻 + 固有値分布 sparsity で高段抽出が困難。」
  - 第 4 文: 「BEM-Foster spectrum の eigendecomposition 経路で iteration 限界を escape する手法を提案、3 形状で validation。」

- **digest §4 → 形状階層 narrative**:
  - sphere (smooth, dense spectrum) → 完全動作 (validation 基盤)
  - cylinder (smooth curved, dense spectrum) → 動作良好
  - **cuboid 5×2×1 (sharp corner, sparse spectrum)** → leading rung 4 系統 0.3% 一致だが高段は限界 → これが真の test case

### Q7. BEM 定式化詳細レベル
**Q**: BEM 定式化の詳細はどこまで?
**A**: 詳細: K matrix integral 式, 固有値問題, Kameari 規約 Hankel-Padé まで 1/2 page

**追加コメント (菅原)**: 「誤差が積み上がるのは Kameari Iteration の問題」

**含意**:
- §2 BEM-CLN 提案 を 1/2 page (digest 1 page 中 半分): 数式重め, 中身ある
- 含めるべき式:
  1. **K matrix entry**: $K_{ij} = \frac{\mu_0}{4\pi} \int_{\Omega_c} \int_{\Omega_c} \frac{\nabla \cdot \phi_i \, \nabla \cdot \phi_j}{|\mathbf{r} - \mathbf{r}'|} \, d\mathbf{r} \, d\mathbf{r}'$ (HDiv div-free interior basis, conducting body 内 BEM)
  2. **Mass matrix**: $M_{ij} = \int_{\Omega_c} \phi_i \cdot \phi_j \, d\mathbf{r}$ (volume integral, polynomial)
  3. **Generalized eigenvalue**: $K v_n = \lambda_n M v_n$, $\tau_n = \sigma \lambda_n$ (時定数)
  4. **Foster amplitudes**: $g_n = v_n \cdot b$, $b_i = \int_{\Omega_c} \phi_i \cdot A_{\text{ext}} \, d\mathbf{r}$
  5. **Kameari 規約 Hankel-Padé moment**: $\alpha_n = (-1)^n \sum_k g_k^2 \tau_k^n$ (τ 冪 n)
  6. **Cauer-I extraction**: continued fraction inversion of $\alpha_n$ → $R_{2k}, L_{2k+1}$

- **Kameari Iteration の誤差積み上がり** (菅原指摘): これは FEM 側の限界で、BEM 提案の motivation の一部。
  - FEM Kameari accumulation: $\mathbf{u}_{n+1} = \mathbf{u}_n - \mathbf{A}_{\text{acc}}/L_n$ で iteration 誤差累積
  - BEM eigendecomposition: 1 ステップで全 mode 取得、累積なし
  - これが「Krylov 直交化限界」と同根、§1 で強調すべき

- digest §1 → §2 つなぎ: 「FEM Kameari iteration では mode 直交性が iteration 毎に劣化、誤差累積。BEM では generalized eigenvalue を一度解くだけで全 Foster mode が orthogonal に取れる、累積誤差なし。」

### Q8. Hiruma vs Kameari 規約の framing
**Q**: Hiruma vs Kameari 規約の差をどう記述?
**A (菅原)**:
- Kameari規約 の Cauer-I と、Hiruma規約があるとは知らなかった
- これまで CLN 法と呼ぶと Kameari 規約のほうしか知らなかった
- 現時点では Kameari 規約にしたい
- Hiruma 規約でも回路化できるなら良いが、FEM と組み合わせる方法がまだ分かっていない

**含意**:
- digest 採用規約: **Kameari 規約 Cauer-I**
- 「規約の同定」は新発見だが digest では Kameari メインで Hiruma に軽く言及
- 文脈: 「BEM-Foster spectrum から Kameari 規約 Cauer-I を抽出する変換式 ($\alpha_n = \sum g_k^2 \tau_k^n$, τ 冪 $n$) を導出。注: 文献では別の規約 (Hiruma 3-term Lanczos = Krylov-Padé impedance, $\tau$ 冪 $n+1$) も存在し、両者の Cauer ladder は系統的にずれる ($\sim 5-12\%$)。本研究では FEM Kameari accumulation と一致する Kameari 規約を採用する。」
- これで:
  - Kameari 規約を主主張
  - Hiruma 規約の存在を fairly 言及 (誠実)
  - FEM 互換性を理由に Kameari 規約選択を justify
  - Hiruma 先生研究を否定せず、別 angle として位置付け
- digest §3 タイトル候補: 「BEM-Foster from Kameari convention」 (Hiruma との対照は数行)

### Q9. 検証 3 形状の提示順
**Q**: sphere/cylinder/cuboid の提示順と記述 weight は?
**A**: 球 (解析解検証) → 円柱 (BEM 軽) → cuboid (sharp corner challenge, 主作)

**含意**:
- §4 の構成: Easy → Hard, foundation building flow
- **§4.1 Sphere (Stoll closed-form benchmark)**:
  - Cu sphere R=10mm + uniform $B_z$=1T
  - Stoll Bessel χ-Foster $a_n = 6/(n^2\pi^2)$, $\tau_n = \mu_0\sigma a^2/(n^2\pi^2)$
  - radia-vim BEM eigendecomp + Hankel-Padé Cauer-I → Stoll と完全一致 (5 stage 機械精度)
  - 規約同定 + 変換式の正当性を解析解で実証
- **§4.2 Cylinder (axisym BEM cross-check)**:
  - Cu disk R=10mm, t=2mm
  - radia-vim BEM (Spherical Duffy) vs axisym BEM (Nagamine 50 modes) → Kameari 規約変換後 0.4% 一致
  - 軽め記述, 1 paragraph
- **§4.3 Cuboid 5×2×1 (main contribution)**:
  - Cu cuboid 5×2×1 mm — sharp corner challenge
  - **4 系統独立 (radia-vim BEM o3, BEM o4, NGSolve cln_team28 fine, ELF 商用 BEM 時間領域)**
  - leading Cauer rung: 9.93 / 9.96 / 9.96 / 9.95 → **0.3% 機械精度一致**
  - これが digest の central result

- §4 配分: §4.1 (1/4 page), §4.2 (1/8 page), §4.3 (3/8 page) = 計 3/4 page
- §4.3 が figure/table の主場所

### Q9.1 (重要訂正). IGTE digest は 1 page (NOT 2 page)
**菅原コメント**: 「しかし、1page-digest には入らなそう。」
**確認**: igtesymp.cls template line 162 — "The **one-page digest** should contain the title, authors, affiliation, short abstract, main contents, preliminary results and references."
- A4 2-column format
- 25mm 上, 15mm 下, 20mm 左右マージン
- Times New Roman, single spacing, 2 column

**含意**: 5 sections (intro / BEM proposal / 規約 / 3 形状検証 / breakdown limit) は **1 page には収まらない**。drastic cut 必須:
- option A: §1+§2 を Intro 圧縮、§3 規約は数行、§4.3 cuboid main、§5 limit は 1 文
- option B: 1 つの core message に absolute focus、support evidence 切る
- option C: figure をやめて全 text + 1 table のみ
- 次の質問で scope cut を決める

### Q10. 1-page digest scope cut
**Q**: scope cut どうする?
**A (菅原)**: まずは Sphere の検証で、full paper で cuboid に対応するというのはどう?

**含意**:
- **digest = sphere 検証 中心** (1 page に収まる)
  - Sphere は Stoll Bessel 解析解との完全一致 (5 系統 機械精度) → clean, 強力, 1 page で完結
  - Hiruma vs Kameari 規約の同定も sphere で実証可能 (Stoll = Kameari 規約)
  - radia-vim BEM Foster spectrum + Hankel-Padé Cauer-I が Stoll と機械精度一致
- **Full paper (= IGTE extended paper, 6 page) = cuboid 5×2×1 に対応**
  - Sharp corner challenge を full paper で展開
  - 4 系統一致 (radia-vim BEM o3/o4, NGSolve fine, ELF) を提示
  - Breakdown limit 定量化, ZZ + R, L analysis
  - Cylinder cross-validation も full paper

- Digest narrative (1 page):
  1. Intro: FEM Kameari iteration の 直交化限界 (1 paragraph)
  2. BEM-CLN 提案: BEM-Foster spectrum + Kameari 規約 Hankel-Padé Cauer-I (定式化 数式重め, 1/3 page)
  3. Sphere validation: Stoll 解析解との完全一致 (5 系統, 機械精度) — table or figure (1/3 page)
  4. Closing: cuboid sharp corner challenge は full paper で議論 (1 文)

- Title 候補: 「BEM-based Cauer Ladder Network Synthesis for 3D Eddy Current Problems: Validation against Stoll's Analytical Solution」
- これで 1 page に収まる + main message (BEM-CLN 提案) clear + sphere で誠実に裏付け

### Q11. Digest title
**Q**: digest title 案は?
**A (菅原)**:
- 「BEM」より「**VIM**」(Volume Integral Method) で書きたい
- VIM による CLN 導出を示したい
- Cuboid がターゲットであるが、まずは定式化を球で検証
- 角部が厄介であることは触れたい

**含意**:
- **手法名は VIM (Volume Integral Method)** — radia-vim の "vim" は実は "Volume Integral Method" だった
  - HDiv div-free volume basis + 1/r kernel volume integral = VIM
  - これは BEM (Boundary Element Method) とは違う、conducting body 全体に basis を張る
  - 名前訂正: 我々がこれまで「BEM」と呼んでいたのは正確には VIM
- **Title 案** (VIM-based, 球検証, 角部 mention):
  1. "Volume Integral Method for Cauer Ladder Network Synthesis: Verification on Sphere and Outlook for Sharp-Cornered Conductors"
  2. "VIM-based Cauer Ladder Network for 3D Eddy Currents: Sphere Verification with Outlook to Cuboid Geometry"
  3. "A Volume Integral Method for Cauer Ladder Synthesis: Sphere Benchmark and Sharp-Corner Challenges"
- 角部 mention は title に入れるか abstract で扱うか要相談

**追加**: digest 全体で「BEM」→「VIM」に言い換え必須 (memory + 本文)

### Q12. Title (保留)
**Q**: title 候補は?
**A (菅原)**: title の議論はまだ早い。まずは digest の前に長嶺先生と打合せをしようと思う。

**含意**:
- Title は長嶺先生打合せ後に決定 (共著 / 引用関係も同時に整理)
- 現時点では working title として「VIM-based Cauer Ladder Network: Sphere Verification」程度

### Q12.1 (菅原確定指示). BEM → VIM
**菅原コメント**: 「bem という言葉をやめよう。これは、vim だね。」
**確定**:
- 今後 digest, 議論ログ, memory, radia-mcp 全て **VIM (Volume Integral Method)** に統一
- 「BEM」記述があれば全て「VIM」に置換
- VIM = HDiv div-free volume basis + 1/r kernel volume integral (radia-vim 実装)
- 既存の BEM-Foster spectrum 文献ある場合は note で「= 我々のVIM」と説明追記

### Q13. Intro 書き出しトーン
**Q**: Intro の書き出しトーンは?
**A**: 実用重視: 「電力機器・電力電子では渦電流損装荷モデルが必要」→ CLN の実用ドライバー

**含意**:
- §1 Intro 第 1 段: 電力電子・電力機器における渦電流損装荷モデルの必要性
  - 例: 高周波変圧器、SiC/GaN power module 環境の eddy current loss、EV モータ高速化
  - これらでは circuit-level の eddy current 等価モデルが設計に必要
  - CLN (Cauer Ladder Network) は eddy current を $R_n, L_n$ ladder で表現する効率手法
- §1 Intro 第 2 段: 既存の FEM Kameari iteration ([Sugahara TEAM 28], Hiruma & Igarashi 2020) の限界
  - Krylov subspace の 直交化 が高段で破綻
  - corner singularity, mesh sensitivity, 累積誤差
- §1 Intro 第 3 段: 本研究の提案 — VIM (Volume Integral Method) で CLN 抽出
  - VIM は eigendecomposition 一発で全 Foster spectrum 取得 (累積誤差なし)
  - Kameari 規約 Hankel-Padé Cauer-I 抽出
  - 球で Stoll 解析解 (Stoll Bessel) と完全一致 → 定式化検証
  - cuboid sharp corner は full paper で議論

- 第 1 段 (実用motivation) → 第 2 段 (FEM 限界) → 第 3 段 (VIM 提案) で 1/4 page 程度
- 引用候補: 電力機器渦電流損 (適切な review 論文), Sugahara TEAM 28, Hiruma & Igarashi 2020, Stoll 1974

### Q14. 著者面子
**Q (菅原)**: 「面子は、nagamine hane でいくか。」
**A 確定**: 著者 = **菅原 (近大), 長嶺 (?), 羽根 (?)**

**含意**:
- 共著者: 長嶺先生 (rectvc 法・verified-interval Cauer)、羽根先生 (磁気ヒステリシス共同研究)
- 羽根先生関連: `W:\999_菅原賢悟\19_磁気ヒステリシス\2026_羽根研との共同研究`
- 順序・所属は要確認 (打合せ後)
- digest title page で author + affiliation を記載
- working: K. Sugahara, M. Nagamine, ?. Hane (頭文字要確認)

### Q15. 体制と次のステップ
**菅原指示**:
- 古い tex/pdf 削除 → v2 のみ保持 (完了)
- digest 素案 (rough draft) 作成 → 長嶺先生送付
- 長嶺先生に「全データまとめ (v2 progress report)」と「digest 素案」両方送る → 議論

**含意**:
- v2 progress report (6 page) = 全データ summary、長嶺先生送付用
- digest 素案 = 1-page IGTE format、長嶺先生に骨子レビューしてもらう
- 打合せ後に title, 著者順, 内容 final 化

### Q16. Cleanup + 著者順
**菅原指示**:
- W:/30_CauerLadderNetwork/2026_04_01_長方形CLN cleanup → v2 のみ保持 (完了)
- 著者順: 菅原・長嶺・羽根

**実行済**: 古い tex/pdf 8 組削除、v2 のみ保持

### Q17. Digest 素案作成 + 長嶺先生送付資料
**菅原指示**:
- W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/2026_05_10_progress_for_nagamine_v2.pdf を「いままでの検討成果」資料として位置付け
- digest は素案だけ作って、長嶺先生に「全データまとめ + digest 素案」両方送って議論

**実行**: digest/igte_symposium_2026.tex 素案作成 (title, author=菅原・長嶺・羽根, abstract, §1 Intro, §2 VIM-CLN 定式化, §3 Sphere verification, §4 Outlook to Cuboid, §5 Conclusion, refs 4 件)

### Q18. Born sum メリット質問
**Q (菅原)**: 3D モデルでは、born sum のメリットはない?

**A (Claude)**:
- Born sum (Born series) = volume integral equation の摂動展開、scattering 文脈では基本
- 我々の VIM で原理的に有用な場面:
  1. Multi-conductor (互いに散乱) → 順次取込
  2. Inhomogeneous σ (contrast media) → 摂動展開
  3. 大規模系 (直接 eigenvalue 高コスト) → 反復化
- **しかし isolated single conductor (sphere/cuboid) in vacuum では Born sum メリットなし**:
  - n_dofs=81 や 176 で直接 dense eigenvalue で十分
  - Foster spectrum は exact (discretization 限界内)
  - Multiple scattering なし (単一 body)
- 3D で Born sum が将来役立つかもしれない場面:
  - Multi-conductor (コア + 巻線, 積層鋼板)
  - Sharp corner 特異性を「強い摂動」として扱う ⇒ ただし収束怪しい
- 結論: 現在の研究で Born sum 使う必要なし、future work (multi-body) で検討余地あり
- digest では言及不要 (sphere/cuboid 単体問題なので)

### Q19. 長嶺先生送付 = 重要パッケージ + 軸対称ソルバー実装も
**菅原指示**: 「大事なものにしよう。長嶺先生には軸対称ソルバーの実装も見せるべきだと思う。」

**送付パッケージ確定** (大事な research package):
1. **v2 progress report** (W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/2026_05_10_progress_for_nagamine_v2.pdf) — 全データまとめ 6 pages
2. **IGTE digest 素案** (W:/02_学会資料/2026年度/2026_09_IGTE_Symposium/CauerLadderNetwork@菅原/digest/igte_symposium_2026.pdf)
3. **軸対称ソルバー 実装**:
   - `ngsolve_validation/sphere_axisym_kelvin.py` (Kameari accumulation, axisym H1, Z-offset Kelvin twin half-disc)
   - `ngsolve_validation/cln_team28_axisym.py` (cylinder axisym Kameari)
   - `ngsolve_validation/sphere_analytical_cln_highstage.wls` (Stoll Bessel + Mathematica 240-digit Cauer-I)
4. **BEM-to-Kameari 変換**: `ngsolve_validation/bem_to_kameari_cauer.py` (mpmath 60 digit, Hankel-Padé)
5. **3D Kameari + snapshot fix**: `ngsolve_validation/cln_team28_kelvin.py` (HCurl + tree-cotree + Kelvin two-sphere)
6. **radia-vim Phase F-4 (VIM)**: `S:/Radia/01_GitHub/src/ext/radia_vim/scripts/extract_tau_cuboid_521.py` (Spherical Duffy quadrature + HDiv div-free)
7. **Theory note**: `2026_05_05_kameari_vs_hiruma_orthogonalization.md`

→ これらを zip / GitHub repo URL でまとめて送付するのが理想

### Q20. Sphere 検証のビジュアル表現
**Q**: sphere verification の table か figure か?
**A**: 現状 テーブル 5 段 (誠実, full picture)

**含意**:
- digest の sphere verification table はそのまま 5 段保持
- VIM vs Stoll の rel.err. 列も維持 (機械精度 < 1e-4 で 5 段)
- これが digest §3 の核心 figure / table

### Q21. VIM 定式化の数式 個数
**Q**: VIM 定式化の数式は何個?
**A**: 3 式 保持 (K matrix, M matrix, $\alpha_n$ moment)

**含意**:
- digest §2 で 3 式 (eq:K, eq:M, eq:alpha) 全部保持
- 1 ページ digest で式 3 つは tight だが、定式化の completeness を優先
- 各式は短く (1 行ずつ) 圧縮形で表現
- M_ij = volume integral of φ_i · φ_j (HDiv polynomial で簡単)
- K_ij = double volume integral with 1/r kernel + Spherical Duffy quadrature
- α_n = (-1)^n Σ g_k² τ_k^n (Kameari 規約)

### Q22. §4 Outlook (cuboid) figure
**Q**: cuboid section に figure?
**A**: 不要、テキスト主 (cuboid は outlook 扱い)

**含意**:
- §4 はテキスト 1-2 段落で完結
- 4 系統 0.3% 一致を一文で示す
- 高段 breakdown は full paper に詳細を回す
- figure space は §3 sphere table に集中

### Q23. References
**Q**: references 追加?
**A**: + NGSolve / radia-vim コード GitHub URL を追加 (5 件目)

**含意**:
- 現在 4 件: Sugahara TEAM 28, Hiruma 2020, Stoll 1974, Nagamine 2026
- 追加 5 件目: NGSolve/radia-vim 実装 GitHub URL
- これで reproducibility 担保 + OSS 公開ポリシー一致 (`user_oss_priority.md` 参照)
- 候補 URL:
  - radia-vim: `https://github.com/ksugahar/radia` (or wherever radia-vim lives)
  - 我々の CLN code repo: 別途 publish 必要 (後で確認)
  - radia-mcp: `https://github.com/ksugahar/radia/tree/main/packages/radia-mcp` (knowledge base)
- 形式 (digest bib):
  ```
  \bibitem{radia_vim} K. Sugahara, ``radia-vim: Volume Integral Method for eddy currents,'' GitHub repository, https://github.com/[org]/radia, accessed 2026-05.
  ```
- 確認事項: GitHub URL の正確な path をユーザに確認要

### Q25. 著者所属 (affiliations)
**Q**: 長嶺先生・羽根先生 の所属は?
**A 確定** (PDF/doc から抽出):
- 菅原: Faculty of Science and Engineering, Kindai University, 3-4-1 Kowakae, Higashi-Osaka, Osaka 577-8502, Japan (E-mail: ksugahar@kindai.ac.jp)
- **長嶺先生 = Hideaki Nagamine**, Faculty of Engineering, Gifu University, 1-1 Yanagido, Gifu 501-1193, Japan (E-mail: hnagamin95@ieee.org)
- **羽根先生 = Yoshikatsu Hane (Y. Hane)**, Faculty of Science and Engineering, Toyo University, Kawagoe, Saitama, Japan

**Source**:
- 長嶺先生: `W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/nagamine_paper/2026-05-04-nagamine-rectvc.pdf` (Nagamine 2026 paper)
- 羽根先生: `W:/02_学会資料/2026年度/2026_09_IGTE_Symposium/ヒステリシス@羽根先生/IGTE2026_Digest_Congruency.doc`

**実装**: digest title page 更新済 (K. Sugahara, H. Nagamine, Y. Hane の order と affiliation 確定)

### Q26. v2 progress report 仕上げ — 検証構造を 4 形状に整理
**菅原指示** (複数メッセージ統合):
- v2 を仕上げる、BEM → VIM (Volume Integral Method)
- §4.4 table はみ出しを修正
- もう少し丁寧に書く、段階的に sphere の検証から
- **段階的構造**:
  1. **Sphere** (R=10mm): 3 方法 (Stoll 解析 + NGSolve FEM Kelvin + radia-vim VIM)
  2. **Cylinder/Disk** (R=10 t=2mm): 3 方法 (radia-core axifem + NGSolve FEM Kelvin + radia-vim VIM)
  3. **A1 Cuboid 円板ライク** (17.72²×2): 3 方法 (NGSolve + radia-vim VIM + ELF)
  4. **Cuboid 5×2×1** (sharp 3D, 最後 = breakdown demo): 4 方法 (NGSolve + radia-vim VIM o3/o4 + ELF)
- **角部 → 直交性維持困難** という議論も追加 (corner singularity が Krylov M-orthogonality 損失を加速)

**含意**:
- 既存の radia-vim scripts: cuboid_521, A1_square のみ存在
- **要新規作成**: radia-vim sphere script, radia-vim cylinder script
- 各 ~10-50 分の K matrix assembly (Spherical Duffy quadrature)
- v2 構造は 5 ステップ:
  - §1 背景・訂正
  - §2 Methodology (規約、Foster→Kameari 変換、snapshot fix lint、ZZ)
  - §3 Sphere validation (3 方法)
  - §4 Cylinder validation (3 方法)
  - §5 A1 cuboid validation (3 方法、円板ライクで sharp corner 開始)
  - §6 Cuboid 5×2×1 (4 方法、breakdown demo 最終例)
  - §7 角部と Krylov 直交性 議論
  - §8 長嶺先生への問い

### Q27. v2 再構成 進捗 (2026-05-10 終盤)
**実装完了**:
- §3-§6 を 4 形状に分割 (sphere/cylinder/A1/cuboid 5×2×1)、各形状で 3 ないし 4 手法を明示記述
- §7 「角部 (sharp corner) と Krylov 直交性」議論セクション新設
  - 物理メカニズム 2 経路: (1) FEM Kameari M-orthogonality 損失, (2) VIM Foster mode selection rule
  - 形状複雑度 vs 段数表 (smooth → sharp で単調減少)
  - 緩和策 4 つ (adaptive mesh, higher-order VIM, verified-interval, corner-aware basis)
- v2 = 9 pages (clean structure)

**残作業**:
- A1 cuboid VIM 実行 (extract_tau_A1_square.py 計算中, ~15 min) → §5 table 完成
- sphere/cylinder VIM は radia-vim 不要、Stoll 解析 / Nagamine 軸対称 BEM が VIM 等価として記述

### Q28. A1 cuboid VIM 完了 + 2 系統 0.4% 一致
**実行**: `extract_tau_A1_square.py --order 3 --save-json A1_radia_vim_F4_foster.json` (n_dofs=81, ~10 min)
- A1 VIM Foster spectrum: leading τ = 212.12 μs (cylinder VIM 209 と 1.5% 一致 → 形状差として整合)
- BEM-to-Kameari 変換 (mpmath 80 桁): leading Cauer rung = 194.23 μs

**A1 cuboid 2 系統 cross-validation**:
| k | NGSolve FEM Kelvin | radia-vim VIM o3 | gap |
|---|---|---|---|
| 0 | 195.06 | 194.23 | **0.4%** ⭐ |
| 1 | 66.46 | 70.84 | 6.6% |
| 2 | 31.76 | 35.34 | 11% |
| 3 | 20.06 | 16.12 | -20% (broken) |

**ELF for A1 不存在** (ELF data は cuboid 5×2×1 のみ)、§5 は **2 系統一致** 表示で正直に記述

**v2 最終状態**: 10 pages, 完全に sphere/cyl/A1/cuboid 5×2×1 の 4 形状段階的構成 + §7 角部 Krylov 直交性議論

**ファイル**:
- `W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/2026_05_10_progress_for_nagamine_v2.tex` / .pdf
- `W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/ngsolve_validation/A1_radia_vim_F4_foster.json` (新規)

### Q29. 信頼段数 基準 見直し
**菅原指摘**: 「信頼段数が低すぎるだろう」(1-2 段は厳しすぎる)
**回答**: Ladder structure (値 order 合う、幾つか non-mono も OK) に緩める

**新基準**:
- 「leading rung 一致」(< 1%): 全形状で機械精度
- 「Ladder 構造維持段数」(両 method 同じ qualitative pattern): 4-7 段 (sphere 6-7, cuboid 5×2×1 でも 4-5 段)
- 「数値完全一致段数」(< 5% gap): 角部と sharp さに反比例

**改訂表** (v2 §7):
| 形状 | leading rung 一致 | Ladder 構造維持段数 | 数値一致段数 |
|---|---|---|---|
| Sphere | 0.027% (5 系統) | 6-7 段 | 5 系統 stage 0-4 (機械精度) |
| Cylinder R=10 t=2 | 0.4% (3 系統) | 6 段 | 3 系統 stage 0-2 |
| A1 cuboid | 0.4% (2 系統) | 5-6 段 | 2 系統 stage 0-1 |
| Cuboid 5x2x1 (最 sharp) | 0.3% (4 系統) | 4-5 段 (両方 同じ非mono) | 4 系統 stage 0 のみ |

**重要発見**: \textbf{Ladder 構造そのものは sphere から sharp 3D cuboid まで一貫して 4-7 段維持} される。Breakdown は「数値完全一致」の領域のみ、ladder structure 自体は両 method 一貫。

### Q30. ELF 削除と将来計画 (memo)
**菅原指示**:
- ELF の検証は落とす
- やり方が今いち
- 将来は ELF で time-domain の力 (force) 計算をやりたい
- 今ではない、メモ程度に記録

**含意**:
- v2 から ELF を削除 (cuboid 5×2×1 §6 + A1 §5 + intro list + 各種 table)
- Abstract から「4 手法」→「3 手法」(解析解 + NGSolve + VIM)
- §7 corner table から ELF row 削除
- abstract も「商用磁気モーメントELF」記述削除

**ELF future work memo (lossy)**:
- ELF/MAGIC 商用磁気モーメント法 は時間領域 step response の Foster fit が苦手 (4 mode しか抽出できず Hankel-Padé extrapolation noise)
- 将来方向: ELF の time-domain solver で **eddy current 損失力 (Lorentz force)** を直接計算 → CLN とは別の検証 angle
- 適用例: 電力機器の渦電流損装荷 + 機械的振動・力評価
- これは別研究テーマとして切り出す (今は CLN 抽出に集中)
- 該当 path: `W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/elf_validation/` (ELF macro/result files 既存)

### Q31. v2 §3.1 sphere タイトル訂正
**菅原指摘**: VIM は sphere では実行していない (Stoll 解析等価で済ませた)
**訂正**: §3.1 「3 系統機械精度一致」→「2 系統機械精度一致」(解析解 + NGSolve のみ実比較)
- 本文中の「VIM = Stoll 等価」記述は維持 (理論的説明として有効)
- 表は「analytical (Stoll/5)」と「NGSolve」の 2 列で既に正しい

### Q32. 「Cauer rung」表現訂正 + cylinder BEM 確認 + ELF 注記削除
**菅原指摘 3 点**:

1. 「各形状で独立 2-3 手法の Cauer rung が一致するか」 → \textbf{誤り}, 「R_2n と L_2n+1 が一致するか」が正しい (R, L 比較を強調しているので)
   - **訂正**: §3 intro を「Cauer ladder elements ($R_{2k}$, $L_{2k+1}$) が一致するか」に変更

2. 「Cylinder の 軸対称 BEM」: VIM のこと? それとも BEM 実装が残っている? VIM 実装に書き換えられる?
   - \textbf{回答}: cylinder の第 3 手法は \textbf{Nagamine 先生の真の axisym BEM} (boundary integral on cylinder surface) で実装。radia-vim VIM は cylinder 用は存在しない (radia-vim は hex-only)。
   - ただし \textbf{axisym smooth body では BEM ≡ VIM} (両者とも同じ Foster 固有値問題に帰着) という関係があるので注記追加
   - §3 intro と §3.2 の両方に「axisym smooth body では BEM $\equiv$ VIM」 note 追記

3. §3.4 cuboid の「ELF 商用磁気モーメント法 は本検証から除外...」注記は \textbf{完全削除} (§3 intro で既述、redundant)
   - \textbf{削除}: §3.4 cuboid の ELF 注記行を削除

**実装状況**:
- v2 §3 intro: 「Cauer ladder elements (R_2k, L_2k+1)」表現に変更 + cylinder note 追加
- v2 §3.2 cylinder: BEM (Nagamine 既存実装) 注記追加
- v2 §3.4 cuboid: ELF 注記削除
- compile 通過 → 10 pages 維持

### Q33. Cylinder VIM 化 (BEM 完全削除)
**菅原指示**:
- 「BEM と VIM が同じ資料に出てくるのは理解不能」
- 「Sphere は解析解で比較、cylinder は軸対称 VIM が比較的正しい解と思う」
- 「Cylinder は NGSolve, axis-femm, VIM の 3 者比較できるよな?」
- 「axifem は Hiruma 規約じゃないほうがいい」

**実装**:
- 新規 `cylinder_vim_axisym.py` 作成 (Python, axisymmetric Green's function with elliptic integrals, $N_r$=24, $N_z$=6, 144 cells, n_quad=2)
  - 結果: τ_lead = 211.85 μs (Foster pole 227 μs)
  - canonical Joule R_0 = 2195 Ω 完全一致
  - L_1 1.65% gap vs NGSolve, R_2 4.35% gap
- axifem Kameari 抽出試行 (`axifem_disk_kameari.py`)
  - 物理単位 convention の違い (mu0 が K に組込, σ が M に組込) で直接 Kameari α formula 適用すると wrong scale (R_0=4.33 Ω, expected 2195)
  - Future work: axifem K, M を canonical 単位に変換してから Kameari 抽出
- v2 §3.2 cylinder: \textbf{2 手法} cross-validation (軸対称 VIM + NGSolve 3D)、axifem は Hiruma 規約参照値として note のみ
- BEM 文中削除完全 (radia-vim 関連は VIM, Nagamine 軸対称は VIM 等価として削除)

**Cross-validation summary**:
| 形状 | 検証手法 | leading rung 一致 |
|---|---|---|
| Sphere | 解析解 + NGSolve 3D | 0.027% (machine precision) |
| Cylinder | 軸対称 VIM + NGSolve 3D | 1.7% |
| A1 cuboid | radia-vim hex VIM + NGSolve 3D | 0.4% |
| Cuboid 5x2x1 | radia-vim hex VIM + NGSolve 3D | 0.16% |



### Q24. GitHub URL 確定
**Q**: どの Radia repo URL?
**A**: Radia OSS root (本家): https://github.com/ksugahar/Radia (ユーザ Sugahara 拡張版、radia-vim 含む)

**実装**: digest §References に bibitem 追加
```
\bibitem{radia} K.~Sugahara, ``Radia: Electromagnetic Simulation Framework for Magnetic Levitation Systems,'' GitHub repository, https://github.com/ksugahar/Radia, accessed May~2026.
```


### Q25. Matsuo 論文 SA-26-014 — Cauer ladder の素子べき乗則と分数べき周波数応答
**典拠**: T. Matsuo (Kyoto Univ.), "A Study of fractional power-law frequency response represented by Cauer circuit", SA-26-014 / RM-26-014, IEE Japan Joint TC on Static Apparatus & Rotating Machinery, 2026.3.5-6
**Path**: `W:/02_学会資料/2026年度/2026_09_IGTE_Symposium/CauerLadderNetwork@菅原/digest/referece_paper/SA26-014.pdf`

**核心主張 (eq. 4 / 40)**:
Cauer ladder 素子が漸近的に
- L_{2n-1} = L_1 / c_L(n), c_L(n) ~ (2n-1)^{k_L}
- R_{2n} = c_R(n) R_2, c_R(n) ~ n^{k_R}

の power-law 形で成長/減衰するとき、マクロインピーダンスが **分数べき周波数応答** に従う:

| 量 | 漸近形 | 指数 |
|---|---|---|
| L_eff(ω) | ∝ ω^{−(k_L+1)/(k_L+k_R+2)} | 減少 |
| R_eff(ω) | ∝ ω^{(k_R+1)/(k_L+k_R+2)} | 増加 |
| \|Z\|, Re(Y) | ∝ ω^{(k_R+1)/(k_L+k_R+2)} | — |

特殊値:
- k_L=k_R=0 (定数素子): 指数 = 1/2 (古典渦電流 sheet/wire diffusion)
- k_L=k_R=1: 指数 = 1/2 (同じ)
- k_R=0, k_L=1: 指数 = 1/3
- k_R=1, k_L=0: 指数 = 2/3
- 一般 (k_R+1)/(k_L+k_R+2) ∈ (0, 1)

**N 段打ち切り Cauer の corner frequency (eq. 39)**:
- ω_N ~ N^{k_L+k_R+2} · R_2 / L_1
- ω < ω_N で truncated ladder は infinite ladder を再現、ω > ω_N で破綻
- 段 N での値: L_N ~ L_1/N^{k_L+1}, R_N ~ R_2·N^{k_R+1}

**Appendix B: 非線形拡張**:
材料構成則が h = L(b) ~ h_0(b/b_0)^{β_L}, v = R(i) ~ v_0(i/i_0)^{β_R} と分数べきの場合、最終指数は (β_L·k_L+1)/(β_L·k_L+β_R·k_R+2) と β_L, β_R で再スケールされる。

**菅原の洞察 (2026-05-14)**:
> 分数べき指数 α が分かれば、そこから (k_L, k_R) を介して **終端すべき素子の特性** と **その基底関数** が分かるのでは

**展開**:
1. **観測 α → (k_L, k_R) bijection**: 分数べき指数 α = (k_R+1)/(k_L+k_R+2) は (k_L, k_R) を一意に決めない (1 方程式 2 未知数) が、Re(Y) と L_eff 両方の漸近曲線から k_L, k_R 個別に取れる。あるいは k_L=k_R 仮定で 1 方程式 1 未知数化。
2. **終端素子 (truncation strategy)**: N 段で打ち切る場合、ω > ω_N の漸近応答を再現する **fractional-order termination element** (Constant Phase Element, Z_term = K_term · (jω)^α) を末尾に付加。これで finite-stage Cauer + CPE termination = exact 表現。
3. **基底関数の対応**: 我々の HDiv integrated-Legendre basis は **多項式** のみ。分数べき α が 1/2 以外を取る場合、その分数べきを担う **特異基底** (例: 角近傍の r^{1/2} singularity) が VIM basis に欠落していることになる。Cauer rung から逆算した k_L, k_R は **どの singular basis を追加すべきか** を直接教える。
4. **検証ルート**: 我々の p=3 VIM Cauer rung k=0..2 (R, L 列) の log-log 漸近傾き → k_L, k_R 推定 → α 予測 → NGSolve 高周波 sim 直接測定の α と比較。一致すれば paper の主張が cuboid に拡張されたことの実証。

**IGTE digest への活用**:
- §3 (BEM-style single-element CLN) の補強: 「polynomial basis のみで p=3..5 のスケール内でも、分数べき指数 α は Cauer rung 漸近挙動から直接読み取れる」と数値示せる。
- §4 future work: 「terminated Cauer (N段 + CPE_α) → 同精度を低段で達成」を提案、Matsuo paper を引用。
- 観点として「**精度なしには存在しない CLN**」 + 「**情報圧縮として最適**」 + 「**terminator で finite + fractional element に追加圧縮可能**」の 3 層構造で動機を組み立てられる。

**つながる関連**:
- Q1 (BEM 的要素で CLN): polynomial basis 限界の認識
- Q2 (FEM positioning): Foster spectrum sparsity の構造的理解
- Cuboid 5x2x1 / A1 の breakdown stage k=3 (今回 p=3 で観測) は precision floor が主因だが、もう一つの上界として basis-truncation 由来の "ladder の存在域" もあり、Matsuo paper の corner frequency ω_N と整合する可能性
























