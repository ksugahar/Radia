# Obika君向け HDiv-MMM加速器電磁石最適化 研究引継書

最終更新: 2026-08-21

対象リポジトリ: `S:\Radia\01_GitHub`

現在の作業ブランチ: `backup/main-pre-release-20260821`

## 1. 研究の最終目標

入力として

1. 設計軌道
2. 設計軌道まわりで実現したい任意のtransfer matrix
3. 固定非設計領域（ビーム軌道、ギャップ、コイル、支持・加工制約）

を与え、HDiv-MMMの磁性体要素を使って加速器用偏向電磁石の鉄形状を逆設計する。
当面は円形加速器の閉軌道問題ではなく、固定した設計軌道に沿う1-pass問題を対象とする。

HDiv-MMMを選ぶ研究上の狙いは、空気領域を体積メッシュせず、鉄要素の追加・削除・変形が
磁場と軌道光学へ与える影響を重ね合わせとSchur補完で高速に評価できる点にある。最終的には
加速器電磁石だけでなく、モータなどの磁性体形状・トポロジー設計へ発展させる。

## 2. 現在採用している設計思想

最適化は次の4段カスケードとして整理した。

1. 必要なsection transfer matrix
2. 必要な軌道上磁場・多極場応答
3. 必要な要素磁化・鉄充填率
4. 必要な滑らかな鉄境界形状

ただし、段2で磁場の最小ノルム解を先に一意に選んではいけない。同じtransfer matrixを作る
磁場変化は多様体をなすため、材料として作りやすい磁場を段3と同時に選ぶ。実装上は、場から
section仕様への解析・AD Jacobianを `J`、1要素を完全な鉄として励起したときの場応答を `E`
として `J @ E` を直接阿部ACA--QR--TSVDへ渡す。必要磁場は解いた後に `E @ fill` として
読み出す。このため4段の物理的説明は維持されるが、段2と段3の数値解法は材料可製造性を介して
結合される。

設計変数はHDiv DOFではなく、必ず**1要素あたり1個の符号付き鉄充填率**とする。

- 既存鉄要素: 削除可能容量 `[-1, 0]`
- 充填可能な空気セル: 追加可能容量 `[0, 1]`
- 正の解: 鉄を追加
- 負の解: 鉄を削除

HDiv DOFはベクトル場の成分であり、DOFごとに `[0,1]` や `[-1,0]` を課してはいけない。
この誤りは過去に30720 DOF中15504 DOFのclipと反復発散を生んだ。

## 3. 二つの形状更新経路

### 3.1 Lego・Schur経路

候補HEX要素を追加・削除する粗い探索である。全候補をまとめた応答にACA--QR--TSVDを使い、
Schur補完で候補バッチを評価し、完全solveで採否を決める。BESOは廃止した。固定非設計領域、
成長の連結性、空気領域の外部連結性を制約として持たせる。Cubit SculptはこのLego探索の後処理で
あり、Lego状態の探索そのものではない。

### 3.2 阿部・滑らか境界経路

連続した符号付き充填率を物理体積に換算し、面binの面積で割って磁極面高さへ保存的に変換する。
aperture面とpole rootを固定したblendでGetTrafo変形し、品質限界を超える場合だけCubit再生成へ
進む。完全HDiv solve、軌道追跡、native transfer map評価を行い、元のengineering-band目的が
改善し、bend/orbit/mesh guardを通った形状だけを採用する。

両経路は競合ではない。Lego・Schurは高速で広い離散探索、阿部・GetTrafoは滑らかな境界を得る
連続探索として比較・融合する価値がある。

## 4. 恒久化済みの実装

### Python

- `src/radia/accelerator_abe_topopt.py`
  - `measured_element_fill_patterns`
  - `contract_hdiv_element_fill_response`
  - `compose_specification_fill_response`
  - `solve_abe_element_fill_plan`
  - `bin_element_fill_to_interface_height`
  - `blended_interface_displacement`
  - `optimize_abe_section_contour`
- `src/radia/accelerator_section_optics.py`
  - 固定設計軌道上のsection optics逆問題
- `src/radia/accelerator_magnet_topopt.py`
  - 設計軌道・transfer matrix・HDiv-MMM材料更新の上位構成
- `src/radia/accelerator_taylor_topopt.py`
  - `R/T/U` Taylor mapへの拡張
- `src/radia/stream_function.py`
  - 阿部/DUCAS mode selectionと境界付き反復
- `src/radia/topology_optimization.py`
  - VIM線形化、Schur、LP、解析形状微分
  - `hdiv_mmm_all_single_removal_responses`
    - 1回の候補Schur分解から全単一削除応答を厳密評価するbinary oracle
    - 候補ごとにほぼ同じSchur行列を再分解しない
  - `hdiv_mmm_removal_group_responses`
    - 同じSchur逆ブロック公式を任意の削除グループへ一般化
    - 全候補ペアなどを追加H-matrix solveなしで比較する

### C++・MEX

- `src/core/rad_stream_function.h`
- `src/core/rad_stream_function.cpp`
  - HACApK ACA+
  - thin QR + small TSVD
  - box-bounded repeated Abe correction `SolveAbeBounded`
- `src/matlab/radia_mex.cpp`
  - `topopt.abe_element_fill_plan`
- `matlab/+radia/+topopt/solveAbeElementFillPlan.m`
  - standalone MEXの公開MATLAB入口
- `matlab/+radia/+topopt/contractHDivElementFillResponse.m`
- `matlab/+radia/+topopt/composeSpecificationFillResponse.m`
- `matlab/+radia/+topopt/binElementFillToInterfaceHeight.m`
- `matlab/+radia/+topopt/blendedInterfaceDisplacement.m`
- `matlab/+radia/+python/acceleratorAbeTopopt.m`
  - NGSolve callbackを含む外側反復の明示的なbatch fallback

MEXはPythonを起動せず、pybindと同じHACApK ACA+・QR--TSVD C++核を使う。境界付き反復の
C++関数はMEXから利用済みだが、現時点ではまだpybindへ直接公開していない。Pythonの
`solve_abe_element_fill_plan`は既存Python境界反復を使うため、将来、C++ bounded solverを
pybindへ薄く公開してPython/MATLABの反復も完全な単一数値源にする余地がある。

### 試験・検証記録

- `tests/test_accelerator_abe_topopt.py`
- `tests/test_stream_function.py`
- `tests/test_accelerator_section_optics.py`
- `tests/test_accelerator_magnet_topopt.py`
- `tests/test_topology_optimization.py`
- `tests/matlab/test_radia_mex.m`
- `packages/radia-mcp/tests/test_matlab_radia_mex_contract.py`
- `validation_test/ffag_topopt/README.md`
- `validation_test/ffag_topopt/validation_abe_manufactured_edge_cell.py`
- `まとめ.md`

## 5. 再現済みの数値結果

保存済みTURBO磁石データで次を確認した。

- 観測数: 170
- 設計要素数: 2560
- 追加候補空気セル: 713
- 削除候補鉄セル: 1847
- 仕様直接応答の保持rank: 2
- 保持条件数: 123.1
- 仕様の線形到達可能率: 100%
- 段3の厳密map予測達成率: 98.2%
- gross材料体積: 0.2374 cm3
- net材料体積: -0.0397 cm3
- implied field peak: 1.13e-4

恒久APIと旧scratch実装の差は次である。

- 要素応答縮約: `1.626e-19`
- `J @ E` 仕様合成: `1.821e-17`
- 充填率: `3.660e-16`
- 必要磁場: `9.148e-20`

しかし、充填率を滑らかな2D磁極面高さへ変換し、完全HDiv solveとtransfer map再評価まで行った
段4の達成率は**55.3%**に落ちた。段3の98.2%を「最終形状が98.2%達成した」と報告しては
ならない。

比較値:

| 経路 | 達成率 | 材料 | field fit | 形状 |
|---|---:|---:|---:|---|
| 二値Lego/LP | 90.1% | 0.031 cm3 | 8.77e-3 → 1.54e-2 | 9要素のギザギザ |
| 阿部・滑らか輪郭 | 55.3% | 0.237 cm3 | 8.77e-3 → 8.77e-3 | 滑らかな面 |

阿部経路は現在、場の表現を劣化させず滑らかな面を作るが、3D材料需要から2D面厚みへの変換効率が
56.3%しかない。Lego経路はtransfer matrixをよく合わせるが、場fitと加工形状を悪化させる。

### 5.1 答えを製造した単純物理問題での成功

実問題へ進む前のcalibration benchmarkとして、入口端の1 HEXを除去した完全HDiv-MMM解を
target transfer matrixとし、全鉄状態からその要素を逆同定する検証を追加した。これはfake evaluator
ではなく、正解生成・最終照合とも完全HDiv-MMM solveである。optimizer内の設計有限差分は使わない。

96 HEX、BDM1 3456 DOF、入口端16候補のLAB smokeで次を確認した。

- 直接磁化寄与だけの候補列: 正解セルは8位で失敗
- Schur再分布を含む連続ACA--QR--TSVD係数: 正解セルは5位
- 同じ縮約Schur行列上の全単一削除oracle: 正解セルは1位
- 選択セル: target element 49と一致
- 最終完全solveのtarget transfer matrixに対する最大band比: 0.0
- ACA--QR--TSVD時間: 0.237秒
- 新しい全単一削除oracle時間: 1.073秒
- 全体のLAB smoke観測: 141.2秒から68.9秒へ短縮

最後の時間値はLAB上のsmoke観測であり、公開ベンチ値にしてはならない。mdxとhibinoは実行時に
SSH timeoutで到達できなかったため、compute-host benchmarkは未実施である。重要なアルゴリズム上の
結論は、連続TSVD係数の絶対値最大をそのまま1セルへ丸めてはいけないことである。相関列ではTSVDは
正しい低rank方向を複数セルへ分配する。TSVDは符号・候補部分空間を決め、最終の離散セルは縮約Schur
oracleで選ぶ。

同じ96 HEX問題を既知2セル削除targetへ拡張した結果も合格した。

- 設計候補: 16セル、全削除ペア: 120
- 製造正解: element `[24, 36]`
- 連続ACA--QR--TSVDにおける正解セル順位: 1位・2位
- 縮約Schur group oracleにおける正解ペア順位: 1位
- 最終完全solveの最大band比: 0.0
- 120ペアのgroup oracle: 0.0315秒
- LAB smoke全体: 27.6秒

これらの時間も公開ベンチ値ではない。単一削除と複数削除を同じSchur逆行列から評価できること、
2セルの協調効果を単一セル係数の単純丸めに頼らず復元できることが検証上の成果である。

## 6. これまで分かった失敗原因

### 6.1 設計空間が既存鉄だけだった

応答行列を既存鉄だけから作ると、要求の100%が既存鉄内部へ落ち、外側へ成長できない。
設計空間は `removable iron | addable air` の両方を必ず含める。

### 6.2 HDiv DOFを材料変数にしていた

各HDiv DOFへ材料容量を課すのは物理的に誤りである。完全要素の磁化パターンで観測行を縮約し、
1要素1列の応答を作る。

### 6.3 到達不能な残差までTSVD modeを積んだ

小さい特異値まで使うと係数が物理容量を超え、全候補がclipされ、clipが作った誤差を次の反復が
追いかけて発散する。残差目標、保持rank、条件数、clip数、gross/net材料量を必ず同時に監視する。

### 6.4 transfer matrix誤差から一意な磁場目標を先に作った

磁場応答の条件数は約1e10だったが、仕様を材料へ直接問う `J @ E` は123.1まで改善した。
必要磁場は先に固定せず、材料解から事後に読む。

### 6.5 2粒子だけでetaとeta-primeを評価した

最終の軌道光学評価では設計運動量の前後に複数粒子を置き、設計軌道まわりのmapを評価する。
transfer matrixは設計軌道が決まらないと定義できないため、軌道復元とmap評価の順序を崩さない。

## 7. 厳守する実装境界

1. NGSolveが要素orientation、local/global DOF変換、Piola写像、曲面写像、quadrature、
   `CoefficientFunction`、`GridFunction`、GetTrafoを所有する。PythonでFE plumbingを再実装しない。
2. 非鉄セルへ鉄セルのローカル磁化パターンを移すときは、NGSolve側の`pattern_transfer`を使うか、
   構造HEXでローカルDOF順が互換と検証済みであることを明示する。
3. 最適化感度に有限差分を使わない。有限差分は解析・AD微分の回帰照合だけに使う。
4. 各形状の採否は、完全HDiv solveとexact/native transfer map再評価で決める。
5. 近似予測だけで最終達成を宣言しない。
6. Python FE処理の`ngsolve.TaskManager`はcaller-ownedとする。
7. C++/Python内部配列はrow-major契約を維持し、MEX境界だけでMATLAB column-majorへ変換する。
8. 軌道、ギャップ、コイルは固定非設計領域とし、ヨーク・磁極端部だけを動かす。
9. 穴を絶対禁止にはしないが、加速器磁石では連結性・加工性を優先し、穴は積極的な設計自由度に
   しない。鉄連結性と外部空気連結性を検査する。
10. Scratchは`C:\temp`だけを使う。重いvalidationはidleなmdx、またはhibinoで実行する。
11. 共有worktreeの他ユーザー・他agentの未コミット変更を変更・削除・commitしない。

## 8. Obika君に最初に進めてほしい中核タスク

### Task A: 実規模の外側再線形化反復を閉じる

まず上記manufactured benchmarkを次の順に拡張してから、
`optimize_abe_section_contour`を実際のTURBO/偏向電磁石モデルへ接続する。

1. **完了**: 既知の2セル削除targetを全候補Schur group oracleで復元する。
2. **完了**: 固定設計軌道と1個の滑らかなGetTrafo pole-face modeで、既知振幅から作ったtarget mapを復元する。
   BDM1 TET 72要素・546 DOFで4.000 mmの既知係数を3.99967 mmまで回収し、fresh full solveの
   最大band比0.854を得た。`dM+dB+dG+dC+drhs`はすべて解析微分で、optimizer有限差分は0回である。
3. **次に実施**: 入口・出口の2 modeへ拡張し、軌道復元を入れる。
4. 4～8個の滑らかなmodeで任意mapの小摂動を再現する。

単純問題で壊れた時点でTURBO実問題へ進まず、candidate model、離散oracle、形状実現、軌道復元の
どの段で失敗したかを分離する。

必要なcallback:

- `realize_fill(accumulated_fill)`
  - 参照形状から累積fillを絶対値として実現する。
  - trialごとに前trialの変形へ増分を重ねない。
  - fill→保存的面高さ→固定面blend→GetTrafoを第一選択とする。
  - 品質限界を越えたときだけCubitで再生成する。
- `evaluate_exact(realization)`
  - 完全HDiv-MMM solve
  - 設計軌道の復元・確認
  - native transfer matrixまたはR/T/U map評価
  - bend/orbit/mesh品質を含むpayloadを返す。
- `relinearize(exact, realization)`
  - 受理した形状で `E` と `J` を解析・ADで再構築し、新しい `J @ E` を返す。

最初の目標は、段4の55.3%を1回以上の受理済み再線形化で有意に改善することである。
「有意」は推定ではなく、完全solve後のengineering-band比で判定する。

### Task B: 3D材料需要→2D磁極面の損失を分解する

現在最大の損失は段3の98.2%から段4の55.3%への低下である。少なくとも次を別々に測る。

1. 同一bin内の正負fill相殺
2. 上下磁極または入口・出口を同じ面binへ畳むことによる相殺
3. 面平滑化による高空間周波数の消失
4. GetTrafo変形後の磁化再分布
5. 固定aperture/root blendによる有効変位の減衰
6. 線形化後のtransfer map非線形誤差

改善案は、符号・磁極面・長手方向を分離した保存的basis、trust region、複数の滑らかな面basis、
受理後再線形化である。まず損失を数値分解し、思いつきで平滑化係数を調整しない。

### Task C: Lego・Schurとの公平な比較

同じ固定設計軌道、同じtransfer matrix target、同じ非設計領域、同じ完全solve acceptanceで

- Lego・Schur
- 阿部・滑らか輪郭
- Legoで粗探索後に阿部/GetTrafoで仕上げるhybrid

を比較する。達成率だけでなく、field fit、材料量、形状滑らかさ、連結性、solve時間、H-matrix再構築
回数を記録する。

### Task D: bounded solverの完全なC++単一数値源化

`radia::stream_function::SolveAbeBounded`をpybind11へ薄く公開し、Pythonの
`solve_abe_element_fill_plan`からも利用できるようにする。Python固有callbackやNGSolve objectは
C++へ持ち込まない。小規模dense golden、ACA経路、既存Python実装、MEXの4者を同一入力で比較し、
fill、残差履歴、clip履歴、rank、条件数を照合する。

## 9. 完成判定

PoCを完成と呼ぶには、最低限次を満たす。

1. 設計軌道とtarget transfer matrixを入力として受け取る。
2. transfer matrixはその設計軌道まわりで評価される。
3. `J @ E` のrank、条件数、到達可能射影を記録する。
4. 追加・削除の両候補を含む要素充填解を得る。
5. 充填解を実形状へ変換し、完全HDiv-MMM solveを行う。
6. 完全solve後のtransfer matrixが指定bandを満たす。
7. bend angle、基準軌道、ギャップ、コイル非干渉、mesh品質、鉄・空気連結性が合格する。
8. 少なくとも一つの独立経路で磁場または軌道mapをcross-checkする。
9. 全反復のcommit、host、runtime、rank、残差、材料量、形状品質をresult JSON/logへ保存する。
10. 有限差分をoptimizerに使用していない。

現時点では条件1～4の基盤に加え、製造正解1セル・2セル問題と単一GetTrafo pole-face modeでは
条件5～6を完全solveで満たした。ただし、後者は解析`dC/dq`が完成しているBDM1 TETの
calibration benchmarkであり、HEX configured-field微分、複数mode、軌道復元、および実規模で
条件6を満たす外側反復は未完である。

## 10. 再現コマンド

### Python回帰

```powershell
python -m pytest -q `
  tests/test_accelerator_abe_topopt.py `
  tests/test_stream_function.py `
  tests/test_accelerator_section_optics.py `
  tests/test_accelerator_magnet_topopt.py `
  tests/test_topology_optimization.py
```

2026-08-21時点で151/151件合格。

新しいSchur逆ブロックoracleのfocused回帰:

```powershell
python -m pytest tests/test_topology_optimization.py -k `
  "removal_oracle or block_schur_bundle" -q
```

2026-08-21時点で3/3件合格。単一削除、2セル削除、既存block Schur照合を含む。

### 製造正解2セル物理検証

```powershell
python -u validation_test\ffag_topopt\validation_abe_manufactured_edge_cell.py `
  --output C:\temp\abe_edge_96hex.json `
  --nx 8 --ny 6 --nz 2 --segments 16 --threads 16 `
  --candidate-model schur --target-cell-count 2
```

`status=pass`、`target_elements=[24,36]`、`schur_oracle_target_rank=1`、
`exact_max_band_ratio=0.0`を確認する。重い正式実行はidleなmdx/hibinoで行い、hostとruntimeを
result JSONへ残す。

### 製造正解GetTrafo transfer-map検証

```powershell
python validation_test\ffag_topopt\validation_gettrafo_manufactured_transfer_map.py `
  --output C:\temp\manufactured_gettrafo_transfer_map.json
```

LABの小規模実装smokeでは`status=pass`、72 TET、BDM1 546 DOF、既知4.000 mmに対する
回収値3.99967 mm、fresh full solveの`final_exact_max_band_ratio=0.854`を確認した。これは
速度benchmarkではない。multi-orbit観測はrow-major packed contractを使い、通常のbatched row
builderとの差は最大`6.36e-21`、optimizerの有限差分使用は0回である。

### MEX build

```powershell
pwsh -ExecutionPolicy Bypass -File .\Build.ps1 -MatlabMexOnly -Verbose
```

生成先は `matlab/radia_mex.mexw64`。バイナリを手動コピーせず、Build.ps1のpost-buildを使う。

### MATLAB focused回帰

```powershell
matlab -batch "addpath('matlab'); clear mex; s=testsuite('tests/matlab/test_radia_mex.m'); n=string({s.Name}); r=run(s(contains(n,'testAbeElementFillPlanNative'))); disp(r); assert(all([r.Passed]));"
```

2026-08-21時点で1/1件合格。dense参照とACA経路の両方を含む。

### Python/MCP MEX契約

```powershell
python -m pytest -q `
  tests/test_accelerator_abe_topopt.py `
  tests/test_stream_function.py `
  tests/test_simulink_release_package.py `
  packages/radia-mcp/tests/test_matlab_radia_mex_contract.py
```

2026-08-21時点で63/63件合格。公開MEXコマンド数は363。

## 11. 重要なcommit

- `2ae5cbc7d` `feat(topopt): promote Abe element-fill contour optimization`
- `484a3f91d` `feat(matlab): add native Abe element-fill MEX solve`
- `d2eb092a7` `feat(topopt): accelerate exact Schur removal oracle`
- `5553beb62` `feat(topopt): evaluate collaborative Schur removals`
- `86c14d6dd` `feat(topopt): close manufactured GetTrafo map inverse`

この5本は現在のブランチ上で直列になっている。

## 12. Scratchの扱い

過去の試行は `C:\temp\da_topopt` にある場合がある。

- `turbo_abe_realizable.py`
- `stage3_fast.py`
- `turbo_abe_closed_loop.py`

これらは研究経緯の調査には使えるが、恒久APIのsource of truthではない。再利用価値のある計算は
`src/`、高速回帰は`tests/`、重い実問題は`validation_test/ffag_topopt/`へ昇格する。
`examples/`は廃止済みなので作らない。

## 13. Obika君へ渡す依頼文

以下をそのまま次の作業依頼に使える。

> `S:\Radia\01_GitHub\obika_handover.md`を最初から最後まで読み、commit
> `2ae5cbc7d`、`484a3f91d`、`d2eb092a7`、`5553beb62`、`86c14d6dd`の実装を確認してください。まず既存のPython・MEX focused回帰を
> 再現してください。次に`validation_abe_manufactured_edge_cell.py`をcompute hostで再現し、
> `hdiv_mmm_all_single_removal_responses`が全候補ごとのSchur再分解をせず正解セルを1位にすることを
> 確認してください。続いて`validation_gettrafo_manufactured_transfer_map.py`で単一modeの解析形状微分と
> fresh full solveを再現してください。その後、入口・出口2 modeへ拡張して軌道復元を入れてから、
> 保存済みTURBO/偏向電磁石モデルに`optimize_abe_section_contour`の
> `realize_fill`、`evaluate_exact`、`relinearize` callbackを接続してください。目的は、段3の
> 98.2%予測を最終形状の達成と誤認せず、現在55.3%の段4完全閉ループを、受理済み外側再線形化
> により改善することです。各反復でrank、条件数、clip履歴、gross/net材料量、面高さ、mesh品質、
> bend/orbit guard、完全solve後transfer matrixを記録してください。有限差分はoptimizerに使わず、
> NGSolveのFE plumbingをPythonで再実装しないでください。重い計算はidleなmdxまたはhibinoで行い、
> 既存の未コミット変更を壊さないでください。
