# A-V formulation Cauer ladder extraction — 実装計画ドラフト

**作成**: 2026-05-09
**目的**: 3D rectangular conductor (cuboid 5×2×1mm Cu) の vacuum-coupled Cauer ladder を、Hiruma 3-term の cuboid gauge collapse を回避して抽出する algorithm の設計

## 背景

- **Tanimoto A-T**: PEC 内部近似 (vacuum coupling 不可) → drop
- **Hiruma 3-term + Kelvin (v5)**: cylinder で 0.02% 動くが cuboid で gauge collapse (-76%) → 3D 研究課題
- **A-V formulation**: 電気スカラー V を導入して div(σE) = 0 を陽に強制、gauge robust → 候補

## A-V formulation 数式

### 弱形式 (frequency domain, $s = j\omega$)

\textbf{基礎方程式}:
- $\nabla\!\times(\nu\nabla\!\times\mathbf{A}) + \sigma s\mathbf{A} + \sigma\nabla V = \mathbf{J}_\text{src}$ on conductor + air
- $\nabla\!\cdot[\sigma(s\mathbf{A} + \nabla V)] = 0$ on conductor (continuity of current)

\textbf{弱形式}:
\begin{align}
  &\int_\Omega \nu \nabla\!\times\mathbf{A}\!\cdot\!\nabla\!\times\mathbf{A}'\,dV
   + s\int_{\Omega_c} \sigma\mathbf{A}\!\cdot\!\mathbf{A}'\,dV
   + \int_{\Omega_c} \sigma\nabla V\!\cdot\!\mathbf{A}'\,dV = 0 \\
  &s\int_{\Omega_c} \sigma\mathbf{A}\!\cdot\!\nabla V'\,dV
   + \int_{\Omega_c} \sigma\nabla V\!\cdot\!\nabla V'\,dV = 0
\end{align}
($\mathbf{A}' \in V_h^A$, $V' \in V_h^V$ は test function)

### 離散空間

```python
fesA = Periodic(HCurl(mesh, order=ORDER, dirichlet_bbnd="GND", nograds=True))
fesV = H1(mesh, order=ORDER, definedon=mesh.Materials("conductor"))
fes = fesA * fesV   # mixed space
(A, V), (A_test, V_test) = fes.TnT()
```

- fesA: HCurl on full mesh (cond + air + Kelvin), tree-cotree masking 不要 (V が gauge を吸収)
- fesV: H1 on conductor only (V は σ ≠ 0 領域でのみ意味あり)

### 行列構造

ブロック行列 $\begin{pmatrix} K + s M_{AA} & D \\ s D^\top & K_V \end{pmatrix} \begin{pmatrix} \mathbf{A} \\ V \end{pmatrix} = \begin{pmatrix} \mathbf{f}_A \\ \mathbf{f}_V \end{pmatrix}$

- $K = \int_\Omega \nu_\text{cf}\,\nabla\!\times\!\cdot\nabla\!\times\!\cdot\,dV$ (Kelvin 含む)
- $M_{AA} = \int_{\Omega_c} \sigma\,\cdot\!\cdot\,dV$
- $D = \int_{\Omega_c} \sigma\nabla(\cdot)\!\cdot\!\cdot\,dV$ (A と V の coupling)
- $K_V = \int_{\Omega_c} \sigma\nabla(\cdot)\!\cdot\!\nabla(\cdot)\,dV$

### Cauer ladder 抽出

\textbf{方針 1: Hiruma 3-term on (K_full, M_full)}
- 拡張行列 $K_\text{ext} = \begin{pmatrix} K & D \\ D^\top & K_V \end{pmatrix}$, $M_\text{ext} = \begin{pmatrix} M_{AA} & 0 \\ 0 & 0 \end{pmatrix}$
- $K_\text{ext}^{-1} M_\text{ext}$ の Lanczos で Foster poles 抽出
- gauge 問題: $K_\text{ext}$ の kernel に何があるか? V が gauge を吸収するかチェック必要

\textbf{方針 2: V を schur complement で消去**
- A 方程式から V を解析的に消去:
  $V = -K_V^{-1} (s D^\top \mathbf{A})$
- A 方程式に代入: $(K + s M_{AA} - s D K_V^{-1} D^\top) \mathbf{A} = \mathbf{f}_A$
- 実効 K + 実効 M を Lanczos で 3-term 反復
- $M_\text{eff} = M_{AA} - D K_V^{-1} D^\top$ (現状 σuv より構造が違う)

## 実装計画 (3 段階)

### Phase A1: 単純テスト (1-2 日)
\textbf{目標}: A-V mixed space の動作確認、cylinder で v5 と一致するか

```python
# tanimoto_AV_kelvin.py (新規)
- mesh: cyl + air + Kelvin (v5 と同じ build_geo)
- fes: fesA * fesV (mixed)
- 弱形式: 上記
- Direct frequency sweep (ω = 1e3, 1e4, ..., 1e8 rad/s)
- Y(iω) = (induced moment) / B_0
- Foster fit (rational interpolation, scipy)
```

\textbf{基準値}:
- Cylinder leading τ = 218.62 μs (v5 reference)
- Foster modes [τ_0, τ_1, ...] = [218.62, 78.21, 39.59, 23.12, 16.07, 13.30] μs

\textbf{基準点}: cylinder で 1\% 以内に再現できれば A-V formulation の正当性確認。

### Phase A2: Cauer ladder 直接抽出 (3-5 日)
\textbf{目標}: A-V で Cauer 反復を直接実装、$K_\text{ext}^{-1} M_\text{ext}$ の Lanczos

```python
# AV_cauer_lanczos.py
- 拡張 K, M を組み立て
- Lanczos 3-term recurrence on K^{-1}M
- Foster poles {τ_n} と Cauer (R_n, L_n) 抽出
- v5 Hiruma 3-term と直接比較
```

\textbf{checkpoint}:
- cylinder で v5 と $\le 0.5$\% 一致 → \textbf{A-V formulation 確証}
- 一致しない場合: V 方程式の BC, schur complement の取扱い等を debug

### Phase A3: Cuboid extension (1 週間)
\textbf{目標}: cuboid 5×2×1mm で gauge collapse なしに Cauer ladder 抽出

\textbf{予測}:
- A-V は V 方程式が gauge を吸収 → tree-cotree なしでも動作
- cuboid corner 形状でも sharp transition なし (V は H1)
- v5 で破綻していた cuboid 11 μs を再現できるか

\textbf{比較対象}: BEM-Foster Phase F-4 (10.85 μs) + ELF Foster (11.04-11.51 μs) の 4-way consensus 11.0 μs。

## 参考コード

\textbf{者多部 (Tanimoto) thesis 関連}:
- W:/00_CAE/NGSolve/谷本/修論/CLN_AT.ipynb (A-T 元コード、wire skin-effect)
- W:/00_CAE/NGSolve/谷本/修論/CLN_APhi.ipynb (A-Φ formulation、A-V に近い)
- W:/00_CAE/NGSolve/谷本/修論/CLN_T-Omega.ipynb (T-Ω formulation、別アプローチ)

\textbf{v5 Hiruma 3D HCurl + Kelvin (mesh setup re-use)}:
- ngsolve_validation/disk_3d_kameari_hiruma_v5_orderphi.py

\textbf{Kelvin helpers}:
- S:/Radia/01_GitHub/src/radia/kelvin_geometry.py
- S:/Radia/01_GitHub/src/radia/kelvin_material.py

\textbf{NGSolve mixed space example}:
- NGSolve i-tutorial 2.7 (mixed Helmholtz)
- NGSolve i-tutorial 6.1 (Maxwell time-harmonic)

## リスクと開放問題

1. \textbf{V 方程式の BC}: ∂V/∂n = ? on conductor surface. 通常は σA·n (electric current normal continuity) だが Kelvin 領域接続でどう扱う?
2. \textbf{Schur complement 数値安定性}: $K_V^{-1}$ が ill-conditioned だと M_eff も悪化。preconditioner 必要かも
3. \textbf{Lanczos break-down}: A-V mixed space で Krylov 直交性が保たれるか未検証
4. \textbf{Foster pole 抽出時の規格化}: R_0 を B=1T or A_s=Wb/m どちらで規格化するか (菅原指示で要決定)
5. \textbf{初段 R_0 の物理的意味}: ∫σ|A_s|² か別の量か、A-V では V の loading effect が入る

## 工数見積もり

- Phase A1 (単純テスト): 1-2 日
- Phase A2 (Cauer 直接抽出): 3-5 日
- Phase A3 (cuboid extension): 1 週間
- 合計: 約 2 週間

\textbf{本セッションでは Phase A1 まで実施可能}。Phase A2-A3 は別セッションで継続。

## メモ条件 (本実装でも適用)

- Kelvin 必須: Phase A1 から air + Kelvin sphere 含めて mesh build
- L = A·J 積分: Phase A2 で Cauer ladder 抽出時
- A_s に Helmholtz-Hodge 補正: Phase A1 から source 構築時に適用
- A_s 由来エネルギーは CLN に含めない: B_acc, J_acc 累積で induced のみ
- BEM と FEM の初段一致確認: Phase A2 完了後に cylinder で確認
- R_0 規格化方針 (B=1T or A_s=Wb/m): Phase A1 完了時点で決定要
