# T-Ω + Kelvin 検証設計書 (cuboid 5×2×1 in vacuum, uniform B_z applied)

**目的**: A-formulation が構造的に向かない uniform B_z 問題に対し、T-Ω
formulation で reduced-potential + Kelvin を実装、ELF τ_lead = 11.51 μs
との照合を目指す。

## 1. 動機

A-formulation reduced + Kelvin の問題点 (v11-v14 で実証):
- v11/v12/v13: (ν - ν₀) reduced-A 形が Kelvin pullback と非互換
- v14 ((ν - ν₀) 修正): A_phys = (B_0/2)(-y, x, 0) が無限遠で発散
  → Convention A pullback 不可、Convention B 適用しても Ampere 由来の
    inner cancellation が機能しない (∇×(ν_0 ∇×A_s) = 0 で空回り)
- 結果: τ_0 = 5985 μs (ELF 11.51 μs の 520×)、stage 1 sign flip

T-Ω は構造的にこの問題を回避できる:

| 領域 | 変数 | 形式階数 | A_s/Ω_s 表現 |
|------|------|----------|-------------|
| 導体内 | T | 1-form | 不要 (J = curl T を解く) |
| 内部 air | Ω | 0-form | Ω_s = -H_0 z (free space) |
| Kelvin kext | Ω | 0-form | Ω_s' = scalar pullback (trivial) |

**核心**: Ω は 0-form なので Kelvin pullback は単純な座標合成
Ω'(r') = Ω_phys(T(r'))。Convention A 問題 (1-form pullback の
特異性) も Convention B 問題 ((ν-ν₀) と Ampere 不整合) も**生じない**。

## 2. 強形式

magnetic field = H_total = T - ∇Ω (T defined only in conductor)

Maxwell:
- ∇ × H = J → ∇ × (T - ∇Ω) = J (= 0 in air, = σ E in conductor)
- ∇ · B = 0 → ∇ · (μ H) = 0 (everywhere)

Reduced potential:
- Ω = Ω_s + Ω_r where Ω_s = -H_0 z (uniform B_z applied)
- Ω_r decays at infinity (compatible with Kelvin pullback)
- T (conductor only)

## 3. Weak Form (FE 実装)

### 3.1 FE 空間

```python
fes_T = HCurl(mesh, order=order, definedon=mesh.Materials("conductor"),
              nograds=True)
fes_Omega = Periodic(H1(mesh, order=order, dirichlet_bbnd="GND"))
fes = fes_T * fes_Omega   # mixed space
```

- T は conductor 限定 → Kelvin 領域に存在しない
- Ω は inner air + kelvin で定義 (Periodic で kelvin_int <-> kelvin_ext 同一視)
- GND: Kelvin 中心 (offset 点) で Ω = 0 → 物理的無限遠で Ω_r → 0

### 3.2 双線形形式

```python
T, Ω = fes.TrialFunction()
W, v = fes.TestFunction()

a += (1/sigma) * curl(T) * curl(W) * dx("conductor", bonus_intorder=8)
a += mu_cf * (T - grad(Ω)) * (W - grad(v)) * dx(bonus_intorder=8)
```

ここで:
- mu_cf = μ_0 in inner, μ' = (R/ρ')² μ_0 in kelvin (Kelvin 変調)

### 3.3 線形形式 (reduced-Ω = -H_0 z + Ω_r)

Ω = Ω_s + Ω_r とおき、Ω_s = -H_0 z (一様外部 H_z = H_0 を生成)。

注意: Kelvin 領域の Ω_s' は **0-form pullback** で:
```
Ω_s'(r') = Ω_s(r_phys) = -H_0 × z_phys
         = -H_0 × (R²/ρ'²) × (z' - o_z)
```

これは ρ' → 0 で発散する! 1-form pullback と同じ問題が再発する...

**修正案**: 代わりに **reduced-potential 用 0-form Convention B**:
```
Ω_s'_kelvin = -(ρ'/R)² × Ω_s_phys(r' - offset)
            = -(ρ'/R)² × (-H_0 × (z - o_z))
            = (ρ'/R)² × H_0 × (z - o_z)
```

これは offset で消滅 (有界)、boundary で sign flip 整合。

実装: `make_reduced_potential_background_cf` の **scalar 版**を別途用意
する必要がある (現状は vector のみ)。

```python
Omega_s_inner = -H_0 * z   # inner: free space H_z = H_0
# kext: Convention B scalar (TODO: add scalar variant of helper)
# 暫定: 手書き
rho_p_sq = (x - ox)**2 + (y - oy)**2 + (z - oz)**2 + 1e-24
factor = -rho_p_sq / R_K**2
Omega_s_kelvin = factor * (-H_0) * (z - oz)
Omega_s_cf = mesh.MaterialCF({
    "conductor": Omega_s_inner,
    "air": Omega_s_inner,
    "kelvin": Omega_s_kelvin,
})

# linear form: drive A via -∇Ω_s contribution
# (J_imp in conductor = σ × (T_s - ∇Ω_s) where T_s = 0 initially)
f += -mu_cf * grad(Omega_s_cf) * (W - grad(v)) * dx(bonus_intorder=8)
```

## 4. Kameari Iteration (T-Ω 版)

```python
Omega_s = -H_0 * z   # known background
gfT_pot = GridFunction(fes_T)   # accumulator for ladder
gfOmega_pot = GridFunction(fes_Omega)
J_imp = sigma * (T_s_initial - grad(Omega_s))  # 初期 impressed current

for n in range(N_STAGES):
    # Solve mixed (T, Ω) for stage n
    # Linear form: (J_imp, W)_cond + reduced-Ω contributions
    solve(...)
    # Compute R_n, L_n, update orthogonalization
    R_n = 1 / Integrate(J_imp * J_imp / sigma, mesh, def="conductor")
    Apot += R_n * gfA_n   # アクセス用 A_total = T - ∇Ω
    L_n = R_n * Integrate(J_imp * Apot * dx, mesh, def="conductor")
    J_imp = J_imp - sigma * Apot / L_n
```

## 5. 期待される改善点

| 問題 | A-formulation v14 | T-Ω 期待 |
|------|------------------|---------|
| (ν - ν₀) 形破綻 | 修正済 | 該当しない (T conductor only) |
| Kelvin pullback の A_s 特異性 | Convention A 失敗、B も部分的 | Ω は 0-form、特異性なし |
| Ampere inner cancellation 不整合 | uniform B_z で空回り | T-Ω は構造的に対応 |
| τ_lead 抽出精度 | 520× off | ELF 11.51 μs に近づく見込み |

## 6. 実装上の課題

1. **Mixed space FE in NGSolve**: `fes_T * fes_Omega` の coupling 形式
2. **make_reduced_potential_background_cf の 0-form scalar 版**: 新たに
   `make_reduced_potential_scalar_cf` を mcp-server に追加する必要あり
3. **Periodic BC for Ω**: GND は Kelvin 中心 (offset 点) に設定
4. **T 空間の境界条件**: T × n = 0 on conductor surface (current confinement)
5. **Foster pole 同定**: ladder {R_n, L_n} → Cauer-II 連分数 → τ_lead

## 7. 検証ステップ

### Step T-Ω-0 (Kelvin なし、有限 air-box)
- T-Ω formulation を air-box で実装 (truncation BC)
- Kameari ladder 抽出 → Tanimoto 修論の T-Ω 結果と比較
- 期待: 既存の T-Ω 実装 (CLN_T_Omega.py) で動作確認済

### Step T-Ω-1 (Kelvin 追加、まだ Convention B なし)
- Kelvin 領域追加 (offset 球)、ν', μ' 変調
- Ω は scalar pullback (0-form 自然形)
- 期待: Step 0 と同じ τ_lead を再現 (有限 air-box → 無限領域への自然な拡張)

### Step T-Ω-2 (full reduced-Ω + Kelvin + Convention B)
- Ω = Ω_s + Ω_r 分離
- Convention B for Ω_s in kext
- Kameari iteration 完成
- 期待: ELF τ_lead = 11.51 μs に <5% 一致

## 8. 工数見積

| Step | 内容 | 工数 |
|------|------|------|
| Step 0 | T-Ω no-Kelvin 実装、CLN_T_Omega.py から adapt | 半日 |
| Step 1 | Kelvin 追加、scalar Ω pullback 検証 | 半日 |
| Step 2 | Full Kameari + Convention B 実装、4-stage 動作 | 1-2 日 |
| 検証 | ELF cross-check + paper 反映 | 半日 |

合計 **2-3 日** 程度の実装作業。本日のセッションでは設計書まで。

## 9. 関連ドキュメント

- `docs/kelvin/KELVIN_TRANSFORMATION.md` §7.5 — (ν-ν₀) 形破綻と修正
- `docs/cln/CAUER_LADDER_NETWORK.md` §6 — CLN + Kelvin 接続、T-Ω 推奨
- `packages/radia-mcp/.../cln_notebooks/CLN_T_Omega.py` — Tanimoto T-Ω 実装
  (Kelvin なし、ベース)
- `docs/kelvin/kelvin_classic_demos.ipynb` /
  `docs/kelvin/kelvin_classic_demos_results.json` —
  archived H-formulation Kelvin 例 (Ω-style, scalar pullback の参考)

## 10. 結論

A-formulation v11-v14 で確認した「Kelvin + reduced potential の落とし穴」は
構造的に T-Ω では回避できる。実装は 2-3 日程度の作業量で、
ELF reference との照合により論文 (cuboid CLN paper draft) の
"open problem" を一つクローズできる見込み。

本日のセッションでは設計書まで。次セッションで Step 0-2 を順次実装。
