# PEEC + CLN + Dowell連分数展開によるコイルモデル縮約

## 概要

本Exampleでは、PEEC（Partial Element Equivalent Circuit）法とCauer Ladder Network（CLN）を組み合わせ、Dowellの連分数展開を用いたコイルモデル縮約の実装と検証を行う。

**理論的背景**: [docs/PEEC_SURFACE_IMPEDANCE.md](../../docs/PEEC_SURFACE_IMPEDANCE.md) を参照。

## 主要な検証結果

### 1. CLN(DC) + Dowell補正 = Dowell式（完全一致）

`cln_with_dowell_correction.py` による検証結果:

| 周波数 | |Z_Dowell| | |Z_CLN+Dowell| | 誤差 |
|--------|------------|---------------|------|
| 1 Hz | 1.724e-04 | 1.724e-04 | 0% |
| 1 kHz | 1.724e-04 | 1.724e-04 | 0% |
| 100 kHz | 2.044e-04 | 2.044e-04 | 0% |
| 10 MHz | 1.229e-03 | 1.229e-03 | 0% |
| 100 MHz | 3.867e-03 | 3.867e-03 | 0% |

**結論**: 最大誤差 < 1e-10%（数値精度限界）

### 2. z*coth(z)の連分数展開

`derive_dowell_cf_algorithm.py` でViskovatovアルゴリズムによる導出を検証:

```
z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + w/(9 + ...))))

連分数係数: b_n = 2n + 1  (n = 1, 2, 3, ...)
           -> 3, 5, 7, 9, 11, 13, ...
```

**数値検証**:

| Delta | z*coth(z) (exact) | J-fraction (N=10) | 誤差 |
|-------|-------------------|-------------------|------|
| 0.1 | 1.00333 | 1.00333 | 3.9e-15% |
| 1.0 | 1.31304 | 1.31304 | 4.8e-14% |
| 2.0 | 2.16396 | 2.16396 | 1.9e-12% |
| 5.0 | 5.03326 | 5.03326 | 2.7e-08% |

### 3. ラダーネットワークへの変換

連分数 `1 + w/(3 + w/(5 + ...))` は以下のRCラダーに対応:

```
        R1        R2        R3
  o----[  ]------[  ]------[  ]------...
         |         |         |
        L1        L2        L3
         |         |         |
        ===       ===       ===
```

ここで `w = tau * s`、`tau = d^2 * mu * sigma / 2`。

## ファイル構成

### コア検証スクリプト

| ファイル | 説明 |
|----------|------|
| `cln_with_dowell_correction.py` | CLN(DC) + Dowell補正とDowell式の一致検証 |
| `derive_dowell_cf_algorithm.py` | Viskovatovアルゴリズムによる連分数係数導出 |
| `derive_dowell_pade.py` | Pade近似による有理関数展開 |
| `derive_dowell_cf.py` | Taylor展開の解析 |

### PEEC-CLN統合デモ

| ファイル | 説明 |
|----------|------|
| `demo_peec_cln_reduction.py` | PEEC + CLN + Dowellによるモデル次数削減デモ |
| `verify_peec_cln_dowell.py` | PEEC-CLN-Dowell統合の検証 |

### 補助スクリプト

| ファイル | 説明 |
|----------|------|
| `analyze_dowell_cln.py` | DowellファクターとCLNの関係分析 |
| `verify_cln_dowell_correct.py` | CLN vs Dowell詳細比較 |
| `verify_cln_dowell_form.py` | CLNからF_R, F_L抽出 |
| `verify_dowell_diffusion.py` | Dowell vs 1D拡散方程式 |

### PEECコイル解析

| ファイル | 説明 |
|----------|------|
| `coil_impedance_peec.py` | PEECによるコイルインピーダンス解析 |
| `coil_on_magnetic_core_peec.py` | 磁性体コア上のコイル解析 |

## 使用方法

### 基本的な検証の実行

```bash
# CLN + Dowell補正の検証
cd examples/peec_integration
python cln_with_dowell_correction.py

# 連分数係数の導出と検証
python derive_dowell_cf_algorithm.py

# PEEC + CLN統合デモ
python demo_peec_cln_reduction.py
```

### 出力ファイル

各スクリプトは `.png` 画像ファイルを出力:

- `cln_with_dowell_correction.png` - CLN vs Dowell比較
- `derive_dowell_cf_algorithm.png` - 連分数展開の検証
- `demo_peec_cln_reduction.png` - モデル次数削減結果

## 理論的背景

### Dowell式

導体の表皮効果を考慮したインピーダンス:

```
Z(s) = R_dc * F_R(xi) + s * L_int_dc * F_L(xi)
```

ここで:
- `R_dc = 1/(sigma*d)` : DC抵抗
- `L_int_dc = mu*d/3` : DC内部インダクタンス
- `xi = d/delta` : 正規化導体厚さ
- `delta = sqrt(2/(omega*mu*sigma))` : 表皮深さ

### F_R, F_L の定義

```
z = (1+j) * Delta    (Delta = xi/sqrt(2))

F_R = Re[z * coth(z)]
F_L = (3/2) * Im[z * coth(z)] / Delta^2
```

DC極限（xi -> 0）で: F_R(0) = 1, F_L(0) = 1

### 連分数展開（J-fraction）

```
z*coth(z) = 1 + w/(3 + w/(5 + w/(7 + w/(9 + ...))))

w = z^2 = 2j*Delta^2 = tau*s
tau = d^2*mu*sigma/2
```

この展開は**厳密**（近似誤差なし）。

## 適用条件

本手法が厳密となる条件:

| 項目 | 条件 |
|------|------|
| 断面形状 | 長方形（1D表皮効果） |
| 材料分布 | 一様（sigma, mu が定数） |
| 連分数次数 | 5-10項で高精度 |

## 磁性体を含む場合の拡張

磁性体（軟磁性コア）がコイル近傍にある場合、インピーダンスは以下のように分解される：

```
Z_total(s) = Z_cond(s) + Z_mag(s)

Z_cond(s) = R_dc * F_R + s * L_int_dc * F_L + s * L_ext_air  [導体: Dowell適用]
Z_mag(s)  = R_mag(omega) + s * L_mag                         [磁性体: Dowell非適用]
```

### Dowell補正の適用範囲

| 成分 | Dowell補正 | 理由 |
|------|------------|------|
| R_dc, L_int_dc | **適用** | 導体内部の表皮効果 |
| L_ext_air | 適用しない | 空気中の磁束（導体外） |
| L_mag | 適用しない | 磁性体内の磁束（導体外） |
| R_mag | 適用しない | 磁性体損失（別メカニズム） |

### 複素透磁率による損失

```
mu = mu_0 * (mu'_r - j * mu"_r)

L_mag = L_mag_0 * mu'_r    [リアクタンス]
R_mag = omega * L_mag_0 * mu"_r / mu_0  [磁性体損失]
```

### 関連スクリプト

| ファイル | 説明 |
|----------|------|
| `coil_on_magnetic_core_peec.py` | 磁性体コア上のコイル解析（CplMagソルバ） |
| `esim_conductor_model.py` | ESIM表面インピーダンスモデル |

詳細な理論は [docs/PEEC_SURFACE_IMPEDANCE.md](../../docs/PEEC_SURFACE_IMPEDANCE.md) を参照。

## 利点

### 周波数領域
- DCから高周波まで厳密な周波数特性
- 各周波数でF_R, F_Lを直接計算可能

### 時間領域
- FFT/IFFTが不要
- 非線形材料との結合が可能
- SPICEとの連携が容易

### モデル次数削減
- Lanczosで効率的に削減
- リアルタイムシミュレーション向け
- 制御系設計への応用

## 参考文献

1. P.L. Dowell, "Effects of eddy currents in transformer windings," Proc. IEE, 1966.
2. J.A. Ferreira, "Electromagnetic Modelling of Power Electronic Converters," 1989.
3. A. Ruehli, "Equivalent Circuit Models for Three-Dimensional Multiconductor Systems," IEEE Trans. MTT, 1974.
