# SIBC 検証完了: 2D 軸対称 Kelvin FEM で Full-Resolution と Robin BC が一致

**日付**: 2026-04-14  
**担当**: 菅原  
**リポジトリ**: [ksugahar/Radia](https://github.com/ksugahar/Radia)  
**検証スクリプト**: [`examples/eddy_current_analytical_validation/reference_2d_axisym.py`](https://github.com/ksugahar/Radia/blob/main/examples/eddy_current_analytical_validation/reference_2d_axisym.py)

---

## 概要

2D 軸対称 FEM（phi = r*A_phi 定式化、z-offset Kelvin 開境界）で、
渦電流を内部まで完全に解いた **Full-Resolution** と、
導体表面の **SIBC (Robin BC)** が L, P ともに 1-2% 以内で一致することを確認した。

**SIBC は Robin 境界条件そのもの**であり、導体内部は解かない（hole approach）。

---

## 形状

IH パネル `ih_fem_kelvin_sample.jou` と同一（gap 5度以外）：

| 要素 | パラメータ |
|------|-----------|
| コイル | 円形断面、R = 30 mm、a = 3 mm |
| ワークピース | 円筒、r = 25 mm、h = 25 mm |
| 開境界 | z-offset Kelvin 変換 (a = 100 mm) |

---

## 検証結果

### Air-only（渦電流なし）

| 手法 | L [nH] | Neumann 解析解との誤差 |
|------|--------|----------------------|
| 2D FEM (Kelvin) | 99.47 | **+0.2%** |
| Neumann 公式 | 99.23 | (参照値) |

### Full-Resolution vs SIBC (Kelvin 付き)

SIBC の Z_s は円筒 Bessel 式（Dowell flat-slab ではない）：

$$Z_s = \rho \gamma \frac{I_1(\gamma a)}{I_0(\gamma a)}, \quad \gamma = \sqrt{j\omega \mu_r \mu_0 \sigma}$$

| 材料 | R/delta 範囲 | L 誤差 | P 誤差 |
|------|-------------|--------|--------|
| **銅** (sigma=5.8e7) | 3.8 - 84.6 | < 0.3% | < 1% |
| **鋼** (sigma=2e6, mu_r=100) | 7.0 - 157 | < 0.4% | < 2% |
| **アルミ** (sigma=3.5e7) | 2.9 - 65.7 | < 0.8% | < 2% |

全材料・全周波数で **L < 1%、P < 2%** で一致。

### 詳細テーブル

```
     Mat      f   R/d  L_full_K  L_sibc_K   ratio    P_full    P_sibc  Pratio
================================================================================
  copper    100   3.8     75.68     75.42  0.9966  4.14e-06  4.15e-06  1.0014
  copper   1000  12.0     61.76     61.69  0.9988  1.93e-05  1.93e-05  0.9961
  copper   7000  31.7     57.18     57.15  0.9995  5.83e-05  5.79e-05  0.9932
  copper  50000  84.6     55.39     55.37  0.9997  1.64e-04  1.63e-04  0.9958
   steel    100   7.0    170.67    169.99  0.9961  1.82e-06  1.84e-06  1.0129
   steel   1000  22.2    157.25    156.67  0.9963  4.67e-05  4.62e-05  0.9899
   steel   7000  58.8    131.91    131.62  0.9977  5.19e-04  5.14e-04  0.9912
   steel  50000 157.1     98.40     98.33  0.9993  3.62e-03  3.60e-03  0.9955
aluminum    100   2.9     81.00     80.34  0.9918  4.55e-06  4.65e-06  1.0199
aluminum   1000   9.3     63.79     63.70  0.9986  2.35e-05  2.35e-05  0.9987
aluminum   7000  24.6     57.99     57.96  0.9994  7.33e-05  7.28e-05  0.9930
aluminum  50000  65.7     55.70     55.68  0.9997  2.09e-04  2.08e-04  0.9950
```

---

## SIBC 実装の要点

### Robin BC（正しい実装）

```python
# 導体はメッシュに含めない（hole）
# 導体表面に Robin BC を課す
a += (1j * omega / Z_s) / r_safe * u * v * ds("wp_bnd")
```

- **導体内部は FEM として解かない**
- 表面インピーダンス Z_s が導体内部の物理を全て記述
- Hole boundary に Robin BC を課すだけ

### Interface approach は間違い

導体内部を空気としてメッシュし、内部界面に Robin BC を課す方法は**間違い**：
磁束が Robin BC を迂回して透明な内部を素通りする。

| approach | Steel 1kHz L | 正解比 |
|----------|-------------|--------|
| **Hole + Robin** | 156.67 nH | **1.00** |
| Interface + Robin | 96.24 nH | 0.62 |

---

## 修正したバグ

### 1. WorkPlane.Arc の中心位置

`WorkPlane.Arc(a, 180)` はデフォルトの heading (+x) のまま使うと、
円の中心が (R+a, a) になる。正しくは (R, 0)。

**修正**: Arc の前に `Direction(0, 1)` を追加。

修正前: L = 113 nH (+14% 過大)、修正後: L = 99.5 nH (Neumann と 0.2% 一致)

### 2. Dodd-Deeds 式の前係数

無限ロッド解析解の前係数: `2*pi*r_c` → `2*pi*r_c**2`

---

## 今後の方針

この 2D 軸対称 Kelvin + hole + Robin BC の手法は 3D にも同様に適用可能。
3D パネル (`calc_fem_kelvin.py`) との定量比較が次のステップ。
