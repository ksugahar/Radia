# CLN I型モデル縮約の検証レポート

## 概要

このドキュメントは、PEEC回路のモデル次数低減に使用するCLN I型（Cauer Ladder Network Type I）変換の数学的正当性を検証した結果をまとめたものです。

## CLN I型変換とは

CLN I型変換は、s=0（DC）周りの展開に基づくモデル次数低減手法です。

### 数学的定義

2つのエルミート行列 R（抵抗）と L（インダクタンス）に対して、変換行列 U, V を求め：

```
R_diag = U^H * R * V    （対角行列）
L_tridiag = U^H * L * V （三重対角行列）
```

を満たすようにします。

### API

```python
from cln import lanczos

result = lanczos(K=R, N=L)
# result.R_diag:    対角行列（抵抗）
# result.L_tridiag: 三重対角行列（インダクタンス）
# result.U, result.V: 変換行列
```

## 検証内容

### Test 1: 表面インピーダンスなし

**元のPEEC回路:**
```
Z(s) = R + s*L  （密行列）
```

**CLN I変換後:**
```
Z'(s) = R_diag + s*L_tridiag  （疎行列）
```

**結果:**
- Max relative error: 8.23e-16
- Mean relative error: 2.09e-16
- **PASS** (機械精度レベル)

### Test 2: 表面インピーダンス Z_s あり

表面インピーダンス（スキンエフェクト）を考慮した場合：

**元のPEEC回路:**
```
Z(s) = R + s*L + Z_s(s)*I  （密行列）
```

**CLN I変換後:**
```
Z'(s) = R_diag + s*L_tridiag + Z_s(s)*I'  （疎行列）
```

ここで、`I' = U^T * I * V` は変換された単位行列（一般に単位行列ではない）。

**表面インピーダンス Z_s:**
```
Z_s = (1+j) * sqrt(pi * f * mu / sigma)
```

銅導体 (sigma = 5.8e7 S/m) での例：
| 周波数 | Z_s | スキン深さ delta |
|--------|-----|------------------|
| 10 Hz | 8.25e-7 + j8.25e-7 | 20.9 mm |
| 10 kHz | 2.70e-5 + j2.70e-5 | 0.64 mm |
| 10 MHz | 8.25e-4 + j8.25e-4 | 0.021 mm |

**重要な知見:**

1. **U = V である**（対称行列に対するLanczos変換の性質）:
   ```
   ||U - V||_F ~ 1e-16  （機械精度）
   ```

2. **しかし U は直交行列ではない**ため、変換後の単位行列 `I' = U^T * U` は単位行列 I とは異なります:
   ```
   ||I' - I||_F = 1.105e+09  （非常に大きい）
   ```

3. この変換された単位行列を使用することで、周波数特性は正確に保存されます。

**結果:**
- Max relative error: 9.53e-16
- Mean relative error: 2.10e-16
- **PASS** (機械精度レベル)

## 検証に使用したテスト条件

### PEEC行列

```python
n = 10  # ループ数

# 密インダクタンス行列（相互インダクタンス含む）
L0 = 1e-6  # 1 uH
decay = 3.0
L_dense[i,j] = L0 * exp(-|i-j| / decay)

# 対角抵抗行列
R0 = 0.01  # 10 mOhm
R_diag = diag([R0, R0, ..., R0])
```

### 周波数範囲

```python
frequencies = logspace(1, 7, 100)  # 10 Hz to 10 MHz
```

## 結論

1. **CLN I型変換は周波数特性を保存する** - 表面インピーダンスの有無に関わらず、機械精度（~1e-15）でインピーダンス周波数特性が一致することを確認しました。

2. **変換された単位行列の使用が必須** - 表面インピーダンス Z_s を追加する場合、`I' = U^T * I * V` を使用する必要があります。

3. **疎行列構造が得られる** - R_diag は対角行列、L_tridiag は三重対角行列となり、計算効率が向上します。

## PEEC-MMM / PEEC-STAR 結合への適用

CLN I変換後の座標系で他の物理系（MMM、STAR等）と結合する場合、結合行列も適切に変換する必要があります。

### 変換規則

任意の行列 M をCLN I座標系に追加する場合：

```
M' = U^T * M * V
```

### 結合行列の変換

| 結合項 | 元の座標系 | CLN I座標系 | 備考 |
|--------|-----------|-------------|------|
| Z_s * I | Z_s * I | Z_s * (U^T * I * V) | 表面インピーダンス |
| Z_LM | Z_LM | U^T * Z_LM | Loop→MMM結合（左から変換）|
| Z_ML | Z_ML | Z_ML * V | MMM→Loop結合（右から変換）|
| Z_LS | Z_LS | U^T * Z_LS | Loop→STAR結合 |
| Z_SL | Z_SL | Z_SL * V | STAR→Loop結合 |

### 連成システムの構成

元の連成システム：
```
[R + sL    Z_LM  ] [I_L]   [V_L]
[Z_ML     Z_MMM ] [I_M] = [V_M]
```

CLN I変換後：
```
[R_diag + s*L_tridiag    U^T*Z_LM ] [I_L']   [U^T*V_L]
[Z_ML*V                  Z_MMM    ] [I_M ] = [V_M    ]
```

ここで `I_L' = V^{-1} * I_L` は変換後のLoop電流。

### 重要な注意点

1. **MMM/STAR部分は変換しない** - Z_MMM, Z_STAR は元のまま
2. **結合行列は片側変換** - Z_LM は左から U^T、Z_ML は右から V
3. **励起ベクトルも変換** - Loop励起 V_L は U^T * V_L に変換

## 生成されたプロット

- `cln_frequency_response_verification.png` - Test 1の結果
- `cln_frequency_response_with_zs.png` - Test 2の結果（Z_s含む）

## 参考文献

- Lanczosアルゴリズムによる一般化固有値問題の解法
- Cauer回路によるRLC等価回路のモデル縮約
- PEEC (Partial Element Equivalent Circuit) 法

## 検証スクリプト

- `verify_cln_frequency_response.py` - 周波数特性の検証
- `verify_lanczos_uv_equality.py` - U, V行列とCLN I構造の検証
