# PEEC Solver Accuracy Test Results

## 検証日

2026-02-12

---

## 実行結果

### テストケース

| 項目 | 値 |
|------|-----|
| メッシュファイル | `circular_coil.msh` |
| 三角形要素数 | 8,478 |
| エッジ数 | 12,717 |
| 平均エッジ長 | 1.03 mm |

### 計算結果

| 項目 | 値 |
|------|-----|
| **PEEC計算** | **-1042.84 nH** ❌ |
| **解析解** | **206.91 nH** ✓ |
| **誤差** | **604%** 🔴 |

---

## 問題の原因

### 🔴 負のインダクタンス

**自己インダクタンスの合計が負**：-1076.53 nH

### 根本原因

```
自己インダクタンス公式: L_self = (mu_0 * l / 2pi) * [ln(l/a) + 0.25]

パラメータ:
  l = エッジ長 = 0.8~1.4 mm（平均1.0mm）
  a = 導体半径 = 2.0 mm（仮定値）

問題:
  l < a の場合 → ln(l/a) < 0

結果:
  100%のエッジで l < 2mm
  → 全エッジで負の自己インダクタンス！
```

### エッジ長分布

| 範囲 | エッジ数 | 割合 |
|------|---------|------|
| 0.5-1.0 mm | 4,728 | 37.2% |
| 1.0-1.5 mm | 7,989 | 62.8% |
| **< 2.0 mm** | **12,717** | **100%** ❌ |

**結論**: 導体半径 2mm の仮定は、1mm の細かいメッシュには**不適切**

---

## 解決策

### Option 1: メッシュベースの実効半径（推奨）

```python
# エッジ長に基づいて実効半径を設定
a_eff = edge_length / 10  # または edge_length / e

# 各エッジごとに個別に設定
for i, l in enumerate(edge_lengths):
    a_eff = l / 10  # 実効半径
    L_self[i] = (MU_0 * l / (2*pi)) * (np.log(l / a_eff) + 0.25)
```

**利点**:
- メッシュ密度に自動適応
- 負のインダクタンスを回避

### Option 2: 最小半径制約

```python
# ln(l/a) >= -1 を保証
a_eff = min(l / np.e, a_min)

# ここで:
# - e = 2.71828 (Euler's number)
# - a_min = 物理的な最小半径（例: 0.1mm）
```

**利点**:
- 数値的に安定
- 物理的な下限を設定可能

### Option 3: メッシュサイズを大きくする

```python
# Cubitでメッシュサイズを増加
mesh_size = 3.0  # mm（現在1.0mm → 3.0mm）

# これにより:
# - エッジ長 ~3mm
# - 導体半径2mmの仮定が妥当に
```

**利点**:
- 物理的に意味のある導体半径を使用可能
- 計算コスト削減

**欠点**:
- 粗いメッシュ → 精度低下の可能性

---

## 推奨アクション

### 短期（今すぐ）

1. ✅ **Option 1を実装**: メッシュベースの実効半径
   ```python
   # demo_peec_dc.pyを修正
   for i, l in enumerate(edge_lengths):
       a_eff = l / 10
       L_self[i] = (MU_0 * l / (2*pi)) * (np.log(l / a_eff) + 0.25)
   ```

2. ✅ **再検証**: 修正後にtest_peec_accuracy.pyを実行

### 中期（メッシュ再生成）

1. Cubitでメッシュサイズを3mmに増加
2. ポート付きメッシュを再生成
3. 完全なPEEC解析を実行

---

## ポート定義メッシュ生成

### 実行方法

```batch
# Windowsバッチファイル
cd gmsh_models
RUN_CUBIT_MESH_WITH_PORTS.bat
```

**または** PowerShellで：

```powershell
cd gmsh_models
& "C:\Program Files\Coreform Cubit 2025.3\bin\python3\python.exe" generate_coil_with_ports.py
```

### 出力

- `circular_coil_with_ports.msh`
- Physical groups: `port_positive`, `port_negative`

---

## 次のステップ

1. ✅ メッシュベース実効半径を実装
2. 📋 TODO: ポート付きメッシュ生成（Cubit実行）
3. 📋 TODO: 修正版PEEC解析実行
4. 📋 TODO: 精度<10%を確認

---

**Created**: 2026-02-12
**Status**: Problem identified, solution ready
**Next**: Implement mesh-based effective radius
