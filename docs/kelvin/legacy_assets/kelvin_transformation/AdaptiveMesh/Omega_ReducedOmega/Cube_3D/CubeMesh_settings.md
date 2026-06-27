# Archived CubeMesh.py 設定ドキュメント

## 概要

CubeMesh.pyは、3次元立方体モデルのメッシュ生成クラスでした。Kelvin変換を用いた無限領域問題にも対応しています。スタンドアロン Python は `examples` から削除済みで、完全なソースは `docs/kelvin/kelvin_remaining_examples_archive_results.json` に保存されています。

## デフォルトパラメータ

```python
default_values = {
    "name": "Cube",
    "mur": 1000,           # 比透磁率
    "msize": meshsize.moderate,  # メッシュサイズ
    "ndiv": 5,             # 分割数（未使用）
    "type": 0,             # 定式化タイプ (0: Omega法, 1: 全ポテンシャル法)
    "curveOrder": 1,       # 曲線要素次数
    "outerBox": 5,         # 外部ボックスサイズ（Kelvin無効時）
    "rKelvin": 0           # Kelvin変換半径（0=無効）
}
```

## ジオメトリ構成

### 1. 鉄領域 (iron)

- **形状**: Box((0,0,0), (1,1,1)) - 単位立方体
- **材料名**: `iron`
- **境界条件**:
  - X=0面: `Bn0` (法線磁束密度=0)
  - Y=0面: `Bn0`
  - Z=0面: `Ht0` (接線磁界=0)
  - type=1の場合: X=1, Y=1, Z=1面に `A_Omega_boundary`

### 2. A領域 (A_domain)

- **形状**: Box((0,0,0), (1.5,1.5,1.5))
- **材料名**:
  - type=0: `A_domain`
  - type=1: `Omega_domain`
- **境界条件**:
  - X=0面: `Bn0`
  - Y=0面: `Bn0`
  - Z=0面: `Ht0`
  - type=0の場合: X=1.5, Y=1.5, Z=1.5面に `A_Omega_boundary`

### 3. Omega領域 (Omega_domain)

#### rKelvin = 0 の場合（Kelvin変換なし）
- **形状**: Box((0,0,0), (outerBox, outerBox, outerBox))
- **材料名**: `Omega_domain`
- **境界条件**:
  - X=0面: `Bn0`
  - Y=0面: `Bn0`
  - Z=0面: `Ht0`
  - X=outerBox, Y=outerBox, Z=outerBox面: `Omega0`（外部境界）

#### rKelvin > 0 の場合（Kelvin変換あり）
- **形状**: Sphere(原点, r=rKelvin) ∩ Box((0,0,0), (rk,rk,rk)) - 1/8球
- **材料名**: `Omega_domain`
- **境界条件**:
  - X=0面: `Bn0`
  - Y=0面: `Bn0`
  - Z=0面: `Ht0`

### 4. 外部領域 - Kelvin変換 (external_domain)

rKelvin > 0 の場合のみ生成:

- **形状**: Sphere(center=(2*rk,0,0), r=rk) ∩ Box((center,0,0), (center+rk,rk,rk))
- **材料名**: `Kelvin`
- **境界条件**:
  - X=center面: `Bn0`
  - Y=0面: `Bn0`
  - Z=0面: `Ht0`
- **周期境界**: external_domain.faces[0] ↔ Omega_domain.faces[0] (`PERIODIC`)

## 定式化タイプ (type)

### type = 0 (Omega-Reduced Omega法)
- **reduced_region**: `Omega_domain`
- **total_region**: `iron|A_domain`

### type = 1 (全ポテンシャル法)
- **reduced_region**: `A_domain|Omega_domain`
- **total_region**: `iron`

## 物性値

```python
mu0 = 4.e-7 * math.pi  # 真空透磁率
mu = mu0 * mur         # 鉄の透磁率

mu_d = {
    "iron": mu,
    "A_domain": mu0,
    "Omega_domain": mu0,
    "Kelvin": mu0,
    'default': mu0
}
```

## 境界名一覧

| 境界名 | 説明 |
|--------|------|
| `Bn0` | 法線磁束密度 = 0 (対称境界) |
| `Ht0` | 接線磁界 = 0 (対称境界) |
| `A_Omega_boundary` | A法とOmega法の接続境界 |
| `Omega0` | 外部境界（Kelvin無効時のみ） |

## 対称性

1/8モデル（オクタント対称性）:
- X ≥ 0
- Y ≥ 0
- Z ≥ 0

## コイル設定

```python
self.coil = EMPY_UNIF(0, 0, 1, 0)  # Z方向一様磁界
```

## 使用例

```python
# 基本的な使用法
cube = CubeMesh()

# カスタムパラメータ
cube = CubeMesh(mur=500, rKelvin=2.0, curveOrder=2)

# メッシュ情報の表示
cube.Print()
```
