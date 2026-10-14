# PEEC Port Definition Guide

PEECインダクタンス計算のためのポート指定方法

---

## ポート指定の2つのアプローチ

### 方法1: メッシュ生成時に定義（推奨）

**メリット**:
- OK 幾何学的に正確（位置を明示的に指定）
- OK 再現性が高い（メッシュファイルにポート情報が含まれる）
- OK ポート領域の可視化が容易（GMSH GUIで確認可能）

**ワークフロー**:
```
Cubit: ポート領域をSidesetとして定義
  ↓
GMSH v4.1 export: Physical groupとして保存
  ↓
Python: Physical groupからポート要素を抽出
  ↓
PEEC: ポートDOF（エッジ）を特定
```

**実装例**: `generate_coil_with_ports.py`

---

### 方法2: メッシュ読込後にPythonで定義

**メリット**:
- OK 柔軟性が高い（メッシュ再生成不要）
- OK パラメータスタディに適している
- OK 複数ポート配置の試行が容易

**ワークフロー**:
```
GMSH mesh読込
  ↓
Python: 座標条件でポート要素を検索
  ↓
PEEC: ポートDOFを特定
```

**実装例**: 座標ベース検索（後述）

---

## 方法1: Cubitでポート定義（推奨）

### Step 1: Cubitでポート領域を定義

```python
# PORT 1: phi=0付近（正端子、+X方向）
cubit.cmd("sideset 1 add tri in surface with x_coord > 45")
cubit.cmd("sideset 1 name 'port_positive'")

# PORT 2: phi=180付近（負端子、-X方向）
cubit.cmd("sideset 2 add tri in surface with x_coord < -45")
cubit.cmd("sideset 2 name 'port_negative'")
```

### Step 2: GMSH v4.1形式でエクスポート

```python
cubit.cmd('export gmsh "coil_with_ports.msh" overwrite')
```

**出力**: Physical groupsとしてポート情報が保存される

### Step 3: Pythonでポート要素を抽出

```python
import gmsh

gmsh.initialize()
gmsh.open("coil_with_ports.msh")

# Physical groupsを取得
phys_groups = gmsh.model.getPhysicalGroups()
for dim, tag in phys_groups:
    name = gmsh.model.getPhysicalName(dim, tag)
    if 'port' in name.lower():
        # ポート要素を取得
        entities = gmsh.model.getEntitiesForPhysicalGroup(dim, tag)
        # ...
```

**詳細**: `demo_port_extraction.py`参照

---

## 方法2: Pythonで座標ベース検索

### 例: 円形コイルのポート（phi=0, phi=180）

```python
import numpy as np

# メッシュ読込後
node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
coords = node_coords.reshape(-1, 3)

# ポート1: +X方向（phi=0）
port1_nodes = []
for i, coord in enumerate(coords):
    x, y, z = coord
    # 条件: x > 45mm かつ mean_radius付近
    r = np.sqrt(x**2 + y**2)
    if x > 45 and abs(r - 50.0) < 5.0:  # mean_radius=50mm
        port1_nodes.append(node_tags[i])

# ポート2: -X方向（phi=180）
port2_nodes = []
for i, coord in enumerate(coords):
    x, y, z = coord
    r = np.sqrt(x**2 + y**2)
    if x < -45 and abs(r - 50.0) < 5.0:
        port2_nodes.append(node_tags[i])
```

**注意**: この方法は幾何形状に依存するため、汎用性は低い

---

## PEECでのポート処理

### ポートの自由度（DOF）

PEECでは**エッジ（辺）**が電流の自由度となる：

| メッシュ要素 | PEEC DOF |
|-------------|----------|
| ノード（頂点） | 電位（電荷） |
| **エッジ（辺）** | **電流** |
| パネル（面） | 電荷密度 |

**ポート = 電流が流入/流出するエッジ集合**

### ポートエッジの抽出

```python
# 三角形要素からエッジを生成
def get_edges_from_triangle(nodes):
    """三角形の3エッジを取得"""
    n0, n1, n2 = nodes
    edge1 = tuple(sorted([n0, n1]))
    edge2 = tuple(sorted([n1, n2]))
    edge3 = tuple(sorted([n2, n0]))
    return [edge1, edge2, edge3]

# ポート領域の全エッジを収集
port_edges = set()
for elem_tag in port_elements:
    elem_type, elem_nodes = gmsh.model.mesh.getElement(elem_tag)
    edges = get_edges_from_triangle(elem_nodes)
    port_edges.update(edges)
```

### インピーダンス抽出

```python
# PEEC行列方程式: [Z(f)] * I = V
# Z(f) = R(f) + j*omega*L  # SIBC抵抗 + 誘導性リアクタンス

# ポート電圧印加
V = np.zeros(n_edges, dtype=complex)
V[port_positive_edges] = +1.0  # 正端子に+1V
V[port_negative_edges] = 0.0   # 負端子をGND

# 電流を解く
I = solve(Z_matrix, V)

# ポートインピーダンス
I_port = np.sum(I[port_positive_edges])  # 正端子流入電流
Z_port = 1.0 / I_port  # V=1.0を印加したので Z = V/I = 1/I
```

---

## 実装例

### 1. メッシュ生成（ポート定義付き）

```bash
cd gmsh_models
"${CUBIT_PATH:-<Coreform Cubit 2025.8+>/bin}/python3/python.exe" generate_coil_with_ports.py
```

**出力**: `circular_coil_with_ports.msh`

### 2. ポート要素抽出

```bash
cd ..
python demo_port_extraction.py
```

**出力**: ポートエッジリスト、重心座標

### 3. PEEC解析（TODO）

```bash
python demo_peec_impedance.py
```

**出力**: Z(f)、L(f)、R(f)

---

## GMSH GUIでポート確認

```bash
gmsh circular_coil_with_ports.msh
```

**操作**:
1. Tools -> Visibility -> Physical groups
2. "port_positive"と"port_negative"をチェック
3. 異なる色で表示される

---

## まとめ

| 項目 | 方法1: Cubit定義 | 方法2: Python検索 |
|------|-----------------|-----------------|
| 推奨度 | OK 推奨 | △ 場合による |
| 再現性 | 高い | 低い（座標依存） |
| 柔軟性 | 低い（再生成必要） | 高い |
| 可視化 | 容易（GMSH GUI） | 困難 |
| 用途 | 本番解析 | パラメータスタディ |

**推奨**: 本番解析では**方法1（Cubit定義）**を使用し、探索段階では**方法2（Python検索）**を併用

---

**Created**: 2026-02-12
**Author**: Radia Development Team
