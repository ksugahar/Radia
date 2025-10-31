# NGSolve Integration Usage Guide

## 概要

RadiaとNGSolveの統合には2つの実装があります：

| 実装 | 状態 | 推奨度 |
|------|------|--------|
| **Pure Python版** (`rad_ngsolve_py.py`) | ✓ 動作中 | ⭐⭐⭐ 推奨 |
| **C++版** (`rad_ngsolve.pyd`) | △ DLL問題 | Conda環境のみ |

## Pure Python版の使い方（推奨）

### 基本的な使用方法

```python
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\src\python")
sys.path.insert(0, r"S:\radia\01_GitHub\build\lib\Release")

import radia as rad
import rad_ngsolve_py as rad_ngsolve
from ngsolve import *
from netgen.csg import *

# Radiaジオメトリを作成
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1])

# Radiaフィールドオブジェクトを作成
B_radia = rad_ngsolve.RadBfield(magnet)

# メッシュを生成
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=15))

# GridFunctionを作成して可視化
fes = H1(mesh, order=1)
gf_Bz = GridFunction(fes)

# メッシュ頂点でフィールドを評価
for i in range(mesh.nv):
    pnt = mesh.vertices[i].point
    B_vec = B_radia.Evaluate(pnt[0], pnt[1], pnt[2])
    gf_Bz.vec[i] = B_vec[2]  # Z成分

# 可視化
import netgen.gui
netgen.gui.Draw(gf_Bz, mesh, "Bz")

# 積分
integral = Integrate(gf_Bz, mesh)
```

### 完全な例

`examples/ngsolve_integration/visualize_field.py`を参照してください。

## C++版の使い方（Conda環境）

C++版はEMPY_Fieldと同じパターンで実装されていますが、DLL依存関係の問題があります。

### 前提条件

```bash
# Conda環境を作成
conda create -n radia_ngsolve python=3.12
conda activate radia_ngsolve
conda install -c ngsolve ngsolve
```

### ビルド

```powershell
cd S:\radia\01_GitHub
.\build_ngsolve.ps1
```

### 使用方法

```python
# Conda環境内で実行
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\build\Release")

import radia as rad
import rad_ngsolve  # C++版

magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1])
B = rad_ngsolve.RadBfield(magnet)

# NGSolveで直接使用可能（理論上）
from ngsolve import *
Bz = B[2]
```

**注意**: C++版は現在DLL問題により動作しません。Pure Python版の使用を推奨します。

## 実装の比較

### Pure Python版の特徴

✓ **利点**:
- DLL依存なし
- インストール不要
- デバッグが容易
- クロスプラットフォーム

✗ **欠点**:
- GridFunctionへの手動変換が必要
- C++版より若干遅い（実用上は問題なし）

### C++版の特徴

✓ **利点**:
- NGSolve CoefficientFunctionとして直接使用可能
- 若干高速
- EMPY_Fieldと同じパターン

✗ **欠点**:
- DLL依存関係の問題
- Conda環境が必要
- ビルドが複雑

## 技術的詳細

### Pure Python版の仕組み

`rad_ngsolve_py.py`は以下のように動作します：

1. `RadiaBFieldCF`クラスが`Evaluate(x, y, z)`メソッドを提供
2. `rad.Fld()`を呼び出してRadiaフィールドを計算
3. ユーザーがGridFunctionを作成し、メッシュ頂点で評価
4. NGSolveの`Draw()`や`Integrate()`で使用

### C++版の仕組み

`rad_ngsolve.cpp`は以下のように実装されています：

```cpp
class RadiaBFieldCF : public CoefficientFunction
{
public:
    int radia_obj;

    RadiaBFieldCF(int obj)
        : CoefficientFunction(3), radia_obj(obj) {}

    virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
                         FlatVector<> result) const override
    {
        auto pnt = mip.GetPoint();
        double coords[3] = {pnt[0], pnt[1], pnt[2]};
        double B[3];
        int nB = 3;

        RadFld(B, &nB, radia_obj, "b", coords, 1);

        result(0) = B[0];
        result(1) = B[1];
        result(2) = B[2];
    }
};
```

これはEMPY_Fieldの`B_MAGNET_CF`と同じパターンです。

## トラブルシューティング

### Pure Python版

**問題**: `No module named 'rad_ngsolve_py'`

**解決**:
```python
import sys
sys.path.insert(0, r"S:\radia\01_GitHub\src\python")
```

**問題**: `No module named 'radia'`

**解決**:
```python
sys.path.insert(0, r"S:\radia\01_GitHub\build\lib\Release")
```

### C++版

**問題**: `DLL load failed while importing rad_ngsolve`

**解決**: Pure Python版を使用するか、Conda環境を使用してください。

## サンプルコード

### 例1: 単純な磁場評価

```python
import rad_ngsolve_py as rad_ngsolve

magnet = rad.ObjRecMag([0, 0, 0], [10, 10, 10], [0, 0, 1])
B = rad_ngsolve.RadBfield(magnet)

# 点での評価
Bx, By, Bz = B.Evaluate(0, 0, 20)
print(f"B at (0,0,20): Bz = {Bz:.6f} T")
```

### 例2: NGSolve統合

```python
# メッシュ生成
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=15))

# GridFunction作成
fes = H1(mesh, order=1)
gf_Bz = GridFunction(fes)

# フィールド評価
for i in range(mesh.nv):
    pnt = mesh.vertices[i].point
    B_vec = B.Evaluate(pnt[0], pnt[1], pnt[2])
    gf_Bz.vec[i] = B_vec[2]

# 積分
integral = Integrate(gf_Bz, mesh)
print(f"∫Bz dV = {integral:.6f} T·mm³")
```

### 例3: 可視化

```python
import netgen.gui

# Bz成分の可視化
netgen.gui.Draw(gf_Bz, mesh, "Bz_component")
```

## 参考資料

- [NGSOLVE_PYTHON_SOLUTION.md](NGSOLVE_PYTHON_SOLUTION.md) - Pure Python版の詳細
- [NGSOLVE_DLL_ISSUE.md](NGSOLVE_DLL_ISSUE.md) - C++版のDLL問題
- [EMPY_Field実装](S:/radia/02_EMPY_Field) - C++版の参考実装
- [NGSolve Documentation](https://ngsolve.org/)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
