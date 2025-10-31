# NGSolve Integration Status

## 実装状況

### 完了項目 ✓

1. **ソースコード実装**
   - `src/python/rad_ngsolve.cpp` - 完全実装
   - `RadBfield`, `RadHfield`, `RadAfield` クラス
   - NGSolve CoefficientFunction インターフェース実装

2. **CMake設定**
   - NGSolve自動検出
   - `add_ngsolve_python_module` 使用
   - すべての依存関係設定完了

3. **ビルド**
   - `rad_ngsolve.pyd` ビルド成功
   - 場所: `build/Release/rad_ngsolve.pyd`
   - サイズ: 1.3 MB

4. **ドキュメント**
   - ユーザーガイド作成
   - 技術ドキュメント作成
   - テストスクリプト作成

### 既知の問題 ⚠

**DLL依存関係の問題**

`rad_ngsolve.pyd`をインポートする際、以下のエラーが発生：
```
ImportError: DLL load failed while importing rad_ngsolve: 指定されたモジュールが見つかりません。
```

**原因:**
- `rad_ngsolve.pyd`が`ngsolve.dll`, `libngsolve.dll`, `ngcore.dll`などに依存
- これらのDLLがWindows DLL検索パスに含まれていない

**影響:**
- モジュールはビルド済みだが実行時にロード不可
- NGSolve本体のPythonパッケージ内のDLLを見つけられない

## 解決策

### 方法1: 環境変数設定（推奨）

DLLパスをシステムPATHに追加：

```powershell
# PowerShellで実行
$env:PATH = "C:\Program Files\Python312\Lib\site-packages\ngsolve;$env:PATH"
$env:PATH = "C:\Program Files\Python312\Lib\site-packages\netgen;$env:PATH"
```

その後Pythonを起動：
```python
import sys
sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
import rad_ngsolve  # これで動作するはず
```

### 方法2: DLLディレクトリを動的に追加

Pythonスクリプト内で：
```python
import os
import sys

# Windows 10/11の場合
if hasattr(os, 'add_dll_directory'):
    os.add_dll_directory(r'C:\Program Files\Python312\Lib\site-packages\ngsolve')
    os.add_dll_directory(r'C:\Program Files\Python312\Lib\site-packages\netgen')

sys.path.insert(0, r'S:\radia\01_GitHub\build\Release')
import rad_ngsolve
```

### 方法3: 依存DLLをコピー

必要なDLLを`build/Release`にコピー：
```powershell
Copy-Item "C:\Program Files\Python312\Lib\site-packages\ngsolve\*.dll" `
          "S:\radia\01_GitHub\build\Release\"
Copy-Item "C:\Program Files\Python312\Lib\site-packages\netgen\*.dll" `
          "S:\radia\01_GitHub\build\Release\"
```

### 方法4: NGSolveと同じ場所にインストール

最も確実な方法：
```powershell
# rad_ngsolve.pydをNGSolveのsite-packagesにコピー
Copy-Item "S:\radia\01_GitHub\build\Release\rad_ngsolve.pyd" `
          "C:\Program Files\Python312\Lib\site-packages\"
```

その後：
```python
import rad_ngsolve  # どこからでもimport可能
```

## 実装の技術詳細

### アーキテクチャ

```
Python Script
    ↓
rad_ngsolve.pyd (Pythonモジュール)
    ↓
ngfem::RadiaBFieldCF (CoefficientFunction)
    ↓
RadFld() (Radia C API)
    ↓
Radia Core (磁場計算)
```

### 関数シグネチャ

```python
# Python側API
rad_ngsolve.RadBfield(radia_obj: int) -> CoefficientFunction
rad_ngsolve.RadHfield(radia_obj: int) -> CoefficientFunction
rad_ngsolve.RadAfield(radia_obj: int) -> CoefficientFunction
```

### C++実装の要点

```cpp
namespace ngfem {
    class RadiaBFieldCF : public CoefficientFunction {
        int radia_obj;  // Radiaオブジェクトインデックス

        virtual void Evaluate(
            const BaseMappedIntegrationPoint& mip,
            FlatVector<> result) const override
        {
            auto pnt = mip.GetPoint();
            double coords[3] = {pnt[0], pnt[1], pnt[2]};
            double B[3];
            RadFld(B, &nB, radia_obj, "b", coords, 1);
            result(0) = B[0];
            result(1) = B[1];
            result(2) = B[2];
        }
    };
}
```

## テスト

### 基本テスト（DLL問題解決後）

```python
import radia as rad
import rad_ngsolve

# Radiaジオメトリ
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1000])
rad.Solve(magnet, 0.0001, 10000)

# CoefficientFunction作成
B = rad_ngsolve.RadBfield(magnet)
print(type(B))  # <class 'ngsolve.fem.CoefficientFunction'>
```

### NGSolve統合テスト

```python
from ngsolve import *
from netgen.csg import *

# メッシュ
geo = CSGeometry()
box = OrthoBrick(Pnt(-50,-50,-50), Pnt(50,50,50))
geo.Add(box)
mesh = Mesh(geo.GenerateMesh(maxh=10))

# Radia磁場をNGSolveで使用
B_integral = Integrate(B, mesh)
Draw(B, mesh, "B_field")
```

## ビルド方法

```powershell
# 1. NGSolveとNetgenをインストール
conda install -c ngsolve ngsolve

# 2. ビルドスクリプト実行
.\build_ngsolve.ps1

# 3. 出力確認
ls build\Release\rad_ngsolve.pyd
```

## まとめ

**実装**: ✓ 完了
**ビルド**: ✓ 成功
**実行**: ⚠ DLL依存関係の設定が必要

DLL問題を解決すれば、Radiaの磁場計算をNGSolveのFEMシミュレーションで直接使用可能になります。

## 参考リンク

- [NGSolve Documentation](https://ngsolve.org/)
- [Radia Manual](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- [Windows DLL Search Order](https://learn.microsoft.com/en-us/windows/win32/dlls/dynamic-link-library-search-order)
