# Files in ngsolve_integration

## Main Files

### Python Scripts

1. **visualize_field.py** - メインの可視化スクリプト
   - CoefficientFunctionを直接使用（GridFunction補間なし）
   - 正確なRadia値を取得（誤差0%）
   - VTKエクスポート機能

2. **simple_example.py** - シンプルな使用例
   - rad_ngsolveの基本的な使い方
   - 最小限のコード例

3. **example_rad_ngsolve_demo.py** - デモスクリプト
   - より詳しい使用例
   - 複数のテストケース

### Test Scripts

4. **test_units_verification.py** - 単位変換の検証
   - NGSolve (m) ↔ Radia (mm)の変換をテスト
   - すべてのテストポイントでRadiaと完全一致を確認

5. **test_complete_verification.py** - 完全な検証
   - CoefficientFunctionの直接評価（正確）
   - GridFunction補間（FEM誤差あり）
   - フィールド統計

### Documentation

6. **README.md** - プロジェクト概要
7. **STATUS.md** - 実装ステータス
8. **EXAMPLES_GUIDE.md** - 使用例ガイド
9. **UNIT_CONVERSION.md** - 単位変換の説明
10. **COMPLETE_FIX_SUMMARY.md** - 修正の完全な要約

### Other

11. **ngsolve.tcl** - NGSolve設定ファイル

## Removed Files

以下のファイルは不要として削除されました：

### /src/python から削除
- `rad_ngsolve.cpp.backup` - バックアップ（不要）
- `rad_ngsolve_gil_safe.cpp` - 中間バージョン（不要）
- `rad_ngsolve_no_units.cpp` - 単位変換なしバージョン（不要）
- `rad_ngsolve_python_callback.cpp` - 中間バージョン（不要）
- `rad_ngsolve_v2.cpp` - バージョン2（不要）

### examples/ngsolve_integration から削除
- `test_cf_direct.py` - デバッグ用テスト（不要）
- `test_coordinates.py` - デバッグ用（不要）
- `test_rad_ngsolve_full.py` - 古いテスト（不要）
- `test_radfld_call.py` - デバッグ用（不要）
- `test_radfld_debug.py` - デバッグ用（不要）
- `test_radfld_direct.py` - デバッグ用（不要）
- `test_simple_cf.py` - デバッグ用（不要）
- `test_cf_vs_gf.py` - デバッグ用（不要）
- `visualize_field_2d.py.backup` - 2Dバージョンバックアップ（不要）
- `visualize_field_2d_backup.py` - 2Dバージョンバックアップ（不要）
- `visualize_field_3d.py` - 中間バージョン（不要）
- `FIX_SUMMARY.md` - 部分的なドキュメント（COMPLETE_FIX_SUMMARY.mdに統合）
- `RADFLD_FIX_NOTES.md` - 部分的なドキュメント（COMPLETE_FIX_SUMMARY.mdに統合）
- `VISUALIZE_FIELD_FIX.md` - 部分的なドキュメント（COMPLETE_FIX_SUMMARY.mdに統合）
- `radia_field.vtu` - 古いVTK出力（GridFunction版、誤差あり）
- `radia_field_3d.vtu` - 古いVTK出力（不要）

## Key Points

### 重要な変更点

1. **CoefficientFunctionを直接使用**
   - `GridFunction.Set()`は補間誤差が大きい（60%の誤差）
   - CoefficientFunctionの直接評価は正確（誤差0%）

2. **単位変換は自動**
   - NGSolve: メートル (m)
   - Radia: ミリメートル (mm)
   - `rad_ngsolve.cpp`で自動変換（x1000）

3. **3Dメッシュを使用**
   - 2Dメッシュはz座標を無視する
   - 3D磁場評価には3Dメッシュが必要

### 推奨される使用方法

```python
import rad_ngsolve

# Radiaオブジェクト作成（mm単位）
magnet = rad.ObjRecMag([0, 0, 0], [20, 20, 30], [0, 0, 1.2])

# NGSolveメッシュ作成（m単位）
geo = CSGeometry()
geo.Add(OrthoBrick(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05)))
mesh = Mesh(geo.GenerateMesh(maxh=0.02))

# CoefficientFunction作成（自動単位変換）
B_cf = rad_ngsolve.RadBfield(magnet)

# 直接評価（正確）
mesh_pt = mesh(0, 0, 0.02)  # 20mm（m単位で指定）
B = B_cf(mesh_pt)  # 誤差0%

# VTKエクスポート（正確な値）
vtk = VTKOutput(mesh, coefs=[B_cf], names=['B_field'], filename="radia_field")
vtk.Do()
```

### 避けるべき方法

```python
# ✗ GridFunction補間を使用（大きな誤差）
gfB = GridFunction(VectorH1(mesh, order=2))
gfB.Set(B_cf)  # 補間により最大60%の誤差
B = gfB(mesh_pt)  # 不正確
```

Date: 2025-10-31
