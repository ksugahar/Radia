# Radia-NGSolve Viewer Comparison

## Executive Summary

**推奨ワークフロー**:
- **開発・デバッグ**: PyVista (デフォルト) - 通常の.pyスクリプトから使用
- **論文・プレゼン**: ParaView (高品質)
- **インタラクティブ探索**: NGSolve webgui - ブラウザベース、.pyスクリプトから使用
- **メッシュ確認**: Netgen GUI

**重要**: Jupyterノートブック（.ipynb）はClaude Codeとの相性が悪いため、通常のPythonスクリプト（.py）を推奨します。PyVistaもwebguiも.pyスクリプトから完全に使用可能です。

---

## Jupyter Notebook vs Python Script

### Python Script（.py）推奨の理由

| 観点 | Python Script (.py) | Jupyter Notebook (.ipynb) |
|------|-------------------|-------------------------|
| **Claude Code編集** | ✅ 容易（通常のテキスト） | ❌ 困難（JSON構造） |
| **バージョン管理** | ✅ Git差分が明確 | ❌ JSON差分が読みにくい |
| **ファイルサイズ** | ✅ 小さい | ❌ 実行結果で肥大化 |
| **デバッグ** | ✅ 標準デバッガ使用可 | ❌ セル単位のみ |
| **自動化** | ✅ CI/CDで実行容易 | ❌ 追加設定必要 |
| **PyVista** | ✅ インタラクティブウィンドウ | ⚠️ インライン表示（制約あり） |
| **webgui** | ✅ ブラウザ自動起動 | ⚠️ ノートブック内表示 |

**結論**: 特別な理由がない限り、通常のPythonスクリプト（.py）を使用してください。

### PyVistaとwebguiは.pyスクリプトで完全動作

```python
# demo.py - 通常のPythonスクリプト
import pyvista as pv

# PyVista: インタラクティブウィンドウが開く（Jupyter不要）
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')  # ← ウィンドウが開く
```

```python
# webgui_demo.py - 通常のPythonスクリプト
from ngsolve.webgui import Draw

# NGSolve webgui: ブラウザが自動的に開く（Jupyter不要）
Draw(B_gf, mesh, 'B_field')  # ← ブラウザタブが開く
```

---

## 詳細比較

### PyVista (推奨デフォルト)

**長所**:
- ✅ Pythonネイティブ - スクリプト統合が容易
- ✅ Jupyter Notebook/Lab完全対応
- ✅ 迅速な可視化（開発イテレーション高速）
- ✅ VTK全機能へのPythonicアクセス
- ✅ アニメーション・GIF出力が簡単
- ✅ インタラクティブウィジェット（スライダー、チェックボックス）
- ✅ ヘッドレス実行可能（CI/CDパイプライン）

**短所**:
- ❌ 論文品質の微調整はParaViewより劣る
- ❌ ベクトルグラフィックス出力は限定的
- ❌ 大規模データ（>10M cells）では遅い

**Radia-NGSolve統合での使用例**:
```python
import pyvista as pv
import radia as rad
from radia_ngsolve import RadiaField

rad.FldUnits('m')

# RadiaフィールドをVTS出力
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)

# PyVistaで可視化
grid = pv.read('field.vts')
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='B_magnitude', cmap='coolwarm', opacity=0.8)
plotter.add_arrows(grid.points, grid['B_field'], mag=0.01, color='black')
plotter.show()

# Jupyter統合
grid.plot(scalars='B_magnitude', cmap='coolwarm', jupyter_backend='static')
```

**ベストプラクティス**:
- 開発中の可視化確認
- パラメータスタディの自動化
- バッチ処理（複数ケースの比較）
- CI/CDでの回帰テスト可視化

---

### ParaView (論文品質)

**長所**:
- ✅ 最高品質のレンダリング
- ✅ ベクトルグラフィックス出力（SVG, PDF）
- ✅ 高解像度ラスタ画像（300+ DPI）
- ✅ 複雑なフィルタチェーン（Glyph, Contour, StreamTracer）
- ✅ アニメーション・キーフレーム
- ✅ 大規模データ対応（分散並列）
- ✅ カメラ・照明の完全制御

**短所**:
- ❌ GUI操作が必要（スクリプト化は可能だが複雑）
- ❌ 迅速なイテレーションには不向き
- ❌ Jupyter統合は限定的（pvpythonは別プロセス）

**Radia-NGSolve統合での使用例**:
```bash
# 1. Radiaフィールド → VTS出力（Python）
python generate_field_vts.py

# 2. ParaViewで開く
paraview field.vts

# 3. ParaView GUI操作:
#    - Glyph filter: ベクトル矢印
#    - Contour: 等値面
#    - Slice: 断面
#    - Camera: アングル調整
#    - Lighting: 照明設定

# 4. 高解像度エクスポート:
#    File > Save Screenshot
#    - Resolution: 3000x2000 (300 DPI at 10x6.67 inch)
#    - Format: PNG (raster) or SVG (vector)
```

**ParaViewスクリプト自動化** (pvpython):
```python
# publication_figure.py
from paraview.simple import *

# Load VTS
reader = XMLStructuredGridReader(FileName=['field.vts'])

# Add glyph
glyph = Glyph(Input=reader, GlyphType='Arrow')
glyph.ScaleFactor = 0.01
glyph.GlyphMode = 'Every Nth Point'
glyph.Stride = 2

# Render
Show(glyph)
view = GetActiveView()
view.CameraPosition = [0.3, 0.2, 0.5]
view.CameraFocalPoint = [0, 0, 0]

# Save high-res image
SaveScreenshot('figure.png', view, ImageResolution=[3000, 2000])
```

**ベストプラクティス**:
- 論文投稿用図表
- プレゼンテーション資料
- 高解像度ポスター
- アニメーション動画（MP4, AVI）

---

### NGSolve webgui (インタラクティブ探索)

**長所**:
- ✅ NGSolve完全統合
- ✅ メッシュ + フィールド同時表示
- ✅ WebGL - ブラウザで動作
- ✅ Jupyter統合（同じノートブック内）
- ✅ Radia CoefficientFunctionと連携
- ✅ リアルタイム更新（パラメータ変更）

**短所**:
- ❌ VTSファイル直接読込不可（NGSolve GridFunctionのみ）
- ❌ 高度なフィルタ機能なし
- ❌ 論文品質の微調整困難
- ❌ エクスポート形式が限定的

**Radia-NGSolve統合での使用例**:
```python
from ngsolve import *
from ngsolve.webgui import Draw
from radia_ngsolve import RadiaField
import radia as rad

rad.FldUnits('m')

# Cubitメッシュ → NGSolve
mesh = Mesh('model.msh')

# Radia磁石
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# Radia CoefficientFunction
B_cf = RadiaField(magnet, 'b')

# GridFunctionに投影
fes = HDiv(mesh, order=2)
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

# インタラクティブ表示（Jupyter内）
Draw(B_gf, mesh, name='B_field', vectors={'grid_size': 10})

# メッシュも同時表示
Draw(mesh)
```

**ベストプラクティス**:
- フィールド分布のクイック確認
- メッシュ品質チェック
- パラメータ最適化時のリアルタイムフィードバック
- 教育・デモンストレーション

---

### Netgen GUI (形状・メッシュ確認 - 推奨)

**長所**:
- ✅ Netgen/NGSolveネイティブ（Tcl/Tk GUI）
- ✅ 形状（OCC）を**正確に**表示（近似なし）
- ✅ メッシュ品質可視化（アスペクト比、角度など）
- ✅ 軽量・高速起動（ブラウザ不要）
- ✅ STL/STEP/IGES読込
- ✅ 統合ワークフロー（形状確認→メッシュ生成→品質チェック）
- ✅ 通常のPythonスクリプト（.py）から使用可能

**短所**:
- ❌ フィールドデータ表示は限定的
- ❌ 論文品質レンダリング不可
- ⚠️ 古いGUI（Tcl/Tk）- ただし軽量で安定

**Radia-NGSolve統合での使用例**:
```python
from netgen.occ import OCCGeometry, Box, Pnt
from netgen.gui import StartGUI

# Radia磁石をOCC形状に変換
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
geo = OCCGeometry(occ_magnet)

# Netgen GUIで確認（ネイティブウィンドウが開く）
StartGUI()
geo.Draw()  # 形状を正確に表示

# メッシュ生成
mesh = geo.GenerateMesh(maxh=0.005)
mesh.Draw()  # メッシュ品質確認
```

**ベストプラクティス**:
- メッシュ生成前の形状確認（**webguiより正確**）
- メッシュ品質チェック（アスペクト比、角度）
- 境界条件ラベルの確認
- CADインポート後の形状検証
- Cubit → Netgen ワークフローの確認

**webguiとの使い分け**:
| 用途 | netgen.gui | ngsolve.webgui |
|------|-----------|---------------|
| 形状確認 | **✅ 推奨** | ⚠️ ブラウザオーバーヘッド |
| メッシュ品質 | **✅ 推奨** | ❌ 限定的 |
| フィールド可視化 | ❌ 不可 | **✅ 推奨** |
| 軽量・高速 | **✅ ネイティブGUI** | ⚠️ ブラウザ必要 |

**Windowsファイル関連付け** (TODO: 要検証):

.volファイルをダブルクリックでNetgen GUIで開く設定：

```powershell
# Netgen実行ファイルのパス（環境に応じて調整）
$netgenExe = "C:\Program Files\NGSolve\bin\netgen.exe"

# Windowsファイル関連付け設定
# 方法1: 右クリックメニューから
# 1. .volファイルを右クリック
# 2. "プログラムから開く" → "別のプログラムを選択"
# 3. "その他のアプリ" → "このPCで別のアプリを探す"
# 4. netgen.exe を選択
# 5. "常にこのアプリを使って.volファイルを開く" にチェック

# 方法2: コマンドライン（管理者権限必要）
assoc .vol=NetgenMeshFile
ftype NetgenMeshFile="C:\Program Files\NGSolve\bin\netgen.exe" "%1"
```

**注意**: Netgen GUIがコマンドライン引数から.volファイルを開けるかは要検証。
代替案として、Pythonスクリプトをラッパーとして使用することも可能：

```python
# netgen_viewer.py
import sys
from netgen.meshing import Mesh
from netgen.gui import StartGUI

if len(sys.argv) > 1:
    vol_file = sys.argv[1]
    mesh = Mesh()
    mesh.Load(vol_file)
    StartGUI()
    mesh.Draw()
```

このスクリプトを.volファイルに関連付ければ、ダブルクリックで開ける。

---

## Radia-NGSolve統合における推奨フロー

### フロー1: 開発・デバッグ（デフォルト）

```
Cubit → Netgen mesh → Radia solve → VTS export → **PyVista**
                                                     ↓
                                              迅速な確認・修正
```

### フロー2: 論文投稿

```
Cubit → Netgen mesh → Radia solve → VTS export → **ParaView**
                                                     ↓
                                         高解像度PNG/SVG出力
```

### フロー3: インタラクティブ探索

```
Cubit → Netgen mesh → Radia CF → NGSolve GridFunction → **webgui**
                                                            ↓
                                                  ブラウザで即座に確認
```

---

## 実装優先順位

### Phase 1: 基本可視化（現在完成）
- ✅ `rad.FldVTS()` - VTSエクスポート
- ✅ PyVistaによる基本プロット
- ✅ ParaViewでの手動可視化

### Phase 2: NGSolve統合強化（進行中）
- ✅ `RadiaField` CoefficientFunction
- ⏳ `radia_ngsolve` Python API改善
- ⏳ webgui連携スクリプト集

### Phase 3: 高度な可視化（TODO）
- ⏳ ParaView自動化スクリプト（pvpython）
- ⏳ PyVistaアニメーション生成
- ⏳ Radia解析オブジェクト → OCC形状変換

### Phase 4: 統合ビューワー（将来）
- ⏳ カスタムPyVistaインターフェース
- ⏳ Jupyter Widgetによるパラメータ制御
- ⏳ WebアプリケーションUI（Dash/Streamlit）

---

## ユースケース別推奨

| ユースケース | 推奨ビューワー | 理由 |
|-------------|--------------|------|
| **開発中の確認** | PyVista | 迅速、スクリプト統合 |
| **パラメータスタディ** | PyVista | バッチ処理、自動化 |
| **メッシュ品質チェック** | Netgen GUI | 軽量、専用機能 |
| **フィールド分布探索** | NGSolve webgui | インタラクティブ |
| **論文投稿図** | ParaView | 高品質、ベクトル出力 |
| **プレゼン資料** | ParaView | 高解像度、美しい |
| **アニメーション動画** | ParaView | キーフレーム機能 |
| **Jupyter解析** | PyVista + webgui | ノートブック統合 |
| **CI/CD回帰テスト** | PyVista (headless) | 自動化可能 |

---

## まとめ

### デフォルト推奨: **PyVista**
- 開発効率最優先
- 80%のユースケースをカバー
- Jupyter統合で解析フローが自然

### 論文・プレゼン: **ParaView**
- 品質が最重要な場合のみ
- 手動調整の手間を許容できる時

### NGSolve特化: **webgui**
- NGSolve GridFunctionを使う場合
- ブラウザベースが好ましい場合

### メッシュ確認: **Netgen GUI**
- 形状・メッシュの迅速チェック
- フィールド可視化は不要

---

**結論**: Radia-NGSolveフレームワークでは、**PyVistaをデフォルト**とし、必要に応じてParaView（論文）、webgui（インタラクティブ）、Netgen GUI（メッシュ）を使い分ける戦略が最適です。
