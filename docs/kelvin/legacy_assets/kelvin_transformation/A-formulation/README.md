# NGSolve Periodic Boundary Condition Issue Report

## 概要

NGSolve バージョン 6.2.2406 以降で、OCC ジオメトリの `Identify` メソッドを使用した周期境界条件が正しく適用されない後退バグが発生しています。

## 問題の詳細

### 症状
- `external_domain.faces[0].Identify(internal_domain.faces[0], "id")` を実行してもメッシュ生成時に識別情報が失われる
- `mesh.ngmesh.GetIdentifications()` の結果が空（0個）になる
- 周期境界条件が適用されない

### 影響を受けるコード
```python
from ngsolve import *
from netgen.occ import *

iron = Box((0,0,0),(1,1,1))
iron.name = "iron"
internal_domain = Sphere(Pnt(0,0,0.0), r=2)*Box((0,0,0), (5,5,5))
internal_domain.name = "air"
external_domain = Sphere(Pnt(3,0,0.0), r=2)*Box((3,0,0), (3+5,5,5))
external_domain.name = "air"

# 周期境界条件の設定
external_domain.faces[0].Identify(internal_domain.faces[0], "id")
domain = Glue([iron, internal_domain, external_domain])

geo = OCCGeometry(domain, dim=3)
mesh = Mesh(geo.GenerateMesh(mp)).Curve(2)

# 識別情報の確認
identifications = mesh.ngmesh.GetIdentifications()
print(f"Number of identifications: {len(identifications)}")  # 6.2.2406以降: 0
```

## バージョン別動作状況

| バージョン | 識別情報数 | 状態 | 備考 |
|----------|----------|------|------|
| 6.2.2404 | 20個 | ✓ 動作 | 正常動作確認済み |
| 6.2.2405 | 20個 | ✓ 動作 | 正常動作確認済み |
| 6.2.2406 | 0個 | ✗ 不具合 | `Identify`が機能しない（バグ発生） |
| 6.2.2501 | 0個 | ✗ 不具合 | |
| 6.2.2502 | 0個 | ✗ 不具合 | |
| 6.2.2503 | 0個 | ✗ 不具合 | |
| 6.2.2504 | 0個 | ✗ 不具合 | |
| 6.2.2505 | 0個 | ✗ 不具合 | |
| 6.2.2506 | 0個 | ✗ 不具合 | |
| 6.2.2601 | 20個 | ✓ 動作 | 推奨バージョン（バグ修正済み） |

## バージョン別の動作比較

### ✓ 正常動作: NGSolve 6.2.2404 および 6.2.2405

**出力（6.2.2404）:**
```
Number of identifications: 20
[(3.0, 2.0, 0.0), (0.0, 2.0, 0.0)]
[(3.0, 0.0, 2.0), (1.2246467991473532e-16, -2.999519565323715e-32, 2.0)]
[(5.0, -2.4492935982947064e-16, 8.057930683288193e-14), (2.0, -2.4492935982947064e-16, 8.057930683288193e-14)]
[(3.0, 1.868120247864882, 0.7142315727530197), (0.0, 1.868120247864882, 0.7142315727530197)]
...（全20ペア）...
```

![NGSolve 6.2.2404での動作](ngsolve_6.2.2404_working.png)
*図1: NGSolve 6.2.2404での周期境界条件の可視化*

![NGSolve 6.2.2405での動作](2025_09_07_3D_Kelvin変換の練習.png)
*図2: NGSolve 6.2.2405での周期境界条件の可視化（赤線が識別された点のペアを示す）*

上図では、2つの球形領域の境界面（outer）が正しく識別され、周期境界条件が適用されていることが確認できます。左側と右側の領域の対応する点が赤い線で結ばれており、これらの点で周期境界条件が課されています。

### ✗ 不具合: NGSolve 6.2.2406 以降

**出力（6.2.2406）:**
```
Number of identifications: 0

✗ No identifications found
```

![NGSolve 6.2.2406での不具合](ngsolve_6.2.2406_broken.png)
*図3: NGSolve 6.2.2406での不具合 - 識別情報が0個で赤い線が表示されない*

識別情報が0個となり、周期境界条件が全く適用されません。メッシュは生成されますが、2つの領域が独立したままです。図3では、図1・図2と異なり、赤い識別線が一切表示されていないことが確認できます。

## 回避策

### 推奨: バージョン 6.2.2601 以降にアップグレード

```bash
pip install --force-reinstall ngsolve>=6.2.2601
```

バージョン 6.2.2601 で周期境界条件のバグが修正されています。インストール後、Jupyter カーネルを再起動してください。

### 代替案

1. **CSG ジオメトリを使用する**
   - OCC の代わりに CSG ジオメトリと `PeriodicSurfaces` メソッドを使用
   - より確実に動作するが、ジオメトリの柔軟性が低い

2. **単純な形状を使用する**
   - 球との交差などの複雑なブーリアン演算を避ける
   - 単純な Box ジオメトリでは動作する可能性がある

## 原因の推測

バージョン 6.2.2406 で `Glue` メソッドまたは `Identify` メソッドの内部実装が変更され、識別情報がメッシュ生成プロセスで失われるようになったと考えられます。

## 今後の対応

このバグは NGSolve の開発チームに報告することを推奨します：
- リポジトリ: https://github.com/NGSolve/ngsolve
- フォーラム: https://ngsolve.org/forum

## 可視化結果の保存方法

### 方法1: VTK ファイルにエクスポート（推奨）

```python
from ngsolve import VTKOutput

# メッシュと解析結果をVTKファイルとして保存
vtk = VTKOutput(ma=mesh, coefs=[gfu], names=["solution"], filename="output", subdivision=2)
vtk.Do()

# ParaView で開いて高品質な画像を保存可能
# ParaView ダウンロード: https://www.paraview.org/download/
```

### 方法2: OS のスクリーンショット機能を使用（最も簡単）

Jupyter Notebook で `Draw` を実行した後：
1. **Windows**: `Win + Shift + S` でスクリーンショットツールを起動
2. 表示された WebGuiWidget の領域を選択
3. ペイント等に貼り付けて保存

**注意**: WebGL キャンバスは JavaScript からの自動スクリーンショット取得が制限されているため、手動でのキャプチャが最も確実です。

## 参考資料

- [NGSolve Periodicity Documentation](https://docu.ngsolve.org/latest/how_to/periodic.html)
- [NGSolve Forum: Saving netgen output](https://ngsolve.org/forum/ngspy-forum/780-saving-netgen-output)
- [NGSolve GitHub Repository](https://github.com/NGSolve/ngsolve)

## 検証環境

- OS: Windows
- Python: 3.12
- NGSolve: 6.2.2601以降（推奨）、6.2.2405（動作確認済み）
- 検証日: 2025年11月24日（6.2.2405）、2026年3月（6.2.2601）

## ライセンス

このドキュメントは情報共有を目的としています。
