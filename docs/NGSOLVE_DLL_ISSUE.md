# NGSolve Integration - DLL依存関係の問題

## 現状

`rad_ngsolve.pyd`は正常にビルドされましたが、インポート時に以下のエラーが発生します：

```
ImportError: DLL load failed while importing rad_ngsolve: 指定されたモジュールが見つかりません。
```

## 依存関係の確認結果

### rad_ngsolve.pydの直接依存

```
- python312.dll     [OK] - Pythonランタイム
- libngsolve.dll    [コピー済み] - NGSolveコアライブラリ
- KERNEL32.dll      [OK] - Windows API
- VCOMP140.DLL      [コピー済み] - OpenMPランタイム
```

### libngsolve.dllの依存（推定）

```
- ngcore.dll        [コピー済み]
- nglib.dll         [コピー済み]
- mkl_rt.lib        [不明] - Intel MKL
- その他のMKL DLL  [不明]
```

## 試行した解決策

### ✗ 失敗した方法

1. **DLLのコピー**: `netgen`フォルダから必要なDLLをコピー
   - 結果: 依然として同じエラー

2. **PATH環境変数**: netgenディレクトリをPATHに追加
   - 結果: 依然として同じエラー

3. **site-packagesへのコピー**: rad_ngsolve.pydをsite-packagesに直接コピー
   - 結果: 依然として同じエラー

## 根本原因の推定

NGSolveがIntel MKL（Math Kernel Library）に依存しており、これらのDLLが見つからない可能性が高いです。

NGSolveConfig.cmakeの内容：
```cmake
set(NGSOLVE_MKL_LIBRARIES "C:/gitlabci/tools/builds/3zsqG5ns/0/ngsolve/venv_ngs/Library/lib/mkl_rt.lib")
```

このパスは**ビルド環境専用**のパスで、実際のインストール環境には存在しません。

## 推奨される解決策

### オプション1: Condaを使用（最も確実）

Conda環境を使用すると、すべての依存関係が自動的に解決されます：

```bash
# NGSolve環境を作成
conda create -n ngsolve_env python=3.12
conda activate ngsolve_env
conda install -c ngsolve ngsolve

# Radiaをビルド（この環境内で）
cd S:\radia\01_GitHub
.\build_ngsolve.ps1

# テスト
python -c "import rad_ngsolve; print('OK')"
```

### オプション2: スタンドアロンビルド（MKL不要版）

NGSolveへのリンクを避け、ヘッダーのみを使用する実装に変更：

1. `rad_ngsolve.cpp`を修正してNGSolveライブラリへのリンクを削除
2. 必要な定義のみをヘッダーから取得
3. 実行時にngsolveモジュールから必要な機能をインポート

### オプション3: MKLをインストール

Intel oneMKLをシステムにインストール：

```powershell
# Intel oneMKL (無料)をダウンロードしてインストール
# https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html

# インストール後、PATH環境変数に追加
$env:PATH += ";C:\Program Files (x86)\Intel\oneAPI\mkl\latest\redist\intel64"
```

### オプション4: 動的ロード方式

Pythonレベルで実装を変更：

```python
# rad_ngsolve_wrapper.py
import ctypes
import os
import radia as rad

# DLLパスを動的に設定
ngsolve_path = r"C:\Program Files\Python312\Lib\site-packages\netgen"
os.add_dll_directory(ngsolve_path)

# rad_ngsolve.pydをロード
import rad_ngsolve as _rad_ngsolve

# ラッパー関数
def RadBfield(radia_obj):
	return _rad_ngsolve.RadBfield(radia_obj)
```

## 暫定的な回避策

現時点で最も現実的な解決策は、NGSolve統合を**純粋なPython実装**に変更することです：

```python
# rad_ngsolve_py.py - Pure Python implementation
import radia as rad
from ngsolve import CoefficientFunction

class RadiaBFieldPy(CoefficientFunction):
	def __init__(self, radia_obj):
	    super().__init__(3)
	    self.radia_obj = radia_obj

	def Evaluate(self, mip, result):
	    pnt = mip.GetPoint()
	    B = rad.Fld(self.radia_obj, 'b', [pnt[0], pnt[1], pnt[2]])
	    result[0] = B[0]
	    result[1] = B[1]
	    result[2] = B[2]

# 使用方法
magnet = rad.ObjRecMag([0,0,0], [20,20,30], [0,0,1000])
B = RadiaBFieldPy(magnet)
```

この方法の利点：
- ✓ DLL依存関係の問題なし
- ✓ インストール不要
- ✓ デバッグが容易

欠点：
- ✗ C++実装より遅い可能性
- ✗ NGSolveのCoefficientFunction継承が正しく動作しない可能性

## 次のステップ

1. **Conda環境でのテスト**を推奨（最も確実）
2. 失敗した場合は**純粋なPython実装**に切り替え
3. または**MKLをシステムにインストール**

## 参考情報

- NGSolve公式ドキュメント: https://ngsolve.org/
- Intel oneMKL: https://www.intel.com/content/www/us/en/developer/tools/oneapi/onemkl.html
- Windows DLL検索順序: https://learn.microsoft.com/en-us/windows/win32/dlls/dynamic-link-library-search-order
