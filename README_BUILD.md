# Radia Python モジュール ビルドガイド

## Pythonバージョンの互換性について

### 重要な制約

**Python .pydファイルは、Pythonのメジャー・マイナーバージョンごとに異なるバイナリが必要です。**

これは、Python C APIの仕様により、以下の理由があります：

1. **ABI互換性がない**: Python 3.8, 3.9, 3.10, 3.11, 3.12 は、それぞれ異なるバイナリインターフェース
2. **構造体サイズの違い**: 内部データ構造のサイズやレイアウトがバージョン間で異なる
3. **関数シグネチャの変更**: API関数の引数や戻り値が変更される場合がある

### 解決策

**各Pythonバージョン用に個別のpydファイルをビルドする**

## ビルド方法

### 方法1: 単一バージョン用ビルド

現在インストールされているPythonバージョン用にビルド：

```powershell
# 基本ビルド
.\Build.ps1

# Debugビルド
.\Build.ps1 -BuildType Debug

# クリーンビルド
.\Build.ps1 -Rebuild
```

**出力:**
- `build/lib/Release/radia.pyd` - 現在のPython用

### 方法2: 複数バージョン一括ビルド（推奨）

複数のPythonバージョン用に一括ビルド：

```powershell
# 見つかったすべてのPythonバージョン用にビルド
.\build_multi_python.ps1

# ビルド後に各バージョンでテストも実行
.\build_multi_python.ps1 -TestAll

# クリーンビルド
.\build_multi_python.ps1 -Clean
```

**出力例:**
```
output/
├── radia.cp38-win_amd64.pyd    # Python 3.8用
├── radia.cp39-win_amd64.pyd    # Python 3.9用
├── radia.cp310-win_amd64.pyd   # Python 3.10用
├── radia.cp311-win_amd64.pyd   # Python 3.11用
└── radia.cp312-win_amd64.pyd   # Python 3.12用
```

### ファイル命名規則

PEP 3149に従った命名：

```
radia.cp<version>-<platform>.pyd

例:
- radia.cp312-win_amd64.pyd  → Python 3.12 (64-bit Windows)
- radia.cp311-win_amd64.pyd  → Python 3.11 (64-bit Windows)
```

## 必要な環境

### 必須

- Windows 10/11 (64-bit)
- Visual Studio 2022 (Community以上)
- CMake (Visual Studio付属)
- Python 3.x (64-bit) - 1つ以上

### 対応Pythonバージョン

スクリプトはデフォルトで以下のバージョンに対応：

- Python 3.8
- Python 3.9
- Python 3.10
- Python 3.11
- Python 3.12

※ インストールされているバージョンのみビルドされます

## Python のインストール

### 推奨: 公式インストーラー

https://www.python.org/downloads/windows/

**重要: 必ず64-bit版をインストールしてください**

### インストール先の確認

デフォルトのインストール先：
```
C:\Python38\
C:\Python39\
C:\Python310\
C:\Python311\
C:\Program Files\Python312\
```

カスタムパスを使用する場合は、`build_multi_python.ps1`を編集：

```powershell
$PythonVersions = @(
	@{
	    Version = "3.11"
	    Path = "D:\MyPython\Python311"  # カスタムパス
	    MinorVersion = "3.11"
	}
)
```

## 使用方法

### 配布

各Pythonバージョン用の`.pyd`ファイルを配布してください。

### インストール

ユーザーは自分のPythonバージョンに合った`.pyd`ファイルを使用：

```python
# Python 3.12の場合
# radia.cp312-win_amd64.pyd を radia.pyd にリネームして
# site-packages/ または作業ディレクトリに配置

import radia as rad
print(rad.UtiVer())
```

### 自動選択

以下のようなラッパーで自動選択も可能：

```python
# radia_loader.py
import sys
import os
import importlib.util

def load_radia():
	version = f"cp{sys.version_info.major}{sys.version_info.minor}"
	pyd_name = f"radia.{version}-win_amd64.pyd"

	# 同じディレクトリから探す
	module_dir = os.path.dirname(__file__)
	pyd_path = os.path.join(module_dir, pyd_name)

	if not os.path.exists(pyd_path):
	    raise ImportError(f"Radia not available for Python {sys.version_info.major}.{sys.version_info.minor}")

	spec = importlib.util.spec_from_file_location("radia", pyd_path)
	radia = importlib.util.module_from_spec(spec)
	spec.loader.exec_module(radia)
	return radia

rad = load_radia()
```

## トラブルシューティング

### Python not found

```
[WARN] Python 3.X not found at: C:\PythonXX
```

**解決策:**
1. Pythonをインストール
2. または `build_multi_python.ps1` のパス設定を変更

### 32-bit Python detected

```
[WARN] Python 3.X is not 64-bit
```

**解決策:**
64-bit版Pythonを再インストール

### CMake not found

```
ERROR: CMake not found in Visual Studio 2022
```

**解決策:**
Visual Studio Installerで「C++によるデスクトップ開発」をインストール

## ビルド成果物

### ディレクトリ構造

```
04_Radia/
├── Build.ps1                    # 単一バージョンビルド
├── build_multi_python.ps1       # マルチバージョンビルド
├── CMakeLists.txt               # CMake設定
├── output/                      # マルチビルド出力
│   ├── radia.cp38-win_amd64.pyd
│   ├── radia.cp39-win_amd64.pyd
│   ├── radia.cp310-win_amd64.pyd
│   ├── radia.cp311-win_amd64.pyd
│   └── radia.cp312-win_amd64.pyd
├── build/                       # 単一ビルド出力
│   └── lib/Release/radia.pyd
├── build_py38/                  # Python 3.8ビルドディレクトリ
├── build_py39/                  # Python 3.9ビルドディレクトリ
└── ...
```

## テスト

### 単一バージョンテスト

```powershell
python test_simple.py
python test_radia.py
```

### 全バージョン自動テスト

```powershell
.\build_multi_python.ps1 -TestAll
```

これにより、ビルドされた各バージョンで自動的に動作確認が行われます。

## まとめ

**Pythonバージョンごとに異なるpydが必要**という制約がありますが、
`build_multi_python.ps1` を使用することで、複数バージョンを効率的にビルドできます。

配布時は、サポートする各Pythonバージョン用のpydファイルを用意してください。
