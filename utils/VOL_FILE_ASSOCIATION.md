# .vol File Association for Netgen GUI (Windows)

## 概要

.volファイルをダブルクリックでNetgen GUIで開けるようにするWindowsファイル関連付け設定。

---

## 方法1: バッチファイルで自動設定（推奨）

### 手順

1. **管理者権限でコマンドプロンプトを開く**:
   - スタートメニュー → "cmd" と入力
   - "コマンドプロンプト" を右クリック → "管理者として実行"

2. **セットアップスクリプトを実行**:
   ```cmd
   cd S:\Radia\01_GitHub\utils
   setup_vol_file_association.bat
   ```

3. **完了**:
   - レジストリが設定されます
   - Windows Explorerを再起動（またはPC再起動）

4. **動作確認**:
   - 任意の.volファイルをダブルクリック
   - Netgen GUIが自動的に開くはず

---

## 方法2: 右クリックメニューから手動設定

### 手順

1. **.volファイルを右クリック**:
   - 任意の.volファイル（例: `cube.vol`）を右クリック

2. **"プログラムから開く" → "別のプログラムを選択"**

3. **"その他のアプリ"をクリック**

4. **"このPCで別のアプリを探す"をクリック**

5. **Pythonスクリプトを選択**:
   - `S:\Radia\01_GitHub\utils\netgen_vol_viewer.py` を選択
   - または、Pythonインタープリタ（`python.exe`）を選択し、引数を手動設定

6. **"常にこのアプリを使って.volファイルを開く" にチェック**

7. **OKをクリック**

---

## 方法3: レジストリ直接編集（上級者向け）

### レジストリ設定

```reg
Windows Registry Editor Version 5.00

; .vol拡張子をNetgenMeshFileに関連付け
[HKEY_CLASSES_ROOT\.vol]
@="NetgenMeshFile"

; NetgenMeshFileの開くコマンドを定義
[HKEY_CLASSES_ROOT\NetgenMeshFile]
@="Netgen Mesh File"

[HKEY_CLASSES_ROOT\NetgenMeshFile\shell]

[HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open]

[HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open\command]
@="python.exe \"S:\\Radia\\01_GitHub\\utils\\netgen_vol_viewer.py\" \"%1\""
```

### 適用手順

1. 上記内容を `netgen_vol_association.reg` として保存
2. ファイルをダブルクリック
3. "レジストリに追加しますか？" → "はい"

---

## トラブルシューティング

### 問題1: ダブルクリックしても何も起こらない

**原因**: レジストリが正しく設定されていない

**解決策**:
1. レジストリエディタを開く（`regedit`）
2. `HKEY_CLASSES_ROOT\.vol` を確認
3. 値が `NetgenMeshFile` になっているか確認
4. `HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open\command` を確認
5. コマンドが正しいか確認

### 問題2: Pythonが見つからないエラー

**原因**: PythonがPATHに設定されていない

**解決策**:
1. Pythonのフルパスを使用:
   ```
   C:\Users\{username}\AppData\Local\Programs\Python\Python312\python.exe
   ```
2. または、`pythonw.exe` を使用（コンソールウィンドウが出ない）

### 問題3: Netgen GUIが開かない

**原因**: NGSolveがインストールされていない

**解決策**:
```cmd
pip install ngsolve
```

### 問題4: 複数のPython環境がある

**原因**: 複数のPython環境で、NGSolveがインストールされていない環境が使われている

**解決策**:
1. NGSolveがインストールされているPython環境を特定:
   ```cmd
   where python
   python -c "import ngsolve; print(ngsolve.__file__)"
   ```
2. そのPythonのフルパスをレジストリに設定

---

## 代替案: Pythonwを使用（コンソールなし）

コンソールウィンドウを表示したくない場合：

```reg
[HKEY_CLASSES_ROOT\NetgenMeshFile\shell\open\command]
@="pythonw.exe \"S:\\Radia\\01_GitHub\\utils\\netgen_vol_viewer.py\" \"%1\""
```

**注意**: `pythonw.exe` はエラーメッセージが表示されないため、デバッグが困難。

---

## 動作確認

### テスト方法

1. NGSolveサンプルメッシュを使用:
   ```
   S:\NGsolve\01_GitHub\install_ksugahar\share\ngsolve\cube.vol
   ```

2. ダブルクリック

3. Netgen GUIが開き、cubeメッシュが表示されるはず

### 期待される動作

- Netgen GUIウィンドウが開く
- メッシュが3D表示される
- マウスドラッグで回転可能
- マウスホイールでズーム可能

---

## 関連ファイル

| ファイル | 説明 |
|---------|------|
| `netgen_vol_viewer.py` | .volファイルを開くPythonスクリプト |
| `setup_vol_file_association.bat` | 自動設定バッチファイル |
| `VOL_FILE_ASSOCIATION.md` | このドキュメント |

---

## 参考: 他のファイル形式

同様の方法で、他のメッシュ形式も関連付け可能：

| 拡張子 | フォーマット | ビューワー |
|-------|------------|----------|
| `.vol` | Netgen mesh | Netgen GUI |
| `.msh` | Gmsh mesh | Gmsh, ParaView |
| `.vtu` | VTK unstructured | ParaView, PyVista |
| `.stl` | STL surface | Netgen GUI, ParaView |

---

**最終更新**: 2026-02-12
**対象OS**: Windows 10/11
**要件**: Python 3.x, NGSolve
