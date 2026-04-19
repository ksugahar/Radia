# Troubleshooting playbook

## 問題切り分けの基本方針

層別に切り分ける (CLAUDE.md の 4 層):
- Layer 1 (C++ Cubit plugin): ccm/ccl の挙動
- Layer 2 (Cubit Python): panel 起動、register_toolbar
- Layer 3 (PySide6 panel): IHWindow の GUI + subprocess Popen
- Layer 4 (Computation): calc_*.py の physics / solve

## Panel が起動しない

### 症状: Cubit で Solve → Radia-NGSolve を押しても反応なし

```powershell
# 1. Cubit の journal log 確認
# Tools → Show Journal Window
# エラーメッセージ (ModuleNotFoundError, ImportError 等) があるか?
```

原因候補:
1. register_toolbar.py が Cubit の自動 start scripts に入ってない  
   → `cubit-plugin-install` 再実行
2. Python 3.12 が見つからない  
   → `where python` で確認、install_panels.py 側の path resolution を確認
3. 古い panel が残っている  
   → `C:\Users\<user>\.cubit\` を削除、`cubit-plugin-install` 再実行

### 症状: Window 出たが method combo が空

```python
# Python console で直接:
import sys
sys.path.insert(0, r'C:\Program Files\Python312\Lib\site-packages\radia')
from radia_ih import IHPanel
p = IHPanel()
print([p._method_combo.itemText(i) for i in range(p._method_combo.count())])
```

期待値: `['PEEC+BEM (1-way)', 'FEM (coil meshed + SIBC + Kelvin)']`  
空なら radia_ih.py 読み込みが古い。site-packages の radia_ih.py を確認:
```bash
ls -la 'C:\Program Files\Python312\Lib\site-packages\radia\radia_ih.py'
# mtime / size を LAB 側と比較
```

## Calc subprocess がクラッシュする

### 症状: Run 押すと stderr に traceback

Panel は subprocess stderr を capture するので、Run ボタン横のログに出る。出てなければ:
- 手動実行で再現:
  ```bash
  python "C:\Program Files\Python312\Lib\site-packages\radia\panels\calc_peec_bem.py" \
    --peec-step "C:\path\to\coil.step" \
    --vol "C:\path\to\model.vol" \
    --frequency 7000 --current 1.0 --coil-sigma 5.8e7 \
    --sigma 5.8e7 --half-thickness 0.0125 --mu-r 1.0 \
    --peec-nwinc 3 --peec-nhinc 3 --wp-label sibc
  ```
  
- よくある原因:
  - **DLL load failed peec_matrices**: `import radia` が calc_*.py 冒頭に無い。MKL DLL path 未設定。(現版は修正済みだが過去の calc に注意)
  - **wp_label 'sibc' not found**: .vol に sideset `sibc` が無い。.jou を見直し
  - **calc_common ImportError**: `sys.path` に panel folder が入ってない
  - **pardiso out of memory**: ndof 大 + メモリ不足。coarser mesh か別 solver (`--solver bddc`) へ

### 症状: 結果数字が怪しい (P_wp が 2D ref と桁外れ)

Golden test を回して再現確認:
```bash
pytest tests/panels/test_peec_bem_golden.py -v -s
pytest tests/panels/test_fem_coilmesh_golden.py -v -s
```

Golden が通れば **panel の受け渡しに問題**、通らなければ **calc 側の physics bug**。

Physics bug の切り分け:
- **P_wp が 0 に近い**: Normal orientation 反転 (BND 抽出バグ)。`_extract_bnd_only` の `n_flipped` 数値を stderr でチェック
- **P_wp が無限大 / NaN**: source/sink 周辺の singularity。Gap face が `coil_surface` sideset に入ってる?
- **L が 10x 小さい**: Robin SIBC が K_imposed を相殺してる (Neumann-K 案の再発、2026-04-19 retired)
- **L がマイナス**: sign convention の反転、Ampere loop 向き間違い

## .vol 読み込みエラー

### 症状: `Mesh("my.vol")` で ValueError

- **version mismatch**: Netgen fork (post29) vs 公式 6.2.2603 で .vol format 差異  
  → `radia_export netgen "..." order 1` で再 export
- **無名 boundary がある**: 古い Cubit export で Surface_N 自動命名  
  → sideset で明示的に名前付ける
- **NaN 座標**: Cubit の subtract/imprint が退化面作る  
  → .jou の subtract 順、keep_tool 指定を見直し

### 症状: `Boundary label 'sibc' not found`

```python
# 実際のラベルを確認
import radia
from ngsolve import Mesh
m = Mesh("my.vol")
print("Materials:", sorted(set(m.GetMaterials())))
print("Boundaries:", sorted(set(m.GetBoundaries())))
```

対応:
- `sibc` 別名になってる (例: `wp_surface`, `Surface_3`) → panel で --wp-label を変更 or .jou で sideset name 修正
- sibc sideset 自体が無い → .jou に `sideset N add surface in volume {wp_vid}; sideset N name "sibc"` 追加

## MCP server が更新されない

### 症状: Claude Code に聞いても古い情報が返る

1. **MCP server の reload**:
   - Claude Code を restart、または
   - `/mcp` で reload コマンド
   
2. **LAB S: と 100号機 W: が同一 share か確認**:
   ```bash
   ssh 100 'powershell.exe -Command "(Get-Item W:\00_CAE\Radia\01_GitHub\src\radia\mcp_server\ih\ih_knowledge.py).LastWriteTime"'
   stat -c '%y' s:/Radia/01_GitHub/src/radia/mcp_server/ih/ih_knowledge.py
   ```
   Mtime 同じなら share、違うなら 100 で git pull 必要

3. **site-packages の radia/mcp_server/ が stale**:
   MCP は repo path 直接参照なのでここは無関係のはずだが、念のため確認

## Cubit plugin .ccm/.ccl が認識されない

### 症状: Cubit 起動時に "Plugin not loaded: radia_cubit.ccm"

- Path 確認: `C:\Program Files\Coreform Cubit 2025.3\bin\plugins\radia_cubit.ccm`
- サイズ確認: 689,664 bytes (2026-04-19 build)
- Cubit version: 2025.3 限定。古い version では ABI 違う
- **cubit.cmd("radia_export ...")** で試すと invalid command エラー出る

対応:
- `cubit-plugin-install` で再配布
- または SMB で直接コピー:
  ```bash
  cp src/radia/radia_cubit.ccm '//192.168.11.100/c$/Program Files/Coreform Cubit 2025.3/bin/plugins/'
  ```

## Panel の save state が壊れる

### 症状: 起動するたびに method combo が空 or 違う

- Settings file: `C:\Users\<user>\AppData\Roaming\Radia\ih.ini` (or similar)
- Text-based save なので古い method name (例: "PEEC+FEM") が残っていると `findText` で −1 → 空
- 対応: settings file を削除して default にリセット

## Deploy 絡みの罠

### 症状: 100号機で panel 動作が古い

- **SMB share mount の状態確認**:
  ```bash
  ssh 100 'Test-Path "C:\Program Files\Python312\Lib\site-packages\radia\panels\calc_peec_bem.py"'
  ```
- **site-packages と LAB S: の両方を見る**: 100号機 site-packages は LAB repo の snapshot
- **Cubit がロックしてる**: Cubit 開きっぱなしだと .ccm/.ccl 書き換わらない
  ```bash
  ssh 100 'Stop-Process -Name coreform_cubit -Force'
  ```
  再 deploy

### 症状: MCP が古い内容を返す (LAB 側)

- LAB の MCP server はどこから読んでる?
  Claude Code の設定 `.mcp.json` を確認
- Repo 直接参照 or pip installed?
- pip installed (`site-packages/radia/mcp_server/...`) なら editable install で最新反映:
  ```bash
  pip install -e . --no-deps
  ```

## よくある "これが原因でした" 事例

| 症状 | 真の原因 | Fix |
|---|---|---|
| BEM P_wp が 20x 小さい | sibc FD が複数あるうち 1 つしか使ってない | `target_indices = {...}` で全 FD 拾う |
| BEM P_wp が 4x 大きい | BND 三角形の winding が FD ごとに逆向き | Centroid-based orientation check |
| FEM L が 10x 小さい | Robin SIBC on coil が K_imposed 相殺 | Coil 側 Robin は外す (Neumann-K) |
| FEM P_coil が 75% 高い | Coil mesh h > delta | `volume size` を delta/3 に |
| `Gimbal lock detected` warning | CoilBuilder の Euler 角 decomposition | Z 軸対称 torus で正常、無視 OK |
| `コイル h_max > delta` WARN | 内部 tet が大きい | mesh 改善 or P_coil の精度諦める |

## Escalate

切り分け不能 or 重大な physics 疑義:
- Slack #radia-dev
- 菅原 (ksugahar@ele.kindai.ac.jp)
- 解析結果 JSON + stderr + golden 結果 添付
