# Radia v4.7.0 リリースノート (2026-04-22)

3 パッケージ同時リリース:

| Package | Version | PyPI |
|---|---|---|
| `radia` | **4.7.0** | https://pypi.org/project/radia/4.7.0/ |
| `cubit-mesh-export` | **0.6.0** | https://pypi.org/project/cubit-mesh-export/0.6.0/ |
| `radia-mcp` | **0.32.0** | https://pypi.org/project/radia-mcp/0.32.0/ |

## 今回の目玉

### 1. PEEC-inductance が任意 STEP で動くように

Cubit で作ったコイル STEP をそのまま "PEEC inductance (coil only,
STEP)" パネルに食わせれば L / R が出ます。`.jou` の同梱は不要。

| コイル形状 | before | after |
|---|---|---|
| 1 周 円形トーラス (gapped) | 100 倍ずれ | ✅ L = 85 nH (解析式 -3.9%) |
| 1 周 矩形断面トーラス | 動作せず | ✅ L = 88 nH |
| 3 ターン loft (Kubota の 10 MB STEP) | walker ハング >5 分 | ✅ 12.9 秒で L = 431 nH (`.jou` 経路 426 nH と +1.1 %) |
| `.jou` 直接 (明示 centerline) | | ✅ 継続サポート、最速 |

サイドバーに `3turncoil.stp` + `3turnCoil.jou` 両方ある場合、自動的に
`.jou` 経路へ切替 (メッセージ `found sibling .jou, preferring it`)。

### 2. 日本語パス対応

Cubit plugin の全 6 エクスポーター (Netgen / GMSH / Nastran / VTK /
MEG / FEMEEM) で `radia_export netgen "C:/temp/日本語/coil.vol"` が
書き込めます。UTF-8 → wide API 経由、システムコードページに依存しません。

パネル経由 Python 処理 (calc_peec_inductance, calc_fem_coilmesh etc)
も日本語パスを透過的に扱います。

### 3. Cubit 起動が速くなる

- **初回 GUI 起動**: 30-60 秒 → **3 秒** (ライセンス事前更新)
- **VSCode 再起動後**: 6 秒再生成 → **0.01 秒 アタッチ** (Phase 1 daemon 保持)
- インポート済の STEP や作ったメッシュも VSCode 再起動で失われません

---

## 100号機 (lab 共有機) 使用ユーザーへの連絡事項

### 一度ログアウト + 再ログインをお願いします

背景: v4.7.0 で Windows ログオン時にライセンスを自動更新する scheduled
task (`\Coreform\CubitLicenseRefresh`) を登録しました。次回ログインで
自動発火するので、Cubit 起動が常に 3 秒以内になります。

**既ログイン中のユーザー** (Kubota / keiko / yano / 他 11 人): 一度
Windows からログアウトして再ログインしてください。以下のどちらか
速い方で代用できます:

1. **デスクトップのショートカット**:
   `Coreform Cubit (warm launch)` をダブルクリック
   → キャッシュが古ければ自動更新 + Cubit 起動

2. **即時更新 (Cubit はまだ起動しない)**:
   ```
   C:\ProgramData\CoreformCubit\cubit_refresh.cmd
   ```
   をダブルクリック (5秒以内で完了)

どちらでも自分のライセンスキャッシュが更新され、以降 3 日間は
Cubit 起動 3 秒です。

### VSCode を使っている方へ

VSCode 内の Claude Code MCP サーバーは起動時にロードしたコードを
キャッシュします。今回の v4.7.0 の Cubit daemon 高速化を反映するには
**一度 VSCode を終了 → 再起動** してください。

再起動後:
- MCP `open_in_cubit` / `cubit_show` / `cubit_exec` で Cubit を呼ぶ = **3 秒**
- 別の作業で VSCode を再度再起動した場合も Cubit は生存 → **0.01 秒で再接続**

### Kubota さんの `3turncoil.stp` ケース

- `W:\kubota\3turncoil.stp` を "PEEC inductance" パネルに Browse → 自動で
  sibling `3turnCoil.jou` を採用 → **L = 426.245 nH (約 6 秒)**
- `.jou` なしでも動きます (**L = 430.86 nH、約 13 秒** 別ディレクトリ試験済)

---

## 技術者向け変更点サマリー

- `src/radia/coil_from_cad.py` — revolution-sweep extractor、
  cross-section centroid chain、unit auto-detection
- `src/radia/panels/calc_peec_inductance.py` — sibling `.jou` auto-
  preference (PEEC format guard 付き)
- `src/cubit_plugin/utf8_path.hpp` — 全 6 エクスポーターで UTF-8 path 対応
- `packages/radia-mcp/src/radia_mcp/cubit/cubit_license_warmup.py` —
  `rlm_activate` キャッシュチェック + 必要時ログイン
- `packages/radia-mcp/src/radia_mcp/cubit/cubit_session.py` — Phase 1
  daemon persistence (per-user stable drop-dir + `pid.lock` discovery +
  `CREATE_NEW_PROCESS_GROUP` + `DEVNULL` pipes)
- `tests/panels/test_peec_inductance_golden.py` — 3 path golden locks
- `C:\ProgramData\CoreformCubit\` — warmup script + `.cmd` + Public
  Desktop shortcut + scheduled task (any-user logon trigger)

## アップグレード手順

**mdx や他マシンで PyPI 経由**:
```
pip install --upgrade radia cubit-mesh-export radia-mcp
cubit-plugin-install
```

**100号機 (lab 共有 LAB 経由)**:
既に `release_triple.py phase8` で deploy 済、次回各ユーザーの次回
logon か上記 warmup 起動で反映。

**LAB (開発機)**: editable install 済。

---

問題があれば Sugahara ksugahar@ele.kindai.ac.jp まで。
