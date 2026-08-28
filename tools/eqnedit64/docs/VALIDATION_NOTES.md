# Eqnedit64 検証記録

この文書の2026-08-25/26記録は、公開リポジトリ移行前の内部配布ゲートで得た
履歴である。現在の公開リリース手順は `build\accept_release.ps1` と
[`CANONICAL_OPERATION.md`](CANONICAL_OPERATION.md) を正とする。

## 2026-08-29 フォントホスト回帰の再現

- 対象: 公開済み`O:\Eqnedit64.exe` 3.0.1、SHA-256
  `E7CC8811312AB6C8E42AD5038936E50D4A64B72692224C1602C48095F20B94D0`
- 旧`--ui-interaction-test`は終了コード0だったが、同時刻に同一セッションの
  `fontdrvhost.exe`がPID 20808から7344へ入れ替わった。
- Application Event ID 1000/1001は`fontdrvhost.exe`、例外`0xc0000005`、
  offset`0x366a2`を2回記録していた。したがって旧CIのPASSは誤判定だった。
- `EQNEDIT64_NO_FONT_REG=1`では新規フォントホストクラッシュ0件、ファイル経由の
  `AddFontResourceExW(FR_PRIVATE | FR_NOT_ENUM)`単独プローブも終了コード0、PID維持、
  新規クラッシュ0件だった。
- 修正版の単一`--ui-interaction-test`は、分類タブの本文相当画素高を含め終了コード0、
  PID 7344を維持し、新規クラッシュ0件。32回のライフサイクル試験は対話中LABでは
  実行せず、GitHub-hosted Windows CIの必須ゲートとする。

## 2026-08-25 正本移行

- 対象コミット: `bdde3d450e6831931bd42df9446963b25cce7706`
- 実行: `pwsh -NoProfile -File build\accept_release.ps1 -Deploy`
- 結果: 全ゲート合格、登録配布先更新完了
- 配布物: `Eqnedit64.exe` 1ファイル、1,587,496 bytes
- SHA-256: `A9D0A05ED7C0A9B7DD946780252990092958B8E26C8934A0ADC144AC768669A6`
- Authenticode: `Valid`、`CN=ksugahar`
- build / dist / 登録配布先: サイズ、SHA-256、署名が一致
- 通常GUI耐久: 24 seed×3,000操作に加え、100 seed×5,000操作、合計500,000操作
- ASan: 自己試験、レイアウト、96/120/144/192 dpi相当画素試験、GUI操作、
  24 seed×3,000操作が合格
- PowerPoint: 左寄せ、編集可能なインラインOffice Math、24 pt
- IrfanView: 94×78の非空・不透明・黒文字画像
- Google Slides: 294×243 px、300 dpi、70.56×58.32 pt（24 pt基準）
- `--texclip`: Unicode / LaTeXから294×243 pxのPNG / DIBV5生成、全画素α=255
- 再描画実測（キャッシュなし / あり）:
  - `E=mc^2`: 0.158 / 0.090 ms
  - 代表複合式: 0.403 / 0.280 ms
  - 20分数: 0.819 / 0.612 ms
- Eqnedit32参照原本: SHA-256
  `3C4A68070F3D7F14E488AE4F7EDE8E7ADD0F8029995DC800833126CA062A2C6C`、
  Microsoft署名`Valid`、変更・削除なし
- 残留`Eqnedit64` / `Eqnedit64_asan`プロセス: なし
- 機械可読な実行報告:
  `C:\temp\Eqnedit64-final-acceptance.json`

この記録を追加する文書コミットは、上記の署名済み配布バイナリを変更しない。
次に製品コードまたは配布物を更新するときは、同じ最終ゲートを再実行する。

## 2026-08-26 任意サイズ行列

- 対象コミット: `7f6223fb66e84192ba2c2a10326ca1bffc7d5765`
- 実行: `pwsh -NoProfile -File build\accept_release.ps1 -Deploy`
- 結果: 全ゲート合格、登録配布先更新完了（817.939秒）
- 配布物: `Eqnedit64.exe` 1ファイル、1,597,736 bytes
- SHA-256: `A967B735C429DF2FBB18760CADE6FD6379D3FCEC59912E1D4406E6C92F288B6F`
- Authenticode: `Valid`、`CN=ksugahar`
- 行列モデル: 1×1～99×99、7×9直接生成、長方形の行／列追加・削除、
  既存セル保持、上下移動、名前付きUndo/Redo、最小／最大境界が合格
- TeX往復: 1×3と3×1の空の端セルを固定点で保持し、
  218記号／テンプレートのLaTeX実コンパイルが合格
- パレット: 18パレット、228セル、重複・未定義コマンド0件、
  実`WM_COMMAND`で全セルが状態変更
- メニュ監査: 253項目、241操作成功、12モーダル除外、操作可能項目のno-op 0件
- 構造ファズ: 100 seed×250操作、合計25,000操作でTeX固定点・有限寸法が合格
- 通常GUI耐久: 24 seed×3,000操作に加え、100 seed×5,000操作、合計500,000操作
- ASan: 自己試験、レイアウト、96/120/144/192 dpi相当画素試験、GUI操作、
  24 seed×3,000操作が合格
- PowerPoint: 左寄せ、編集可能なインラインOffice Math、24 pt
- IrfanView: 94×78の非空・不透明・黒文字画像
- Google Slides: 294×243 px、300 dpi、70.56×58.32 pt（24 pt基準）
- `--texclip`: Unicode / LaTeXから294×243 pxのPNG / DIBV5生成、全画素α=255
- 再描画実測（キャッシュなし / あり）:
  - `E=mc^2`: 0.148 / 0.085 ms
  - 代表複合式: 0.385 / 0.274 ms
  - 20分数: 0.824 / 0.633 ms
- ビルド前の実行中リポジトリ版終了: 解決済み絶対パスだけを対象にする試験が合格
- Eqnedit32参照原本: SHA-256
  `3C4A68070F3D7F14E488AE4F7EDE8E7ADD0F8029995DC800833126CA062A2C6C`、
  Microsoft署名`Valid`、変更・削除なし
- 機械可読な実行報告: `C:\temp\Eqnedit64-final-acceptance.json`

この記録追加は署名済み配布バイナリを変更しない。
