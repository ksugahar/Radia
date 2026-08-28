# Eqnedit64 正本運用

## 正本の定義

2026-08-25以降、研究室の数式入力ツールは **Eqnedit64を正** とする。

- 文書資産の正本はUTF-8の `.tex`。新規文書は `equation`、複数行はその内側の
  `aligned` で保存する。
- アプリケーションの正本は、最終受入ゲートを通った
  `dist\Eqnedit64.exe`。公開物はGitHub Releaseへ署名済みEXEとSHA-256を登録する。
- GUI動作の正本は `docs\GUI_SPEC.md`、ショートカットは
  `docs\SHORTCUTS.md`、採用済みの操作判断は
  `docs\USABILITY_PREFERENCES.jsonl` とする。
- インストーラー、Python、VC++再頒布ランタイム、レジストリ登録を運用要件に
  しない。署名済み `Eqnedit64.exe` 一つで動かす。

`.eqn`とMTEFは新しい研究室資産として作成・保存・受け渡ししない。既存資産の
変換が必要な場合だけ、明示的な移行作業として扱い、変換後の `.tex` を正本にする。

## Eqnedit32の位置づけ

Eqnedit32は日常編集、配布、ファイル関連付けの正本ではない。Microsoft由来の
バイナリ、MTEF変換コード、逆アセンブリ資料は公開Eqnedit64ソースに含めず、
研究室の調査済み資産として別管理する。Eqnedit64のビルド・試験・配布は旧版に
依存しない。

## リリースと配布

正本を更新するときは、画面・マウス・実キーボードを操作しない最終ゲートを実行する。
リポジトリ内の `build\Eqnedit64.exe` または `dist\Eqnedit64.exe` が起動中の場合、
ビルドはファイルロック解除のためその実行中プロセスを強制終了する。
別の配布先と無関係な同名プロセスは対象にしない。

```powershell
pwsh -NoProfile -File build\accept_release.ps1
```

この一つのコマンドが、通常回帰、PowerPoint COM、IrfanView CLI、Google Slides
形式、ASan、100 seed×5,000回の非表示Win32操作、96/120/144/192 dpi相当の
オフスクリーン視認性、最大式でも1再描画5 ms未満の性能ゲート、
`CN=ksugahar`署名、単体exeを確認する。途中で一つでも失敗した場合は公開しない。結果は既定で
`C:\temp\Eqnedit64-final-acceptance.json` に残る。

公開順は固定する。まずリリースコミットを`main`へpushしてmain CIを通し、その
`origin/main`と同じコミットからLABで署名済みEXEをビルドする。次に
`.agents\skills\release-eqnedit64\scripts\sync_to_o.ps1`で
`O:\Eqnedit64.exe`を更新し、隣の`O:\Eqnedit64.release.json`へ予定tag、
source SHA、EXE SHA-256、version、signerを記録する。ここまで成功してから最後に
`eqnedit64-v<version>`tagをpushする。

tag CIはO:のマニフェストに記録されたsource SHAがtag SHAと完全一致する場合だけ、
その署名済みEXEからPython 3.10--3.13用wheelを作り、PyPIへ公開した後に
GitHub Releaseへ`Eqnedit64.exe`と`SHA256SUMS.txt`を添付する。O:更新前にtagを
pushしてはならない。旧版へ戻す必要が生じた場合も、対象Gitコミットを再ビルドして
同じ順序と最終ゲートを通す。

## 正本移行後の品質運用

人間にしか判断できないIME候補の見た目や長時間入力の感触は、完成前の曖昧な
一回確認ではなく、正本運用中も継続して観察する。違和感の直後に `F12` を記録し、
重大・高優先度の問題が見つかった場合は回帰試験を追加して修正し、再び最終ゲートを
通す。現時点で既知の重大・高優先度の未解決項目はない。
