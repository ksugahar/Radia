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

`O:\Eqnedit64.exe`は、菅原が候補版を実IME・実PowerPointなどで確認するための
**正規ハンドテスト入口**でもある。候補をコミットし、署名済み単体EXEのローカル
バックグラウンド試験が通ったら、テストを依頼する前に次を実行する。候補版を
`C:\temp`やworktreeだけに置いたまま、O:上の旧版をテスト対象にしてはならない。

```powershell
pwsh -NoProfile -File `
  .agents\skills\release-eqnedit64\scripts\sync_handtest_to_o.ps1 `
  -SourceExe tools\eqnedit64\dist\Eqnedit64.exe `
  -SourceSha (git rev-parse HEAD)
```

ハンドテスト同期にはversion更新、main push、tagを要求しない。ただしsourceはcleanな
commit、EXEは同じbuild stamp、有効な`CN=ksugahar`署名でなければならない。同期時は
`O:\Eqnedit64.handtest.json`を作り、直前の正式マニフェストを
`O:\Eqnedit64.last-release.json`へ保存し、`O:\Eqnedit64.release.json`を外す。
これにより候補版を正式公開物と誤認したtag公開をfail-closedにする。

正本を更新するときは、画面・マウス・実キーボードを操作しない最終ゲートを
GitHub-hosted Windowsの隔離CIで実行する。
リポジトリ内の `build\Eqnedit64.exe` または `dist\Eqnedit64.exe` が起動中の場合、
ビルドはファイルロック解除のためその実行中プロセスを強制終了する。
別の配布先と無関係な同名プロセスは対象にしない。

`accept_release.ps1`は使い捨てVMまたは専用の隔離ユーザーセッションでCI相当を再現する
場合だけ使う。対話中LABでは実行せず、ガードを通すために
`EQNEDIT64_ISOLATED_TEST_SESSION=1`を設定してはならない。

```powershell
pwsh -NoProfile -File build\accept_release.ps1
```

この一つのコマンドが、通常回帰、PowerPoint COM、IrfanView CLI、Google Slides
形式、ASan、100 seed×5,000回の非表示Win32操作、96/120/144/192 dpi相当の
オフスクリーン視認性、最大式でも1再描画5 ms未満の性能ゲート、
`CN=ksugahar`署名、単体exeを確認する。途中で一つでも失敗した場合は公開しない。結果は既定で
`C:\temp\Eqnedit64-final-acceptance.json` に残る。
`fontdrvhost.exe`交代の事前対照窓にも同じApplication Errorがある、または交代にcrash
eventが伴わない場合は製品FAILではなく`INCONCLUSIVE`として停止する。そのrunnerを
greenになるまで回し直して合格扱いにせず、健全な新しい隔離sessionで候補を再検査する。

公開順は固定する。まずリリースコミットを`main`へpushして隔離main CIを通し、その
`origin/main`と同じコミットからLABではコンパイルと署名だけを行う。次に
`.agents\skills\release-eqnedit64\scripts\sync_to_o.ps1`で
`O:\Eqnedit64.exe`を更新し、隣の`O:\Eqnedit64.release.json`へ予定tag、
source SHA、EXE SHA-256、version、signerを記録する。ここまで成功してから最後に
`eqnedit64-v<version>`tagをpushする。

正式同期の成功時は`O:\Eqnedit64.handtest.json`を除去し、現在のEXEと一致する
`O:\Eqnedit64.release.json`を再作成する。同じトランザクションで、この2ファイルを
`Eqnedit64-<version>.exe`と`Eqnedit64-<version>.release.json`という名前で
`eqnedit64-staging` releaseへもアップロードする。

**O:はCIの入力ではない。** runner serviceは`NETWORK SERVICE`として動くため、
対話ユーザーのドライブ文字も、ワークグループ構成の研究室SMB共有も見えない
（2026-09-05に実測。`\\163.51.64.136\work`はrunnerからaccess denied）。よって
tag CIは`eqnedit64-staging` releaseのassetから署名済みEXEとマニフェストを取得し、
GitHub-hosted windowsで検証する。O:はあくまで人手ハンドテストの入口として残す。

tag CIはマニフェストに記録されたsource SHAがtag SHAと完全一致し、EXEの
SHA-256・ProductVersion・`CN=ksugahar`署名がすべて一致する場合だけ、その署名済み
EXEからPython 3.10--3.13用wheelを作り、PyPIへ公開した後にGitHub Releaseへ
`Eqnedit64.exe`と`SHA256SUMS.txt`を添付する。staging更新前にtagをpushしてはならない
（pushした場合はassetが無いため、待機ではなく明示的なエラーで停止する）。旧版へ戻す
必要が生じた場合も、対象Gitコミットを再ビルドして同じ順序と最終ゲートを通す。

## 正本移行後の品質運用

IME未確定文字列と候補のキャレット追従座標は自動回帰で固定する。人間にしか判断
できない候補の美観や長時間入力の感触は、完成前の曖昧な一回確認ではなく、正本
運用中も継続して観察する。違和感の直後に `F12` を記録し、
重大・高優先度の問題が見つかった場合は回帰試験を追加して修正し、再び最終ゲートを
通す。現時点で既知の重大・高優先度の未解決項目はない。
