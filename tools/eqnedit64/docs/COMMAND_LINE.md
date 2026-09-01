# Eqnedit64.exe コマンドライン リファレンス

GUI エディタだが、クリップボード変換・ファイル変換・自己検査はすべて
コマンドラインから起動できる。本書はソース
（`src/eqnedt64_app.cpp` の引数解析部）と一対一の全フラグ一覧。

## まず知っておく約束事

- **GUI サブシステムの実行ファイル**なので、起動してもコンソールに
  自動接続しない。標準出力を取るにはリダイレクトが必要：
  `Eqnedit64.exe --menu-list > menus.txt`。PowerShell で終了コードを
  確実に取るには `&` 呼び出しではなく
  `(Start-Process .\Eqnedit64.exe -ArgumentList '--self-test' -PassThru -Wait -WindowStyle Hidden).ExitCode`
  を使う（`$LASTEXITCODE` は空になることがある）。
- `--menu-list` と `--menu-responds-all` の標準出力は **UTF-16**
  （`_O_U16TEXT`）。PowerShell で読むときは `-Encoding Unicode`。
- **終了コード 0 が成功**。0 以外の番号はソースの `return <番号>;` と
  一対一なので、検索すれば失敗箇所に直行できる。番号は再利用しない規約。
- `--` で始まらない最初の引数は**起動時に開く `.tex` ファイル**。
- TeX 入力（引数・ファイル・クリップボードとも）は貼り付けと同じ正規化を
  通る：外側の `\[...\]`、`$$...$$`、`equation` / `equation*` 環境などの
  外皮は剥がされる。正規化後に空なら終了コード 83。

## 利用者向け：クリップボードへの変換

| コマンド | 動作 |
|---|---|
| `--texclip` | **クリップボード上の TeX を読み、300 dpi PNG（+ DIBV5）に置き換える**。TeXclip の代替。空・空白のみなら 83、書き込み失敗は 84 |
| `--clipboard-tex-to-png` | `--texclip` の長い別名（同一動作） |
| `--copy-tex "E=mc^2"` | TeX を **Office 用の全部入り**でクリップボードへ：編集可能 Office Math（MathML＋後置 NBSP の CF_HTML）+ LaTeX テキスト + EMF + 不透明 DIBV5。Word / PowerPoint に Ctrl+V でネイティブ数式になる |
| `--copy-tex-file equation.tex` | 同上、TeX をファイルから読む |
| `--copy-google-slides "E=mc^2"` | Google Slides 用：300 dpi / 24 pt の PNG + HTML |
| `--copy-google-slides-file equation.tex` | 同上、ファイル入力 |
| `--copy-png "E=mc^2"` | **画像だけ**（PNG + DIBV5）をクリップボードへ。数式を受け付けないアプリ向け |
| `--copy-png-file equation.tex` | 同上、ファイル入力 |

失敗コード：82 = 入力ファイルが読めない、83 = 正規化後の TeX が空、
84 = クリップボード書き込み失敗、94 = 引数不足。

## 利用者向け：ファイルへの変換

| コマンド | 動作 |
|---|---|
| `--render-png "TeX" out.png` | 実カンバスと同じ GDI 経路で描画して PNG 保存。成功時、標準出力に `幅 高さ`（ピクセル）を1行出す |
| `--render-emf "TeX" out.emf` | クリップボードに載せるものと同じ拡張メタファイルを保存（グリフはアウトライン化済み、フォント不要で再生可能） |
| `--render-png-file in.tex out.png` | 同上、TeX をファイルから |
| `--render-emf-file in.tex out.emf` | 同上、TeX をファイルから |

失敗コード：82 / 83 / 94 は上と同じ。95 = 描画失敗、96 = 出力ファイル
作成失敗、97 = 書き込み不完全、98 = PNG データ取得失敗。

## 情報表示

| コマンド | 動作 |
|---|---|
| `--version` | 版・ビルドスタンプ（コミット / 日時）・実行ファイルのフルパスを**ダイアログで表示**して終了。ヘッドレス運用でスクリプトから版を取りたい場合はダイアログが出る点に注意 |

## 開発者向け：自己検査（終了コードで判定）

| コマンド | 動作 |
|---|---|
| `--self-test` | モデル層の自己検査。出力なし・終了コードのみ |
| `--ui-interaction-test` | 実キー入力・実 WM_COMMAND 経路の対話検査。**終了コード＝失敗した検査の番号**（`src/eqnedt64_app.cpp` を `return <番号>;` で検索） |
| `--operation-test out.tex` | ヘッドレスで編集操作列（整列タブ・複数行・スロット移動・TeX 貼り付け）を実行し、保存結果を out.tex へ。20 = 引数不足、21–27 = 各段階の失敗、28 = 書き込み失敗 |
| `--operation-log-test out.log` | 操作ログ（v2 形式：caret / latex / elapsed_ms / focus / input_style …）の記録・読み戻し検査。30 = 引数不足、31–36 = 各失敗 |
| `--clipboard-publish-test` | 選択なし Ctrl+C と同じコマンド経路でコピーし、CF_ENHMETAFILE と CF_DIBV5 が実際に載ったか検査（81 = 失敗） |
| `--google-slides-clipboard-test` | Google Slides 形式の書き込み検査（82 = 失敗） |
| `--status-layout-test` | ステータスバー配置の検査 |
| `--visual-scale-test` | 4 段階のデバイススケールで描画し高 DPI 受け入れの機械判定（ウィンドウ非表示） |
| `--paint-bench` | フルリペイントの所要時間計測。1 ケース 5 ms 超で失敗 |

## 開発者向け：メニューとファズ

| コマンド | 動作 |
|---|---|
| `--menu-list` | 全メニュー項目を `ID<TAB>ラベル` で列挙（UTF-16 出力） |
| `--menu-press <id>` | 指定 ID を 1 プロセスで 1 回押して返る（ハング検出は外側のタイムアウトで） |
| `--menu-audit` | 全項目をディスパッチし「効果が無かった項目」を報告 |
| `--menu-responds-all` | **全メニューを 1 プロセスで**順に押す。項目ごとに `BEGIN<TAB>id<TAB>ラベル` を flush してから実行するので、外側タイムアウトでも止まった項目名が分かる。旧方式（253 プロセス起動）はフォント登録の繰り返しでセッションのフォント基盤を不安定化し得たため廃止 |
| `--ui-fuzz <seed> <ops>` | シード付きランダム操作を実ウィンドウプロシージャ経由で流す。再現は同じ 2 引数で |
| `--ui-fuzz-batch <seeds> <ops>` | シード 1..seeds を **1 プロセスで**連続実行（進捗は stderr）。最初の失敗コードで停止。64 = 引数不正 |

## デバッグ起動オプション（GUI と併用）

| オプション | 動作 |
|---|---|
| `--debug-operations` | 操作ログを記録しながら通常の GUI を起動 |
| `--debug-colors` | 領域を判別色で塗る（カンバス緑・数式箱黄）。スクリーンショットの機械判定用 |

## 環境変数

| 変数 | 効果 |
|---|---|
| `EQNEDIT64_NO_FONT_REG` | 定義されているとフォント登録を一切行わない。描画は代替フォントになり実用不可。フォント汚染切り分け実験（登録あり／なしで同一のプロセス起動列を比較する）専用 |

## 典型レシピ

```powershell
# クリップボードの TeX を画像化（TeXclip 代替）
.\Eqnedit64.exe --texclip

# 式を PowerPoint 用にコピーしてから貼り付け
.\Eqnedit64.exe --copy-tex "\frac{\partial \vec{D}}{\partial t}"

# バッチで数式画像を量産
Get-ChildItem *.tex | ForEach-Object {
  .\Eqnedit64.exe --render-png-file $_.FullName ($_.BaseName + '.png')
}

# ファズの一括実行と、失敗シードの単独再現
.\Eqnedit64.exe --ui-fuzz-batch 24 3000
.\Eqnedit64.exe --ui-fuzz 14 3000
```
