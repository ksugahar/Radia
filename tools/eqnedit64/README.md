# Eqnedit64 — 数式入力に特化したTeXエディタ

Eqnedit64は、軽快な構造編集とTeXファイルを直接つないだ64-bit Windows
数式エディタです。数式は常にTeXとして開き、編集し、保存します。
`Eqnedit64.exe` 単体で動作するポータブルアプリで、インストールや
レジストリ登録は不要です。

公開ソースはRadiaリポジトリの`tools/eqnedit64`、署名済み単体EXEは
[Eqnedit64 GitHub Release](https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.11)
で配布します。旧Eqnedit32バイナリ、MTEF変換コード、逆アセンブリ資料は
Eqnedit64のソース・ビルド・配布物に含めません。

ブラウザ版も [`web/`](web/) を正本として同じディレクトリで保守します。
研究室ホームページは現在の公開先ですが、ホームページ側のコピーを正本には
しません。共通の操作・貼り付け・TeX学習契約は
[`docs/PRODUCT_PARITY.md`](docs/PRODUCT_PARITY.md) で管理します。

GUIの画面構成、状態遷移、TeX/クリップボード契約、自動・手動の受入条件は
[`docs/GUI_SPEC.md`](docs/GUI_SPEC.md) を基準とします。
自動試験の合格と製品の完成を混同しないため、未完了項目は
[`docs/COMPLETION_GAPS.md`](docs/COMPLETION_GAPS.md) に明記します。

## 主な機能

- 分数、根号、上下付き、積分、総和、任意の長方形行列、場合分けを空欄へ順に入力
- 行列セル内から行・列を追加／削除し、1×1から99×99まで構造編集
- 構造キャンバスとTeXソースを同時表示し、どちらからもリアルタイム編集
- `Tab` / `Shift+Tab` で構造内の入力欄を移動
- `Enter` で複数行化、`&` で揃え位置を追加、上下キーで行を移動
- 入力中に数式が動かない左寄せ表示を既定とし、中央・右寄せへ即時切替
- ドラッグ選択、構造単位の切り取り・コピー・置換貼り付け
- PowerPointへ表示可能で編集できるOfficeネイティブ数式、不透明なEMF/DIBV5画像、TeXを同時に提供するコピー
- 300 dpi・24 pt基準でGoogle スライドへ貼る専用コピー（`Ctrl+Alt+C`）
- クリップボード上のTeXを300 dpi PNGへ置き換える `--texclip`
- Unicode文字列とTeX数式の貼り付け
- `sinx`の逐次入力を立体の`sin`と斜体の`x`へ自動判定
- TeX基準の積分上下限位置と、添字内で通常どおり働くBackspace
- 新規文書を`equation`環境のUTF-8 `.tex`として保存
  （既存の`equation*`文書を開いた場合はその外枠を維持）
- SVG書き出し
- 旧版互換ショートカットと、メニュー操作の直後に次のキーを教える
  ショートカット・コーチ
- 基本／解析／集合・記号／幾何／ギリシャの5分類タブで全パレットへ到達
- 分類を切り替えても常に見える `R x` / `I x` / `B x` から、立体・変数斜体・
  ベクトルを `\mathrm` / `\mathit` / `\mathbf` として適用
- パレット選択中は対応TeXだけを短く表示し、挿入後は常時表示ソースで綴りを学習

貼り付け時は `$...$`、`$$...$$`、`\(...\)`、`\[...\]`、
`equation` / `equation*` の外側だけを除きます。`aligned` は複数行構造として
保持します。プレーンテキストもそのまま入力できます。

キャンバスからのコピーは、一つの数式を複数形式でクリップボードへ登録します。
選択範囲があればその範囲を、選択がなければ数式全体を `Ctrl+C` でコピーします。
PowerPoint、Word、Excel向けには、inline 18 pt MathMLと18 ptを明示した末尾NBSPを
1つのCF_HTML断片として提供し、編集可能なOffice Mathへ変換します。Web版も同じ
経路を使います。PowerPointの通常Ctrl+Vでは左寄せを優先し、数式本体、末尾文字、
次の挿入点を標準18 ptでそろえます。
通常コピーでは、PowerPointが中央寄せの数式段落として優先する登録`MathML` /
`MathML Presentation`や、Web版だけの直接OMMLを混在させません。Office保存後は
OOXML内のインライン`m:oMath`となり、`m:oMathPara`にはなりません。
旧Office向けには区切り付きLaTeXも残します。IrfanViewなどの画像ソフトは
EMFまたは全画素不透明の32-bit DIBV5を選べます。TeX対応ソフト向けには生の断片も
`LaTeX` 形式で保持します。

「編集 → Google スライド用コピー」または `Ctrl+Alt+C` は、24 pt基準で描画した
300 dpi PNGと、同じ画像の表示寸法をpt単位で指定したHTMLをクリップボードへ
登録します。通常コピーと分けることで、PowerPointの通常貼り付けは引き続き
編集可能なOffice Mathを選びます。

クリップボード上のTeXをTeXclipと同じ流れでPNGへ置き換える場合は、次を実行します。
ウィンドウは表示せず、処理後すぐ終了します。

```powershell
Eqnedit64.exe --texclip
```

登録 `LaTeX` 形式があればそれを優先し、なければ通常のUnicode文字列を読みます。
`$...$`、`\[...\]`、`equation` などの外側は通常の貼り付けと同じ規則で除き、
白背景の300 dpi PNG（24 pt基準）でクリップボードを置き換えます。登録 `PNG` を
正本とし、IrfanViewなどのWindowsアプリ用に白背景・黒文字・全画素α=255の
`CF_DIBV5` も併記します。
成功時の終了コードは0、入力が空なら83、画像化またはクリップボード登録に失敗した
場合は84です。長い別名 `--clipboard-tex-to-png` も同じ動作です。

### コマンドラインAPI

GUIを開かず、任意TeXをPowerPoint/Word、Google Slides、画像ファイルへ渡せます。
長い式や自動処理では、引用符やWindowsのコマンドライン長に依存しないUTF-8
ファイル入力を推奨します。

```powershell
# PowerPoint/Word: 編集可能Office Math + TeX + EMF + 不透明DIBV5
Eqnedit64.exe --copy-tex-file equation.tex

# Google Slides: 300 dpi / 24 ptのPNG + HTML
Eqnedit64.exe --copy-google-slides-file equation.tex

# 画像だけをクリップボードへ
Eqnedit64.exe --copy-png-file equation.tex

# ファイルへ画像化
Eqnedit64.exe --render-png-file equation.tex equation.png
Eqnedit64.exe --render-emf-file equation.tex equation.emf
```

短い式は `--copy-tex "E=mc^2"`、`--copy-google-slides "E=mc^2"`、
`--copy-png "E=mc^2"`、`--render-png <TeX> <出力>`、
`--render-emf <TeX> <出力>` でも渡せます。外側の `\[...\]`、`equation` 等は
GUI貼り付けと同じ規則で除きます。入力ファイルを読めない場合は82、空入力は83、
クリップボード登録失敗は84、引数不足は94です。

## 操作デバッグ

「ヘルプ → 操作ログを記録」または次のコマンドで、利用者自身の操作を
UTF-8ログに記録できます。外部からマウスやキーボードを操作しません。

```powershell
cmd /c build\run_operation_debug.bat
```

不具合に気づいた直後に `F12` を押すと、ログへ目印が残ります。ログは
`%LOCALAPPDATA%\Eqnedit64\logs` に保存されます。入力した数式も含まれるため、
共有範囲には注意してください。詳しい確認順序は
[`docs/OPERATION_DEBUG_GUIDE.md`](docs/OPERATION_DEBUG_GUIDE.md) にあります。

記録後は次のコマンドで、最新ログからLLM操作感レビュー束を作れます。

```powershell
pwsh -NoProfile -File build\analyze_last_operation_log.ps1
```

既定では数式本文を伏せ、F12、即時Undo、移動の反転、未登録キー、no-op、
修正の集中を証拠列として抽出します。LLMは改善候補と非表示回帰試験を提案し、
利用者が採否を決めます。採用済みの判断は機械可読な選好台帳へ蓄積します。
設計と判断規則は [`docs/LLM_USABILITY_ENGINE.md`](docs/LLM_USABILITY_ENGINE.md)
を参照してください。この分析器は開発・研究用であり、配布exeの実行には
PythonもLLM接続も必要ありません。

## 操作

パレットバーの1段目は「基本／解析／集合・記号／幾何／ギリシャ」の分類タブ、
2段目は選択中の分類に属するパレットだけです。19パレットを一度に探す必要はなく、
全セルへの到達性と既存ショートカットは維持しています。メニューバーも
「挿入」1本の同じ5分類です。

2段目右側の `R x` / `I x` / `B x` は分類共通なので常時表示します。キャンバスで
範囲を選んで押すとその場で書体を変更し、未選択なら以後の入力スタイルになります。
TeXソースでは選択範囲を包み、未選択なら空の `{}` の内側へ入ります。保存には旧式の
`\rm` / `\it` / `\bf` ではなく `\mathrm{...}` / `\mathit{...}` /
`\mathbf{...}` を使います。
追加書体は基本→装飾パレットから `\mathsf`、`\mathtt`、`\mathcal`、
`\mathbb`、`\mathfrak`、`\bm` を選べます。`\boldsymbol` は読み込み時に
同義の `\bm` へ正規化し、`\mathnormal` は標準数式書体へのリセットとして
読み込めます。`\mathbb` / `\mathfrak` を外部文書で使う場合は通常
`amsfonts` または `amssymb`、`\bm` には `bm` パッケージが必要です。

画面下の「TeXソース」は常時表示され、上の構造キャンバスと同じ数式を
双方向に編集します。TeXを1文字変更するたびキャンバスへ反映され、閉じ括弧を
まだ打っていない途中状態も入力欄から消しません。キャンバスへ戻るか保存すると、
構造モデルの正規化済みTeXへ揃います。ソースでの連続入力はまとめて1回のUndoです。
TeXソースは常時表示します。キャンバスまたはTeXソースをクリックすれば、
その場で直接編集できます。GUI／TeXの入力モード切替はありません。
TeXコマンドはソース欄へ直接入力します。キャンバスには
`\command`＋空白確定という別の入力モードを持たせません。
どちらの入力面でも `Tab` / `Shift+Tab` が次／前の空欄へ移動します。TeXソースでは
空の `{}` の内側を順にたどるため、パレットで構造を入れた直後も同じ手癖で埋められます。

パレットのセルを選ぶとステータス欄に対応する正規化TeXを表示し、挿入後は
常時表示のTeXソースで新しく入った綴りだけを強調します。構造ショートカットでも
同じ強調を出しますが、入力フォーカスとキャレットはキャンバスに残るため連続入力を
止めません。TeXソースを編集し始めると強調は解除され、誤って置換しません。
日本語説明を重ねず、記号とTeXを直接対応させます。ショートカットがある項目は
TeXの後ろに次回のキーも表示します。
ヘルプメニューの
「ショートカット・コーチ」でいつでも有効／停止を切り替えられます。
完全な一覧は [`docs/SHORTCUTS.md`](docs/SHORTCUTS.md) を参照してください。
編集キャンバスは既定で左寄せです。これは入力位置を安定させる表示設定だけで、
保存したTeXの `equation` 環境はLaTeX標準どおり中央に組版されます。
複数行も、`&`を置かない間は各行の左端をそろえて入力できます。`&`を入力すると
その位置を境にTeX標準の右列／左列整列へ切り替わり、保存ソースにもその`&`だけを
残します。
数学入力で`sinx`と続けて打つと、`sin`を認識した時点で立体へ固定し、後続の`x`だけを
変数斜体にします。保存TeXは`\sin x`です。積分上下限もTeXのdisplay積分と同じ
右上・右下位置へ配置します。

| キー | 動作 |
|---|---|
| `Ctrl+F` | 分数 |
| `Ctrl+R` | 根号 |
| `Ctrl+L` / `Ctrl+H` / `Ctrl+J` | 下付き / 上付き / 両方 |
| `Ctrl+I` / `Ctrl+9` / `Ctrl+0` | 積分 / 丸括弧 |
| `Ctrl+[` / `Ctrl+]` / `Ctrl+{` / `Ctrl+}` | 角括弧 / 波括弧 |
| `Ctrl+T,S` / `Ctrl+T,P` / `Ctrl+T,M` | 総和 / 総乗 / 3×3行列 |
| `Ctrl+Alt+↓` / `Ctrl+Alt+↑` | 現在セルの下へ行を追加 / 現在行を削除 |
| `Ctrl+Alt+→` / `Ctrl+Alt+←` | 現在セルの右へ列を追加 / 現在列を削除 |
| `Ctrl+K` の後に文字 | 数学記号（例: `I` → ∞、`D` → ∂） |
| `Ctrl+G` の後に文字 | ギリシャ文字（例: `A` → α） |
| `Ctrl+B` の後に英字 | ベクトル太字 |
| `Enter` / `&` | 改行 / 揃え位置 |
| `Tab` / `Shift+Tab` | 次 / 前の入力欄 |

## ビルド

Visual Studio 2022 C++ Build ToolsとPython 3.12を使用します。

開発機ごとに最初の1回だけ、コード署名専用の開発者証明書を準備します。
秘密鍵はCurrentUserストア内の非エクスポート鍵です。

```powershell
pwsh -NoProfile -File build\setup_developer_signing.ps1
```

すべてのビルドとバックグラウンド試験を一括実行できます。
ビルド時にリポジトリ内の `build` / `dist` 版が起動中なら、差し替えのため
そのプロセスは強制終了します。ビルド前に必要な文書を保存してください。

3.0.2以降は、内蔵数式フォントを検証済みユーザーキャッシュへ展開し、
ファイルベースで私有登録します。旧メモリ登録で発生したWindowsフォントホストの
クラッシュを再発させないため、短命プロセスの連続起動を含む総合試験は
GitHub-hosted Windows CI、VM、または専用Windowsユーザーセッションだけで実行し、
普段使う対話セッションでは実行しません。隔離済みセッションでは次のマーカーを
設定します。

```powershell
$env:EQNEDIT64_ISOLATED_TEST_SESSION = '1'
cmd /c build\build_tests.bat
```

リリース候補の通常回帰、ASan、拡張GUI耐久、署名・単体配布を一括判定する
最終ゲートは次です。公開配布は、このゲートが生成した `dist\Eqnedit64.exe`
だけをGitHub Releaseへ添付します。

```powershell
$env:EQNEDIT64_ISOLATED_TEST_SESSION = '1'
pwsh -NoProfile -File build\accept_release.ps1
```

菅原のハンドテストに渡す候補は、署名・ローカルバックグラウンド試験後に
`sync_handtest_to_o.ps1`で `O:\Eqnedit64.exe` へ先に同期します。O:が正規の
ハンドテスト入口であり、`C:\temp`だけに候補を置いてテストを依頼しません。

リリース時は、リリースコミットを先に `main` へpushしてCIを通し、同じ
`origin/main`から作った署名済みEXEを `O:\Eqnedit64.exe` と
`O:\Eqnedit64.release.json` へ同期します。その成功後にだけ
`eqnedit64-v<version>` tagをpushします。tag CIはO:上のEXEとtag SHAの一致を
検証し、PyPI、GitHub Releaseの順に公開します。具体的な手順は
[`release-eqnedit64` skill](../../.agents/skills/release-eqnedit64/SKILL.md)と
[`CANONICAL_OPERATION.md`](docs/CANONICAL_OPERATION.md)を正本とします。

個別に実行する場合:

```powershell
cmd /c build\build_eqnedt64.bat
cmd /c build\build_pymodule.bat
build\test_tex_document.exe
python tests\run_model_tests.py
pwsh -NoProfile -File build\test_background.ps1
pwsh -NoProfile -File build\test_ui_fuzz.ps1
pwsh -Sta -NoProfile -File build\test_external_paste.ps1
```

`run_model_tests.py` は数式フォントの私有登録・解除を繰り返さないよう、7つの
ヘッドレスモデル試験を1つのPythonプロセスで実行します。各試験は個別にも
実行できます。

メモリ破壊を検査するASan版も、同じ非表示GUIファズへ渡せます。

```powershell
cmd /c build\build_asan.bat
pwsh -NoProfile -File build\test_asan.ps1
```

配布物は `dist\Eqnedit64.exe` の1ファイルだけです。VC++ランタイムは
静的リンクされ、Pythonや外付けDLLを必要としません。アイコンのMIT表示も
「ヘルプ → 第三者ライセンス」としてexe内に収録しています。ビルドは
`ksugahar` のAuthenticode署名が `Valid` に
ならなければ失敗します。これは研究室内の自己署名開発者証明書であり、公開CAの
信頼やSmartScreen評価を得るための商用コード署名証明書とは別です。

`test_background.ps1` はウィンドウを表示せず、貼り付け、構造編集、複数行保存、
18 pt inline MathML CF_HTML、Office認識用LaTeX、EMF/DIBV5生成、x64実行形式、
製品バージョン、開発者署名を検査します。
マウスやキーボードは操作しません。

GUIも、最初から人の画面を自動操作するのではなく、非表示のWin32ウィンドウで
コントロール種別、矩形の重なり、ステータスバー区画、フォント高、表示モードを
検査します。さらに `WM_CHAR`、`WM_KEYDOWN`、`WM_COMMAND`、TeXソースの
`EN_CHANGE`を直接送り、文字入力、未完成TeX、双方向同期、1回のUndoへまとまる
連続ソース編集、テンプレート、常時表示・焦点切替をデスクトップに触れず検査します。
GUIファザーにも完成／未完成TeXのソース編集を混ぜます。数式描画はメモリ上のGDI画像で
検査し、PowerPointなど他製品との貼り付けも `test_external_paste.ps1` が非表示の
API経由で検査します。この試験は内部生成関数を迂回路として使わず、選択なしの
通常GUIコピー命令を送ってから、PowerPointの組み込みUI `Paste` が
左寄せ18 ptの編集可能なOffice Mathと分数・根号の構造を作ることに加え、貼り付けた
図形をPowerPoint自身でPNG化し、空白や豆腐文字ではなく分数線を含む数式の輪郭が
描画されること、
IrfanViewの `/clippaste` が非空画像を作ること、DIBV5の全画素がα=255であること、
登録 `LaTeX` 形式のraw断片を確認します。
Google スライド用コピーは、登録PNGとHTML内PNGの同一性、300 dpiメタデータ、
24 pt基準の物理寸法を同じ試験で確認します。
試験前のクリップボードは全形式を複製し、終了時に復元します。PowerPointと
IrfanViewがインストールされた試験機が必要ですが、マウス、キーボード、前面
ウィンドウは操作しません。既存のPowerPointが開いている場合は接続せず試験を中止します。
人間の手動確認には文字の細かな見た目と操作感だけを残します。

`test_tex_fuzz.py` は固定seedから正常・壊れかけのraw TeXを3,072件作り、
既定24 pt MathMLのXML妥当性を含めて、
パース、正規化の固定点、有限な描画寸法、SVG XMLを検査します。失敗はseedで
完全に再現できます。
操作感CIは実exeのoperation log v2を読み、構造プライバシーモードのLLMレビュー束、
摩擦候補、F12周辺の証拠列、選好台帳を検査します。解析中もデスクトップを操作せず、
ネットワークへログを送信しません。
リリース試験は `dist` にexe以外がないこと、ビルド済みexeとSHA-256が
一致すること、VC/Python実行時依存がないこと、配布先から非表示起動できることを
検査します。

## ライセンス

Eqnedit64本体はBSD 2-Clause Licenseです。詳細は[`LICENSE`](LICENSE)を参照して
ください。内蔵Latin Modern Mathフォントとアプリケーションアイコンの第三者素材は
[`assets/THIRD_PARTY_NOTICES.md`](assets/THIRD_PARTY_NOTICES.md)および
[`assets/GUST-FONT-LICENSE.txt`](assets/GUST-FONT-LICENSE.txt)に従います。
