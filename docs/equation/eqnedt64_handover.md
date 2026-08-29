# Eqnedit64 引き継ぎ正本

- 文書状態: 現行
- 対象製品: Eqnedit64 native / Python package / Web editor
- 対象リリース: 3.0.3（公開済み）
- 基準日: 2026-08-29
- リポジトリ: `ksugahar/Radia`

この文書は、Eqnedit64を修正、試験、配布するときに最初に読む引き継ぎ正本である。
過去の `src/ext/equation`、MTEF変換、旧Radiaバージョン、旧インストーラーを前提にした
引き継ぎ文書は廃止する。詳細な画面契約は
`tools/eqnedit64/docs/GUI_SPEC.md`、キー一覧は
`tools/eqnedit64/docs/SHORTCUTS.md`、native/Webの対応方針は
`tools/eqnedit64/docs/PRODUCT_PARITY.md` を参照する。

## 1. 製品の目的と正本

Eqnedit64は、Microsoft 数式3.0（Eqnedit32）の軽さと構造入力の手癖を、標準TeX、
64-bit Windows、Officeの編集可能な数式貼り付けへ移行する数式入力専用ツールである。
長い文書を組むTeXエディタではなく、数式を止まらず入力し、他の文書へ安全に渡すことに
集中する。

研究室の数式資産の正本はUTF-8 TeXである。

- 新規保存は `.tex` の `equation` 環境とする。
- 複数行は `equation` 内の `aligned` で表す。
- `.eqn` とMTEFは新規作成、保存、受け渡しの対象にしない。
- 既存MTEF資産が必要なら、一回限りの移行作業でTeXへ変換し、以後はTeXを正とする。
- Eqnedit32のバイナリや逆アセンブリ資料は調査済み資料であり、Eqnedit64の実行依存にしない。

Eqnedit64が正となった後、旧「数式3.0」は日常運用、関連付け、配布から退役させる。
旧資産をEqnedit64の公開ソースや配布物へ混ぜてはならない。

## 2. 現行の製品面

Eqnedit64は一つの製品系列に三つの操作面を持つ。

| 操作面 | 正本の場所 | 用途 |
|---|---|---|
| native EXE | `tools/eqnedit64` | 高速な構造GUI、ショートカット、TeX保存、Office/画像クリップボード、CLI |
| Python package | `packages/eqnedit64` | `pip install eqnedit64` で自動処理、描画、クリップボード、同梱EXEを提供 |
| Web/JS | `tools/eqnedit64/web` | 研究室ホームページ上のTeX入力、パレット、リアルタイム表示、Office貼り付け |

研究室ホームページへ置かれたWebファイルは複製先であり正本ではない。Web版もRadia側で
保守し、ホームページはRadiaの正本を参照または同期する。

### 2.1 描画エンジンの境界

- native EXEはKaTeXでもMathJaxでもない。C++の構造モデルと独自レイアウト／GDI描画を使う。
- Web版はサイト共通のMathJaxを使う。KaTeXが動いているかをnative表示の診断理由にしない。
- 同じTeXに対し、完全に同じ内部実装ではなく、利用者に同等の構造と貼り付け結果を提供する。
- 共通機能を変えたらnativeとWebの両方を同じ変更内で確認する。適用しない場合は理由を
  `PRODUCT_PARITY.md` に残す。

## 3. 配布契約

### 3.1 standalone EXE

- `Eqnedit64.exe` 一つで起動する。
- Python、外付けDLL、インストーラー、Inno Setup、ActiveX、OLE登録を要求しない。
- `.eqn` のファイル関連付けを新設しない。
- 配布EXEは `CN=ksugahar` の有効なAuthenticode開発者署名を持つ。
- アプリのアイコン、バージョン情報、Latin Modern Math、第三者ライセンスをEXEへ収録する。
- 起動引数に `.tex` を渡した場合はそのファイルを開く。

### 3.2 Python package

`pip install eqnedit64` はWindows wheelとしてPython API、CLI、同じ署名済みEXEを配置する。
pipは通常のPythonパッケージ配置を行うだけで、GUIインストーラーを起動しない。
standalone利用者はpipを必要としない。

対応wheelはCPython 3.10、3.11、3.12、3.13の`win_amd64`である。Eqnedit64はRadiaの
モノレポに属するが、`radia`、`radia-mcp`、`cubit-mesh-export`、`radia-optuna`とは
独立した配布単位・CI単位として扱う。全体CIをEqnedit64の合否代わりにしない。

## 4. native GUIの操作契約

### 4.1 二刀流編集

数式キャンバスとTeXソース欄を常時同時に表示する。F6などのモード切替は設けず、
クリックした側をそのまま編集できる。

- キャンバスは構造スロット、選択、キャレット、Backspace、テンプレート挿入を扱う。
- TeXソース欄はraw入力を保ちながらリアルタイムに構造モデルへ反映する。
- ソース欄を単なる「TeXを表示」機能として隠すメニューは設けない。
- キャンバスに `\command` を入力して空白で確定する旧式のTeXコマンド入力は持たない。
- TeXコマンドはパレットから選び、挿入直後のTeX範囲を短時間強調して手癖として学べるようにする。
- 日本語による長い説明はパレットへ載せない。表示するのは記号、構造、TeX、短いキー案内である。

数式キャンバスの既定位置は左寄せ、上寄せである。中央寄せは保存形式やOffice貼り付けを
表す設定ではなく、画面上の表示設定に限る。

### 4.2 パレット

Web版の学習順をnativeへ取り込み、上段を次の5分類に固定する。

1. 基本
2. 解析
3. 集合・記号
4. 幾何
5. ギリシャ

下段に選択分類のパレットを最大5個表示し、格子ポップアップから2クリック以内で
すべての収録記号・テンプレートへ到達できるようにする。現在の分類タブとボタン配置は
採用済みであり、今回の文字表示修正を理由に旧ツールバーへ戻さない。

分類共通の書体 `R x` / `I x` / `B x` は下段右側の固定位置へ常時表示する。nativeでは
キャンバス選択へ適用するか、未選択時の継続入力スタイルを切り替える。TeXソースと
Web版では選択を包むか、空の `{}` を挿入して内側へキャレットを置く。正規形は
`\mathrm{...}` / `\mathit{...}` / `\mathbf{...}` であり、旧宣言形へ戻さない。

行列は2×2、3×3だけに制限しない。パレットには頻用の1×2から6×6を置き、
行／列の追加・削除で1×1から99×99まで編集できる。ポップアップの操作セルは
Latin Modern Mathでも確実に読める短いASCII主体の `+R`、`−R`、`+C`、`−C` を使う。
日本語の「行」「列」を数式フォントで描かない。

空欄は埋め込みフォントに無い`▯`（U+25AF）ではなく、同フォントが所有する
`□`（U+25A1）で統一する。n乗根のfaceもフォントに無い上付きnではなく`n√□`とする。

上下付きの空欄見本に数字を使わない。`₅`のような見本は実際の入力値5と誤認されるため、
上、下、上下を示す中立的な矢印と空欄（`□↑`、`□↓`、`□↕`）で表す。

### 4.3 ショートカットとスロット移動

Eqnedit32から回収した操作系列は互換層として維持する。追加キーは拡張層とし、
既存キーを理由なく変更しない。完全な一覧は`SHORTCUTS.md`を正とする。

重要な契約は次のとおり。

- `Tab` / `Shift+Tab`: 構造GUIでは次／前の入力穴、TeXソースでは次／前の空の`{}`へ移動。
- `Enter`: `aligned`の改行。
- `&`: 整列位置。
- 選択中の`Ctrl+B`: 選択を `\mathbf{...}` にする。
- `Ctrl+B`の後に英字: 次の英字をベクトル太字として入力する。
- `Ctrl+Shift+B`: 行列・ベクトル入力スタイル。
- `Ctrl+Alt+-`: 直前の項目へベクトル矢印を付ける。
- `Ctrl+Alt+↓/↑`: 行列の行追加／削除。
- `Ctrl+Alt+→/←`: 行列の列追加／削除。
- `Ctrl+Alt+C`: Google Slides用の300 dpi / 24 ptコピー。

ショートカット・コーチは次に使えるキーを短く提示する。学生にコマンド表の暗記を
要求するのではなく、パレット、直前挿入TeXの強調、反復操作で自然に覚えられる設計とする。

### 4.4 Backspace

Backspaceは一打鍵で利用者に見える一項目だけを削除する。構造の中身が残る間は親構造を
まとめて消さない。例えば `E=mc^{2}|` は `E=mc^{}|`、次に `E=mc|` となる。
上付きの直後、空上付き、母項を持つ上付き、入れ子をそれぞれ別の回帰ケースとして試験する。

## 5. TeX保存契約

- 保存文字コードはUTF-8。
- 新規文書は `\begin{equation}` / `\end{equation}` で保存する。
- 複数行は内側の `\begin{aligned}` / `\end{aligned}` と `\\` で表す。
- 既存の `equation*` を開いた場合だけ、その外枠を保存時に維持する。
- 「番号あり／番号なし」を通常操作で切り替えるメニューは設けない。
- 表示の左／中央／右寄せは保存TeXを変えない。
- `.eqn`、MTEF、MathType固有形式を保存候補に復活させない。

## 6. クリップボードと外部アプリ

### 6.1 通常コピー

通常コピーは同じ数式から、少なくとも次の利用経路を提供する。

- UTF-8/Unicode TeX
- 不透明なEMF
- 不透明なDIBV5/PNG画像
- Officeが編集可能な数式へ変換できるMathML/HTML/登録形式

PowerPointとWordでは、空のテキストボックスや画像への退化ではなく、見えて編集できる
Office Mathになることが合格条件である。MathMLの文字列がクリップボードに存在するだけ、
図形数が1以上だけ、例外が出ないだけでは合格にしない。PowerPoint自身に貼り付けさせ、
数式の輪郭と代表構造を画像化して検査する。

Office向けの設計サイズは24 ptを基準にする。Office側が経路により18--24 ptへ再解釈する
ことはあり得るが、Eqnedit64自身が意図せず18 ptを生成する退行は試験で検出する。

### 6.2 Google Slidesと画像

Google Slides用コピーは300 dpiのPNGとHTMLを登録し、24 pt数式と同じ物理寸法になるように
する。画像は黒い数式が不透明に見え、透明すぎる、黒文字が消える、アルファが反転する状態を
不合格とする。

### 6.3 CLI

GUIを開かずに次の代表経路を使える。

```text
Eqnedit64.exe --copy-tex-file equation.tex
Eqnedit64.exe --copy-google-slides-file equation.tex
Eqnedit64.exe --copy-png-file equation.tex
Eqnedit64.exe --render-png-file equation.tex equation.png
Eqnedit64.exe --render-emf-file equation.tex equation.emf
Eqnedit64.exe --clipboard-tex-to-png
```

最後のコマンドはクリップボード上のTeXを読み、同じクリップボードをPNGへ置き換える。
CLI/APIはPowerPoint自動生成や`radia-mcp.presentation`から呼べる安定した境界とする。

## 7. フォント契約

Eqnedit64はLatin Modern MathをEXEへ埋め込む。実行時は内容ハッシュ付きのファイルを
`%LOCALAPPDATA%\Eqnedit64\fonts\latinmodern-math-<hash>.otf`へ展開し、
`AddFontResourceExW(FR_PRIVATE | FR_NOT_ENUM)`でそのプロセスだけに登録する。

`AddFontMemResourceEx`へ戻してはならない。旧3.0.1以前のメモリフォント登録／解除を
繰り返す経路は、Windows Server 2022の対話セッションで`fontdrvhost.exe`を壊し、
Eqnedit64だけでなくLINE等の文字を消す事故を起こした。3.0.2のファイル-backed private登録と
32回ライフサイクル試験が再発防止策である。

文字種ごとの描画責任を混同しない。

- 数式グリフ: Latin Modern Math。
- メニュー、分類タブ、日本語UI: 明示的なCJKシステムフォント。
- 数式フォントで描くowner-drawパレットセル: そのフォントが持つ記号だけを使う。
- 日本語ラベルをLatin Modern Mathへ渡し、フォントリンクに期待してはならない。

古い版によって既にセッションのフォントホストが壊れた場合、修正版の起動だけでは
OSセッション内の壊れた状態を修復できない。一回のサインアウトまたは再起動が必要な場合が
ある。ただし3.0.2以降は新しい破損を起こさないことをCIで確認する。

## 8. 3.0.3で回収した表示不具合

3.0.2に残っていた次の三件は3.0.3で修正し、再発防止試験へ変換した。

### 8.1 `\vec{\mathbf{E}}`のEが太字に見えない

パーサは`\mathbf`内の文字を`TF_VECTOR`として保持し、MathML emitterも
`mathvariant="bold"`を出す。しかしnativeの`Node::kChar`描画が`TF_VECTOR`を無視して
通常のEを出していた。したがってKaTeXの問題ではない。

修正は、Latin、数字、Greekの対応文字をUnicode Mathematical Boldの設計済みグリフへ
写像して描く。対応グリフがない文字だけGDIのsynthetic boldへフォールバックする。
TeX出力が正しいだけでは合格にせず、plain `E` と `\mathbf{E}` の実描画ピクセルが異なり、
太字側のink量が増えることをhidden native testで検査する。

### 8.2 行列ポップアップの文字が黒い塊になる

`＋行`、`－行`、`＋列`、`－列`をowner-drawの数式フォントで描いていたため、CJKを持たない
Latin Modern Math上で欠落または黒い代替形になった。操作セルを `+R`、`−R`、`+C`、`−C`
へ変更する。フォント選択サンプルにも、空欄、矢印、R/C、プラス／マイナスを含める。

### 8.3 上下付きポップアップの空欄が数字5に見える

下付き位置の見本に実際のUnicode下付き数字`₅`を使用していたため、利用者には入力済みの5に
見えた。文字化けではなく不適切な見本である。`□↑`、`□↓`、`□↕`のような値を持たない
表現へ変更し、上下付きのfaceへ数字が混入しないことを試験する。

## 9. 自動試験の原則

「表示用データを生成できた」ことと「人間に正しく見える」ことを分ける。GUI、パレット、
画像、Office貼り付けは、最終ピクセルまたはOffice自身の描画結果まで検査する。

### 9.1 native/model/document tests

- C++ `ctest`: parser、model、renderer、emitters。
- `tools/eqnedit64/tests/run_model_tests.py`: Python拡張を用いたモデル試験。
- `test_layout.py`: display listとSVGの幾何、ink、太字等。
- `test_palettes.py`: 全記号・テンプレートの所属、重複、face契約。
- TeX compile、MathML/OMML、PNG/EMF/DIBV5、fuzz、性能試験。
- `--self-test`、`--visual-scale-test`、`--ui-interaction-test`、`--menu-audit`等の
  非表示native経路。

3.0.3では次を必須の負の回帰試験にする。

- 旧実装の通常E表示なら落ちる、実ピクセルのbold試験。
- 上下付きfaceに数字が含まれたら落ちる試験。
- 行列操作faceがCJKや欠落グリフへ戻ったら落ちる試験。
- 全ポップアップセルが背景と異なる読めるinkを持つ試験。単に数個のpixelがあるだけを
  可視性合格にしない。

### 9.2 GUI操作試験

ユーザーのマウス、キーボード、前面ウィンドウを奪わない。操作系列は非表示ウィンドウ、
off-screen bitmap、メッセージ送信、CLI、COM APIで試す。

- `--ui-interaction-test`は実際のキー経路、Backspace、選択、Tab、行列操作を通す。
- `--visual-scale-test`は96/120/144/192 dpi相当で文字・ボタン・数式のinkを測る。
- UI fuzzは保存、ダイアログ、利用者クリップボードを対象外にし、再現seedを残す。
- ステータスバーは区画位置と文字のclip/ずれを検査する。
- IMEの最終的な候補窓の感触など、人間しか判断できない項目だけを限定的な手動QAへ残す。

### 9.3 フォント試験の隔離

完全なフォントライフサイクル試験を対話中のLABで繰り返してはならない。
GitHub-hosted Windows runnerまたは破棄可能なVM/ユーザーセッションで
`EQNEDIT64_ISOLATED_TEST_SESSION=1`を設定する。

`build/test_font_session.ps1`は32回起動終了し、次を監視する。

- `fontdrvhost.exe`のPIDが途中で変わらない。
- Windows Application logへ新しいEvent 1000が増えない。
- 各起動のprivate font登録、描画、解除が成功する。

対話中LABでは単発のhidden testは実行できるが、環境変数を偽装して32回試験のguardを
迂回しない。

### 9.4 外部貼り付け試験

貼り付けを手動受入だけにしない。

- PowerPoint/WordはCOM/API経由で一時文書へ貼り、Office Mathとして編集可能であること、
  表示が非空で分数線等を含むことをOffice自身のexport画像から確認する。
- IrfanView等の画像経路はCLIまたはAPIでPNG/EMFを読み、寸法、alpha、黒inkを確認する。
- Google Slides用payloadはPNGとHTML内PNGの一致、300 dpi metadata、24 ptの物理寸法、
  不透明背景と黒inkを確認する。
- 外部試験は既存Office文書へ接続しない。既存PowerPointが開いている場合は安全側で中止する。
- 試験前のクリップボードを全形式で退避し、終了時に復元する。

## 10. CI境界

Eqnedit64専用CIは `.github/workflows/eqnedit64.yml` である。次のパス変更だけで起動し、
Radia全体CIの通過を代用しない。

- `tools/eqnedit64/**`
- `packages/eqnedit64/**`
- `radia-mcp.presentation`のEqnedit連携箇所
- Eqnedit64の二つのworkflow

Windows Server 2022 runnerでnative CMake build、wheel candidate、C++ tests、Web/MCP契約、
JS構文、model suites、hidden EXE checks、32回font-host gate、PNG/EMF artifactを実行する。
PR、main、release tagの各CIは同じ製品試験を通す。

CIを改善するときは、実際に報告された壊れた版が新テストで失敗することを確認する。
生成ファイルの存在、終了コード0、画素1個だけのような弱いassertionを追加して完了扱いにしない。

## 11. ビルドと署名

canonical scriptsは次のとおり。

```text
tools/eqnedit64/build/build_eqnedt64.bat
tools/eqnedit64/build/build_pymodule.bat
tools/eqnedit64/build/accept_release.ps1
```

native出力は `tools/eqnedit64/dist/Eqnedit64.exe`。ビルド前に同じ出力ファイルを
ロックしているEqnedit64プロセスがあれば、利用者の許可済み方針に従ってそのEqnedit64だけを
終了してよい。別の場所の同名EXE、Office、ブラウザ、LINEなどを巻き込まない。

リリースEXEには次を満たさせる。

- ProductVersionとpackage versionが一致。
- build stamp/source SHAが公開対象の正確な`origin/main` commitを指す。
- Authenticode signer subjectが厳密に`CN=ksugahar`。
- 署名後のSHA-256をO: manifest、wheel、GitHub Releaseで同一にする。

## 12. 公開手順（順序固定）

公開順は次のとおりで、入れ替えてはならない。

1. 修正、仕様、試験、version/changelogを一つのrelease commitへまとめる。
2. PR CIを通し、`main`へ統合する。
3. `main`をpushし、Eqnedit64専用main CIがgreenであることを確認する。
4. exact `origin/main`からLABでrelease EXEをビルドし、`CN=ksugahar`で署名する。
5. `accept_release.ps1`の最終受入を通す。
6. `.agents/skills/release-eqnedit64/scripts/sync_to_o.ps1`で
   `O:\Eqnedit64.exe`を更新する。
7. `O:\Eqnedit64.release.json`へ予定tag、version、source SHA、EXE SHA-256、signerを記録し、
   実物と一致することを確認する。
8. ここまで成功してから `eqnedit64-vX.Y.Z` tagをpushする。
9. tagのEqnedit64 CIを通す。
10. `.github/workflows/release-eqnedit64-pypi.yml`がO: gateを照合してwheelをPyPIへ公開し、
    GitHub Releaseへ署名済みEXEと`SHA256SUMS.txt`を添付する。
11. PyPIから各wheel、GitHub ReleaseからEXEを再取得し、同梱EXEとO:のEXEがbyte-identical、
    SHA-256一致、署名有効であることを外側から確認する。

O:は対話中LABでは`C:\Users\Administrator\OneDrive`へのSUBSTである。self-hosted runnerは
LocalSystemのためユーザー固有ドライブ文字O:を見られない場合がある。release workflowは
`Get-PSDrive O`を確認し、無い場合は同じOneDrive backing pathを検査する。O:更新前にtagを
pushするとrelease gateが失敗するのが正しい。

## 13. 現在の公開状態

Eqnedit64 3.0.3が2026-08-29時点の公開済み基準版である。

- product tag: `eqnedit64-v3.0.3`
- product tag source: `c15e626d28245ee2db5b688f52d7756e000780a4`
- O: / GitHub Release EXE SHA-256:
  `C6D88C59227718FFE8A9BDB435632EAC64F5709615520FBA04E2C0A5C77C8A8F`
- signer: `CN=ksugahar`
- GitHub Release: `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.3`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.3/`
- PyPI wheels: CPython 3.10、3.11、3.12、3.13の`win_amd64`計4個。
- O:、GitHub Release、CPython 3.12 wheelから再取得した同梱EXEはbyte-identicalで、
  すべて上記SHA-256と有効な`CN=ksugahar`署名を持つ。

3.0.2は安全なfile-backed private font登録を導入した前版である。3.0.3はその安全策を維持し、
8章のvisual regressionと、同じ種類の見落としを防ぐ全パレット実描画試験を同梱する。

## 14. ソース地図

| 場所 | 責任 |
|---|---|
| `tools/eqnedit64/src/equation_edit.cpp` | 構造編集、キャレット、選択、Backspace、行列編集 |
| `tools/eqnedit64/src/equation_render.cpp` | nativeレイアウト、display list、文字と構造の描画 |
| `tools/eqnedit64/src/eqnedt64_app.cpp` | Win32 GUI、フォント、クリップボード、CLI、hidden UI tests |
| `tools/eqnedit64/src/palettes.cpp` | パレット分類、face、挿入payloadの単一出典 |
| `tools/eqnedit64/src/*_emitter.cpp` | TeX、MathML、SVG等の出力 |
| `tools/eqnedit64/tests` | native/model/Web/clipboard/visual regression tests |
| `tools/eqnedit64/docs` | GUI、shortcuts、parity、validationの仕様 |
| `packages/eqnedit64` | PyPI metadata、Python API、wheel同梱物 |
| `.github/workflows/eqnedit64.yml` | Eqnedit64専用PR/main/tag CI |
| `.github/workflows/release-eqnedit64-pypi.yml` | O: gate後のwheel/GitHub/PyPI公開 |
| `.agents/skills/release-eqnedit64/scripts/sync_to_o.ps1` | 署名済みEXEとmanifestのO:同期 |

## 15. 退役済み・非目標

次を完成条件や互換目標へ戻さない。

- `.eqn`入力または保存。
- MTEFの継続解読、MTEFを研究室資産として残すこと。
- Equation Editor 3.0のOLE/ActiveX互換。
- Inno Setupまたは別のGUIインストーラー。
- standalone EXEにPythonを要求すること。
- nativeをKaTeX/MathJaxへ置き換えること。
- Web版にnativeと同じ構造GUIを無理に移植してMathJaxの単純さを失うこと。
- 手動貼り付けだけをrelease acceptanceにすること。

## 16. 完成判定

次をすべて満たして初めてリリースを合格とする。

- TeX保存、複数行、GUI/TeX二刀流、Tab移動、Backspace、行列、ショートカットが回帰なし。
- nativeの全メニュー、分類タブ、ボタン、ポップアップセル、ステータスバーが読める。
- `\mathbf`およびvector styleがnative画面でも明確な太字になる。
- PNG/DIBV5/EMFが黒い数式を不透明に表示する。
- PowerPoint/Wordへ見えて編集可能なOffice MathとしてAPI貼り付けできる。
- Google Slides用画像が300 dpi / 24 pt契約を満たす。
- font-hostを壊さず、隔離された32回ライフサイクル試験が通る。
- PR、main、tagのEqnedit64専用CIがgreen。
- `origin/main`、O: manifest、tag、PyPI wheel、GitHub Releaseのsource/version/hashが一致。
- GitHub ReleaseとPyPIから再取得した成果物を外側から検証済み。

「こちらでは起動した」「形式が一つ登録された」「手で貼れた」は補助情報であり、完成判定の
代わりにしない。報告された見た目の不具合は、修正と同時に旧版を落とす自動試験へ変換する。
