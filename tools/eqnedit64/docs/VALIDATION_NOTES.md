# Eqnedit64 検証記録

この文書の2026-08-25/26記録は、公開リポジトリ移行前の内部配布ゲートで得た
履歴である。現在の公開リリース手順は
[`CANONICAL_OPERATION.md`](CANONICAL_OPERATION.md) と `release-eqnedit64` skillを
正とする。`build\accept_release.ps1`は隔離CI/VM用であり、対話中LABでは実行しない。

## 2026-08-29 Eqnedit64 3.0.8公開候補検証

- 利用者の実機確認でnative/Webの見た目は一致したが、数式末尾のカーソル位置が
  18 ptと表示された。保存済み3.0.7 PPTXを調べると、`m:oMath`内の数式runは
  `sz="2400"`だった一方、左寄せ維持用の末尾NBSP runには`sz`がなく、PowerPointの
  既定18 ptへ戻っていた。従来試験の`TextRange.Font.Size`は先頭の数式runだけを返し、
  この不一致を見逃していた。
- native/WebともCF_HTML末尾を
  `<span style="font-size:24pt">&#160;</span>`とし、数式本体だけでなく末尾カーソル
  位置も24 ptにした。実PowerPoint COMで全体24.0 pt、末尾文字24.0 pt、shape left
  0 pt、first ink x=10 px、インライン`m:oMath`、`m:oMathPara`なしを確認した。
- 基準式のnative/Web保存段落はともに6,003文字、SHA-256
  `44C7FB601F147BD44D6BA881528A02CB5B868BBB85C9ADA36E23A544C0AA5FC1`で
  byte-identical。`m:oMath`は5,842文字、SHA-256
  `226979325F9D7F119B20145361635FE230A310DA4A509B9439C507A03AC6433F`、
  PowerPoint描画640x98 PNGはSHA-256
  `29C82D3B2DF3AB5C573F610DBF42F7EEB9637EF7D923CD75E2CD44B034B32070`で、
  いずれもbyte-identicalだった。
- `build/test_external_paste.ps1`と研究室ホームページ集中QAは、全体サイズに加えて
  末尾文字のCOMサイズと保存OOXMLの末尾`sz="2400"` runを必須にした。18 ptの
  カーソルが混じる場合はリリース失敗となる。

## 2026-08-29 Eqnedit64 3.0.7公開検証

- tag / source: `eqnedit64-v3.0.7` /
  `16a739f89f73f53798d823ed9ed9eb0f430e8c55`。
- 3.0.6の合格判定は外側テキスト段落の左揃えだけを見ており、数式本体が中央の
  `m:oMathPara`になった回帰を見逃していた。公開候補では通常コピーから登録MathMLと
  Web版の直接OMMLを除き、両版をinline 24 pt MathML + 末尾NBSPの同期CF_HTMLへ統一した。
- 基準式
  `x+\sum_{n=1}^{m} a^3 \int_{a}^{b} \frac{f(x)}{\sqrt{y}}\, dx^3 + \overline{u} + \underline{v}`
  を両版から非表示の実PowerPointへ貼り付けた。両方とも24.0 pt、shape left 0 pt、
  first ink x=10 px、インライン`m:oMath`、n-ary 2個、分数、根号、上線、下線となり、
  `m:oMathPara`は0個だった。
- PowerPoint保存OOXMLから抽出した`m:oMath`は両版とも5,842文字、SHA-256
  `226979325F9D7F119B20145361635FE230A310DA4A509B9439C507A03AC6433F`でbyte-identical。
  PowerPoint自身が書き出した640x98 PNGもSHA-256
  `29C82D3B2DF3AB5C573F610DBF42F7EEB9637EF7D923CD75E2CD44B034B32070`で
  byte-identicalだった。
- native外部試験はPowerPointに加え、不透明DIBV5、IrfanView、Google Slides
  300 dpi / 24 pt、`--texclip`も合格した。
- PR CI `33250598409`、main CI `33250713645`、tag CI `33250870137`、
  release workflow `33250963981`はすべて成功。PR/mainのPolicy Lint
  `33250598355` / `33250713649`も成功した。
- exact mainからLABで再ビルドしたEXEはProductVersion/FileVersion 3.0.7、
  Authenticode `Valid`、signer `CN=ksugahar`。SHA-256は
  `E3A923C583CBA3B8162D4A711F152A7CBBD970794A72D3B3B34D51A0166C356D`。
  `O:\Eqnedit64.exe`、GitHub Release EXE、PyPIのCPython
  3.10/3.11/3.12/3.13 wheel内EXEはすべてbyte-identicalだった。
- 研究室ホームページには`elemag/index.php`と`elemag/equation-editor.js`だけを
  加算更新した。公開JS、生成JS、Radia正本のSHA-256はすべて
  `FAF9CF9F4A9139AFC416B4C8F0ED4CD9F080272484EEE7A6A809ED753CC5A811`、
  cache keyは`faf9cf9f4a91`。公開URLの2 viewport browser QAと実PowerPoint QAも
  24 pt、左端ink 10 px、n-ary 2個、分数、根号、上線、下線、3,393 ink pixel、
  native OMML/PNGとのbyte-identical判定で合格した。
- GitHub Release:
  `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.7`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.7/`

## 2026-08-29 Eqnedit64 3.0.6公開検証

- tag / source: `eqnedit64-v3.0.6` /
  `a9289bf1dd1651d777d1f8599766d194d854f65c`。
- Web版のPowerPointコピーを同期CF_HTML fragmentへ戻し、条件付きOMMLと
  inline MathML fallbackを同じ24 pt・左揃え段落に載せた。HTML documentの入れ子は
  作らない。
- PR CI `33247565457`、main CI `33247680323`、tag CI `33247910366`、
  release workflow `33247996052`はすべて成功。PR/mainのPolicy Lint
  `33247565431` / `33247680305`も成功した。
- exact mainからLABで再ビルドしたEXEはProductVersion/FileVersion 3.0.6、
  Authenticode `Valid`、signer `CN=ksugahar`。SHA-256は
  `406D29566312E645070FAE9CC925B2A4DEC05B74818AFD37E99050FEB4D2ECAD`。
- native EXEの外部貼り付け試験は24.0 pt、shape left 0 pt、編集可能Office Math、
  n-ary 2個、分数、根号、描画640x190 px / ink 308x174 / fraction run 85で合格した。
  IrfanView、Google Slides 300 dpi、`--texclip`も合格した。旧試験の18--24 pt許容を
  廃止し、24.0 ptと左揃えを必須にした。
- `O:\Eqnedit64.exe`、GitHub Releaseから再取得したEXE、PyPIのCPython
  3.10/3.11/3.12/3.13 wheel内EXEはすべてbyte-identicalで、上記version・署名・hashを
  満たした。CPython 3.12 wheelの新規venvインストール試験も合格した。
- 研究室ホームページには`index.php`と`equation-editor.js`だけを公開した。公開JSの
  SHA-256 `4F43994C6212BA24F70488479D1C0779026AF475821F77BEF138CD93496BC50B`は
  exact main由来の生成物と一致し、cache keyは`4f43994c6212`。
- 公開URLのbrowser QAはdesktop/mobileの2 viewport、実PowerPoint貼り付けは
  24.0 pt、paragraph alignment 1、shape left 0 pt、n-ary 2個、分数、根号、
  上線・下線2個、3,847 ink pixelとして合格した。
- GitHub Release:
  `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.6`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.6/`

## 2026-08-29 Eqnedit64 3.0.5公開検証

- tag / source: `eqnedit64-v3.0.5` /
  `a8b8a5a3d4607cb6c1f5d68421d7bc8b829f0a36`。
- nativeの幾何タブをWeb版と
  同じ11項目・同じ順序へ修正し、独立した「上線と下線」パレットは基本タブに保持した。
- `test_palettes.py`は19パレット・245セルと、native/Webの微分幾何face・TeX・labelの
  完全一致を確認した。`test_edit.py`は271件、`test_tex_fuzz.py`は3,075件、
  Web契約は10件が合格した。
- native通常コピーはCF_HTML・`MathML`・`MathML Presentation`へ同じinline 24 pt
  MathMLを発行する。実PowerPointは24 ptの編集可能Office Math、保存OOXMLでn-ary
  2個・分数・根号、描画640×190 px / ink 308×174 / fraction run 85として合格した。
- Web版は同じfallback MathMLから条件付きOMMLを生成する。研究室ホームページの集中
  QAは2 viewport、実PowerPointのn-ary 2個・分数・根号・上線・下線2個、描画
  ink 2,148 pxを確認した。ブラウザはWindows登録MathML/RTFを公開できないため、
  PowerPointのHTML取り込みサイズは18 pt（EXE登録MathML経路は24 pt）である。
- PR CI `33241671050`、main CI `33241773295`、tag CI `33242009119`、
  release workflow `33242090467`はすべて成功。PR/mainのPolicy Lintも成功した。
- exact mainからLABで再ビルドしたEXEはProductVersion/FileVersion 3.0.5、build stamp
  `a8b8a5a3d`、Authenticode `Valid`、signer `CN=ksugahar`。SHA-256は
  `C2221FB9127E5868E5CFC2D34469012BA1E8384CBFFDD7EDB76B266457E16AED`。
- `O:\Eqnedit64.exe`、GitHub Releaseから再取得したEXE、PyPIのCPython
  3.10/3.11/3.12/3.13 wheel内EXEはすべてbyte-identicalで、上記version・署名・hashを
  満たした。CPython 3.12 wheelの新規venvインストール試験も合格した。
- GitHub Release:
  `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.5`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.5/`

## 2026-08-29 Eqnedit64 3.0.4公開検証

- tag / source: `eqnedit64-v3.0.4` /
  `abb66868da37032e353e44061ac5b4d24cc33d5c`。
- PR CI `33236393427`、main CI `33236509923`、tag CI `33236746356`、
  release workflow `33236832046`はすべて成功。PR/mainのPolicy Lintも成功した。
- exact mainからLABでビルドしたEXEはProductVersion/FileVersion 3.0.4、build stamp
  `abb66868d`、Authenticode `Valid`、signer `CN=ksugahar`。
- `O:\Eqnedit64.exe`、GitHub Releaseから再取得したEXE、PyPIのCPython
  3.10/3.11/3.12/3.13 wheel内EXEはすべてbyte-identical。SHA-256は
  `7BA0547D8779F97285C8B20857F7C395688386FBFF38B524DEABE494F3043121`。
- 4 wheel内EXEはいずれもProductVersion 3.0.4、署名`Valid` / `CN=ksugahar`。
  CPython 3.12 wheelを新規venvへインストールした`verify_installed.py`も合格した。
- GitHub Release:
  `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.4`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.4/`

## 2026-08-29 Eqnedit64 3.0.4公開候補

- 対象: `codex/eqnedit64-math-style-palette`。
- EXE/JS共通の常設数式スタイルパレットと、`\mathnormal`、`\mathrm`、
  `\mathit`、`\mathbf`、`\mathsf`、`\mathtt`、`\mathcal`、`\mathbb`、
  `\mathfrak`、`\bm`、`\boldsymbol`の入力・表示・変換経路を追加した。
- native rendererは埋め込みLatin Modern MathのUnicode数式アルファベットを使い、
  `\bm`のギリシャ文字も太字で描画する。JSは入力TeXの`\bm`を保持し、MathJax境界で
  `\boldsymbol`へ変換する。
- 隔離CI run `33233890433`で、32回font-host lifecycle、model suites、Office、
  IrfanView、Google Slides、texclip、ASanを含むEqnedit64 jobが合格した。
- `\bm{\wedge beta}`の先頭空白が二重正規化で増える回帰を修正し、永続C++単体試験を
  追加した。ローカルのCTest 2件、Web契約8件、JS挿入試験、228パレット項目の
  TeXコンパイル試験は合格した。
- 公開時はexact `origin/main`から3.0.4を再ビルド・署名し、O: manifest更新後にのみ
  `eqnedit64-v3.0.4`を付与する。

## 2026-08-29 Eqnedit64 3.0.3公開検証

- tag / source: `eqnedit64-v3.0.3` /
  `c15e626d28245ee2db5b688f52d7756e000780a4`。
- PR CI `33227677769`、main CI `33227853273`、tag CI `33228191426`、
  release workflow `33228278493`はすべて成功。
- exact mainからLABで再ビルドしたEXEはProductVersion 3.0.3、build stamp
  `c15e626d2`、Authenticode `Valid`、signer `CN=ksugahar`。
- exact release EXEによる非表示外部試験は、PowerPointの24 pt編集可能Office Math
  （export 640×92 px、ink 90×74、fraction run 89）、IrfanView 94×78不透明PNG、
  Google Slides 294×243 px / 300 dpi / 70.56×58.32 pt、`--texclip` PNG/DIBV5が合格。
- `O:\Eqnedit64.exe`、GitHub Releaseから再取得したEXE、PyPI CPython 3.12 wheel内EXEは
  byte-identical。SHA-256は
  `C6D88C59227718FFE8A9BDB435632EAC64F5709615520FBA04E2C0A5C77C8A8F`。
- 再取得したGitHub EXEとwheel内EXEの署名はともに`Valid` / `CN=ksugahar`。
- PyPIはCPython 3.10、3.11、3.12、3.13の`win_amd64` wheel計4個を公開。
- GitHub Release:
  `https://github.com/ksugahar/Radia/releases/tag/eqnedit64-v3.0.3`
- PyPI: `https://pypi.org/project/eqnedit64/3.0.3/`

## 2026-08-29 3.0.3表示回帰のローカル事前検証

- 対象: `codex/eqnedit64-3.0.3-bold-palette`の未公開候補。
- `test_layout.py`: `\vec{\mathbf{E}}`が通常EではなくU+1D404
  MATHEMATICAL BOLD CAPITAL Eをdisplay listへ出すことを含め合格。
- `test_palettes.py`: 11件合格。上下付きfaceの実数値0件、埋め込みフォントに無い
  U+25AF 0件、結合文字0件、行列操作faceが`+R`/`−R`/`+C`/`−C`であることを確認。
- `--visual-scale-test`: 終了コード0。96/120/144/192 dpiの各条件でplain/boldの
  GDI画素差とink増加、および228個すべてのowner-drawセルのcmap ownershipと実inkを確認。
- `test_tex_document.exe`: 10件合格。
- `--self-test`: 終了コード0。
- `--ui-interaction-test`: 終了コード0。
- version resource: ProductVersion/FileVersionとも3.0.3。
- Authenticode: `Valid`、signer `CN=ksugahar`。
- `dist`: `Eqnedit64.exe` 1ファイルだけ。事前ビルドSHA-256は
  `84E835965992791B36A082E5294838DDE43D6F4E875C01126165F2959F17C2F4`。
- このSHAは未コミット候補を`22fc5630e` stampでビルドした事前確認値であり、公開値ではない。
  release commitをmainへ統合後、exact `origin/main`から再ビルド・再署名して置き換える。
- 32回font-host lifecycle、全model suite、Office外部貼り付け、ASan、release acceptanceは
  対話中LABでguardを迂回せず、PR/main/tag CIと最終release gateで実施する。

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
- PR CI run `33215826230`では32回のライフサイクル後も`fontdrvhost.exe`のPID
  7712を維持し、同ゲートを含むEqnedit64 job全体が2分3秒で合格した。

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
