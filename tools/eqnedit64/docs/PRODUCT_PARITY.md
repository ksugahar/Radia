# Eqnedit product parity policy

Eqnedit64.exe と Web/JS 数式エディタは、Radia の `tools/eqnedit64` で一緒に
保守する同一製品系列である。研究室ホームページは Web 版の公開先であり、正本では
ない。TeX を正本とし、次の優先順位を守る。

## 製品価値の優先順位

1. **MathML 経由の Office ネイティブ数式貼り付け**
   - Word / PowerPoint へ、目に見え、編集できる Office Math として貼り付くこと。
   - MathML・MathZone・XMLの存在だけでは合格にしない。Office自身の描画結果が
     非空で、分数・根号などの代表構造を保持することをAPI試験で確認する。
2. **GUI と TeX ソースの二刀流入力**
   - パレット／構造GUIから入り、同時に標準TeXの綴りを確認できること。
   - 挿入直後のTeXだけを強調し、説明文の暗記を要求しない。
   - GUIとTeXソースのどちらでも `Tab` / `Shift+Tab` により次／前の空欄へ
     移動でき、同じ手癖で式を埋められること。
3. **Eqnedit64 の高速構造入力**
   - Microsoft 数式3.0（Eqnedit32）から回収した操作系列を互換層として維持する。
   - 互換キーを削除・変更しない。追加ショートカットは拡張層として区別し、
     実キー経路をバックグラウンド試験する。

## 機能区分

| 機能 | Eqnedit64.exe | Web/JS | 区分 |
|---|---|---|---|
| TeXを唯一の正本とする | 必須 | 必須 | 共通中核 |
| MathML経由で編集可能なOffice数式へ貼り付け | inline 18 pt MathML + 18 pt NBSPのCF_HTML | 同じinline 18 pt MathML + 18 pt NBSPのCF_HTML | 共通中核 |
| 画像を混在させないOffice専用コピー | 必須 | 必須 | 共通中核 |
| パレットとTeXソースの対応 | 必須 | 必須 | 共通学習面 |
| 常設書体 `\mathrm` / `\mathit` / `\mathbf` | 選択変更・継続入力 | 選択を包む・空欄挿入 | 共通学習面 |
| 追加書体 `\mathsf` / `\mathtt` / `\mathcal` / `\mathbb` / `\mathfrak` / `\bm` | 装飾パレット・TeX・保存・MathML | 装飾パレット・MathJax | `\boldsymbol` は `\bm` の入力別名 |
| 挿入直後のTeX強調 | ソース内の非アクティブ選択 | キャレットを動かさない強調表示 | 共通学習面 |
| `Tab` / `Shift+Tab`で空欄移動 | 構造GUI・TeXソース | TeXソース | 共通学習面 |
| 構造キャンバス編集 | 必須 | 非該当 | native固有 |
| 数式3.0由来ショートカット | 互換層として必須 | 非該当 | native固有 |
| インストール不要のブラウザ利用 | 非該当 | 必須 | Web固有 |

両製品のMathMLはinline・18 ptとし、18 ptを明示した末尾NBSPを加えた単一CF_HTML
断片だけをOfficeへ渡す。PowerPointの通常Ctrl+Vでは編集可能な数式、末尾文字、
次の挿入点を18 pt・左寄せでそろえる。左寄せを24 pt固定より優先する。
登録MathMLはPowerPointが中央寄せ`m:oMathPara`として優先するため通常コピーへ
載せず、Web版だけの直接OMMLも使わない。生成器（native構造木とMathJax）の空白や
補助要素は異なってよいが、総和、積分、分数、根号、上線、下線を正規化し、保存OOXMLの
インライン`m:oMath`とPowerPoint描画を一致させる。実機基準式ではOMML部分とPNGが
一致することを確認する。
複数行`aligned`はnative/Webとも各行を1つの`mtd`へ正規化する。`&`なしの行頭と
各行頭にMathML `maligngroup`、`&`なしの行頭と明示`&`の位置に`malignmark`を置き、
Officeの`m:eqArr`変換後も長短行の
左端または指定整列点を一致させる。MathML `columnalign="left"`だけではOfficeが
無視するため合格条件にならない。
空のテキストボックス、置換文字だけの描画、画像への退化は不合格とする。

PowerPointの通常貼り付け受入試験は、画面外の一時プレゼンテーションに対して
`Application.CommandBars.ExecuteMso("Paste")`を実行する。`Shapes.Paste()`は
異なるオブジェクトモデル経路で24 ptを保持するため、Ctrl+Vの代用にしない。

## 変更規則

- 共通中核・共通学習面を変更する作業は、同じ作業内で両実装を確認し、適用可能なら
  両方を変更する。片方へ適用しない場合は、この表へプラットフォーム上の理由を記す。
- 両実装が同じコードである必要はない。入力キャレット、ブラウザ権限、Windows
  クリップボード形式などを壊さず、同じ利用者結果を満たす実装を選ぶ。
- native変更では `build/test_background.ps1` と該当する外部貼り付け試験を、Web変更では
  `tests/test_web_contract.py` と研究室ホームページ側のブラウザQAを実行する。Office
  貼り付け変更は両方とも実Office APIによる可視描画試験を省略しない。
- Web版の正本は `web/equation-editor.js` と
  `web/equation-editor.fragment.html`。ホームページへ置いたコピーを直接改変しない。
- MTEFと`.eqn`を互換形式として復活させない。両実装とも入出力の正本はTeXとする。
- 仕様、実装、自動試験を同じ変更で更新する。引き継ぎメモだけを規範にしない。
