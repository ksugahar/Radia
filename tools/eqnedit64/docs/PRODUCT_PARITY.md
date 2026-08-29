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
| MathML経由で編集可能なOffice数式へ貼り付け | inline MathMLを含むHTML + 同一登録MathML | 同じfallback MathML + 条件付きOMMLを含むHTML | 共通中核 |
| 画像を混在させないOffice専用コピー | 必須 | 必須 | 共通中核 |
| パレットとTeXソースの対応 | 必須 | 必須 | 共通学習面 |
| 常設書体 `\mathrm` / `\mathit` / `\mathbf` | 選択変更・継続入力 | 選択を包む・空欄挿入 | 共通学習面 |
| 追加書体 `\mathsf` / `\mathtt` / `\mathcal` / `\mathbb` / `\mathfrak` / `\bm` | 装飾パレット・TeX・保存・MathML | 装飾パレット・MathJax | `\boldsymbol` は `\bm` の入力別名 |
| 挿入直後のTeX強調 | ソース内の非アクティブ選択 | キャレットを動かさない強調表示 | 共通学習面 |
| `Tab` / `Shift+Tab`で空欄移動 | 構造GUI・TeXソース | TeXソース | 共通学習面 |
| 構造キャンバス編集 | 必須 | 非該当 | native固有 |
| 数式3.0由来ショートカット | 互換層として必須 | 非該当 | native固有 |
| インストール不要のブラウザ利用 | 非該当 | 必須 | Web固有 |

両製品のfallback MathMLはinline・24 ptとし、総和は上下限、通常の積分はTeX同様の
右側上下限という共通構造契約を持つ。Web版はブラウザ制約を越えるため、同じMathML
から条件付きOMMLも生成してPowerPointへ渡す。生成器（native構造木とMathJax）が
異なるためXMLの補助要素や直列化バイト列は同一でなくてよいが、保存OOXML上の分数、
根号、n-ary、上線、下線と目視結果は一致させる。ブラウザはWindows登録MathMLを
発行できないためPowerPointのWeb HTML経路は18 pt、EXE登録MathML経路は24 ptである。
空のテキストボックス、置換文字だけの描画、画像への退化は不合格とする。

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
