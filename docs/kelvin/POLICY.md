# ポリシー / Policy

## 動作環境 / Requirements

以下のバージョンのNGSolveを使用してください：

**NGSolve 6.2.2601 以降**

```bash
pip install ngsolve>=6.2.2601
```

注意: バージョン 6.2.2406~6.2.2501 には周期境界条件の後退バグがありましたが、6.2.2601 で修正されています。

---

Use the following version of NGSolve:

**NGSolve 6.2.2601 or later**

```bash
pip install ngsolve>=6.2.2601
```

Note: Versions 6.2.2406~6.2.2501 had a Periodic BC regression bug, which is fixed in 6.2.2601+.

## コーディングスタイル / Coding Style

インデントにはスペース4文字ではなく、タブ文字を使用してください。

---

Use tab characters for indentation, not 4 spaces.

## 論文のBibTeX管理 / BibTeX Management for Papers

- `references.bib` ファイルには、将来使用する可能性のある参考文献も含めておきます
- 本文中で `\cite{}` していなくても、`references.bib` からは削除しないでください
- BibTeX エントリは必要に応じて追加し、不要な場合でも保持します

---

- The `references.bib` file may contain references that could be used in the future
- Do NOT delete entries from `references.bib` even if they are not cited in the main text with `\cite{}`
- BibTeX entries should be added as needed and retained even if not currently used

## Adaptive Mesh スクリプトの出力 / Adaptive Mesh Script Output

- 適応メッシュ refinement スクリプトは、各イテレーションで PNG 画像を生成すること
- 2x2 レイアウト: 内部メッシュ、外部（Kelvin変換）メッシュ、収束曲線、誤差履歴
- ファイル名: `*_iter_01.png`, `*_iter_02.png`, ... の形式

---

- Adaptive mesh refinement scripts should generate PNG images at each iteration
- 2x2 layout: inner mesh, outer (Kelvin-transformed) mesh, convergence curve, error history
- Filename format: `*_iter_01.png`, `*_iter_02.png`, ...

## 図のフォントサイズ / Figure Font Size

- すべてのPNG画像のフォントサイズは、横幅8cmで表示したときに10ポイントになるように設定すること

---

- All PNG images should have font sizes that appear as 10pt when displayed at 8cm width

## ドキュメント管理 / Document Management

- `README.md` ファイルは更新するたびに `md2html` で HTML ファイルに変換してください。
- フォーマット: `md2html input.md > output.html`

---

- `README.md` files must be converted to HTML using `md2html` every time they are updated.
- Format: `md2html input.md > output.html`
