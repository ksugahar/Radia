# Eqnedit product parity policy

## Native rendering acceptance

Every supported decoration must render its own mathematical shape, not merely
nonempty ink. A valid TeX/MathML result does not certify the native canvas.
Keep the decoration pairwise comparisons and complete palette sweep, including
styled selection and script fragments with explicit bases.

Native `Node::Tag` dispatch must contain an explicit case for every enumerator
and no default. `SizeNode` is consumed as list state; reaching scalar dispatch
is a logic error. The unused legacy `FontNode` and `RMNode` are unsupported and
must throw a named error, never produce an empty layout. New tags require an
explicit layout decision under MSVC `/W4 /WX /w14062` in release and CMake CI.
The old generic layout fallback is forbidden. Static coverage and C++ runtime
rejection tests protect separate failure modes.

## Palette key acceptance

A key is a promise about what pressing it produces, so the face is part of the
contract and not decoration. Two rules hold, both checked from the source and
the shipped font rather than from a rendered window.

A symbol key must draw the character its command inserts. `\star` renders as
U+22C6, so a key for it may not show U+2605; `\frown` and `\smile` are U+2322
and U+2323, not the intersection and union signs. A face may be a word only
where the inserted character has no standalone shape, such as the combining
overlay behind `\not`.

Every face character must be in the embedded Eqnedit Math cmap.
`pick_button_font` gathers all faces into one sample and accepts a font only if
it owns the whole sample, so an unavailable character does not blank its own
key — it rejects the math font and redraws every palette in a fallback. That
all-or-nothing gate stays; the sample is what must be kept clean.

The cmap check predicts that choice. `--self-test` observes it: it runs the
production chooser, reads back the physical face, and exits 243 naming the
substitute when the palette is not drawing in Eqnedit Math. Keep both —
the static check says which character is at fault, the runtime check says
whether the shipped binary actually got the font.

These are native-specific implementation checks. Web rendering belongs to
MathJax; shared palette/TeX/Office contracts still require both editions' tests.
Publish the matching Web build on the laboratory homepage, even when the only
Web source change is release identification. Do not claim public completion
from a successful binary upload or from the browser tests alone.

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
| TeXソース文字の可読性 | 物理face・cmap・実inkを検査しstock GUI fontへ退避 | 共有CSS変数によるCJK対応monospace fallback | 共通学習面 |
| `Tab` / `Shift+Tab`で空欄移動 | 構造GUI・TeXソース | TeXソース | 共通学習面 |
| `Enter`の構造的行区切り・`Shift+Enter`のソース改行 | 構造キャンバス・TeXソース | TeXソース | 共通学習面 |
| 挿入テンプレートの整形（`\begin`後・`\end`前・`\\`後で改行し、環境深さで字下げ） | 構造モデルの正規化TeX | 挿入断片へ同じ規則を適用 | 共通学習面 |
| 挿入・行区切りの取り消し | 単一Undoへまとめる | ブラウザのUndo履歴（`Ctrl+Z`）を壊さない | 共通学習面 |
| 論文TeXの意味保全 | native parserで標準環境・コメント・メタ命令を正規化 | MathJaxの標準TeX解釈へ委譲 | `~`を`\sim`にせず、`align`の`&`を失わず、命令名を本文にしない |
| 構造キャンバス編集 | 必須 | 非該当 | native固有 |
| 数式3.0由来ショートカット | 互換層として必須 | 非該当 | native固有 |
| インストール不要のブラウザ利用 | 非該当 | 必須 | Web固有 |

2026-10-07 の UXP-0029 は、PowerPoint貼り付けで積分上下限とスクリプトの保存を
左寄せより優先する。Windows版の単一式（cases・行列など式内の表を含む）は
登録 `MathML` / `MathML Presentation` を発行し、中央寄せの編集可能な
`m:oMathPara` を許容する。CF_HTML・18 pt末尾NBSPも保持する。
独立した複数行はCF_HTMLで各行を渡し、左端と上下限を実PowerPoint試験で検査する。
Web版はブラウザーのCF_HTML経路を使う。現在のPowerPointでのWeb・複数行の上下限は
隔離Office試験で再確認するまで未検証とし、native/WebのOOXML・PNG一致は保証しない。
画像への退化・空のMathZone・上下限の欠落は合格としない。明示 `&` と保存TeXは保持する。

論文からのTeX入力では、両版とも裸の`~`を非改行空白、U+223Cを`\sim`として区別し、
`align` / `align*`の`&`列を保持する。`%`コメント、`\label`、`\nonumber`、`\notag`、
`\tag`を可視数式へ変えない。nativeは未知control wordを無視して後続引数を編集可能に残し、
Web/MathJaxは未対応命令を変換エラーとして通知してよいが、どちらも命令名を通常文字として
数式へ捏造してはならない。この差は構造編集を継続できるnativeと、完全TeX parserを持つ
MathJaxの境界による。
native parserの固定点、optional引数の停止条件、Undo用再parseはnative構造編集固有の
実装責任であり、MathJaxへTeX解釈を委譲するWeb版へC++ parserを移植しない。逆に、
Web版のソース文字はブラウザのfont fallbackへ委ねるため、Win32のGDI font検査を
JavaScriptへ再実装しない。ソース欄と直前挿入表示は一つのCJK対応font stackを共有し、
特定のローカルfont一つだけを前提にしない。
TeXソース欄の`Enter`は両版とも数式の行区切り`\\`であり、単なる空白改行ではない。
新しい行は直前の行の整列列を保ち（先頭に`&`）、`\frac`、上下付き、`\text`、入れ子環境の
途中では群を分割せず、現在の行の後ろへ空の行を開く。行環境がまだ無い入力は`aligned`で
包む。ソースだけを改行する場合は両版とも`Shift+Enter`を使う。Web版はさらに、日本語入力の
確定`Enter`（`isComposing` / key code 229）をIMEへ渡し、数式の行区切りにしない。
挿入テンプレートのセルは`{}`とし、`Tab`が全セルへ到達できることを両版の条件とする。
Web版の挿入・行区切りはブラウザの編集コマンド経由で適用し、`textarea.value`の代入で
Undo履歴を捨てない。ブラウザのUndo履歴はnativeのUndoスタックに対応する共通学習面である。

両版とも`\text{hello world}`と`\operatorname{arg max}`の単語間空白を表示・保存し、
escaped percentを含む文章でも空白を脱落させない。nativeは構造木へ取り込む際に連続する
source空白を一つの単語間空白へ正規化するが、`~`、`\ `、`\,`の明示空白は区別して保持する。

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
- native変換器の公開境界はPythonを介さない
  `Eqnedit64.exe INPUT OUTPUT`である。CLI dispatchをpybind11へ複製しない。Python
  packageのlauncherや補助関数は必要に応じて同じEXEをsubprocess実行する薄い利用者であり、
  CLIの正本ではない。Web版にnative CLIを移植しない。
- 仕様、実装、自動試験を同じ変更で更新する。引き継ぎメモだけを規範にしない。
# Palette semantic acceptance

Every offered key must pass the independent [palette intent contract](PALETTE_INTENT.md)
and its checked expectations in `palette_intent.json`. Catalogue completeness,
distinct drawings and TeX round trips alone do not establish the intended meaning.
Native and Web insertions are checked separately; an unreviewed key fails closed.

## Web first-copy readiness

The Web Office-copy button stays disabled with a preparation explanation until
the primary CHTML typesetter, followed by the MathML converter, has initialized
bold-symbol and cancellation macros. Verify synchronous conversion before
enabling copy. Initialization failure leaves copy disabled with a reload
explanation; it must not masquerade as successful readiness. Defer the initial
expression render until readiness too. Palette preview rendering follows this
preparation but must not delay enabling Office copy.

Keep clipboard writing synchronous inside the user gesture. This readiness
change does not alter the inline MathML payload, 18 pt convention, or native
clipboard path. The cold-browser test delays the cancellation package and holds
all palette previews, checks the disabled state, then captures the actual copy
event before any preview completes. The first payload must contain bold-italic
and cancellation markup, without MathJax errors. It does not modify the user's
system clipboard. Website hand testing remains required before publication.

## Native named colours (3.1.2 candidate, UXP-0030)

Both editions accept named `\color{red}{x}` and `\textcolor{red}{x}`.
`\color` switches the colour of the remainder of the enclosing group, as in
MathJax: `\color{red}{x}+y` colours both x and y; `{\color{red}x}+y`
colours only x. `\textcolor` colours only its braced body. Nested colours
restore the enclosing colour after the inner body.

Native supports 148 CSS named colours and 68 case-sensitive MathJax dvips
names (for example `green` = #008000, `Green` = #00A64F), and `[named]`.
Numeric models (`rgb`, `RGB`, `gray`, `HTML`), colour definitions, CSS colour
expressions and transparent/currentColor are not supported in native 3.1.2;
an unsupported name/model reports an error and leaves the valid model intact.
The Web MathJax package has a wider colour-expression/model vocabulary; parity
is claimed only for the supported named-colour inputs.

Colour remains in the editable group, selection, Undo, normalized TeX save /
reopen, GDI canvas / bitmap / outline EMF and SVG glyphs and rules. Native
CF_HTML and registered MathML emit `mathcolor`; Office TeX retains
`\textcolor`. PowerPoint's actual colour retention on this new candidate
requires a new user hand test; earlier 02f74f70c H5/H6 results are historical,
not evidence for the native colour change. The isolated external-paste suite
also includes a red matrix with integral bounds, but has not been run locally.

Other unknown control words retain the historical permissive parser fallback
and can still disappear silently. Proposed follow-up: strict diagnostics on
load/paste/CLI with recoverable source editing, plus a documented compatibility
allowlist; this broader parser change is not included in UXP-0030.

### PowerPoint colour retention: current evidence

The current user hand test of Web `\color{red}{x}` copies successfully,
but PowerPoint imports it without the red colour. This is a known Web paste
limitation. Native emits colour in MathML and Office TeX, but colour retention
by PowerPoint on the new native colour candidate has not been measured.
Neither edition promises that PowerPoint retains colour. Native canvas,
SVG, EMF/bitmap exports and saved TeX carry the colour independently of
PowerPoint import. Publishing `mathcolor` is not evidence that Office keeps it.
Token-level MathML, mstyle, CF_HTML CSS and Office TeX import routes remain
pending a scoped comparison; do not infer that every Office route loses colour
from the observed Web result alone.
