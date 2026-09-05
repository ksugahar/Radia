/* SPDX-License-Identifier: BSD-2-Clause
 * Canonical source: ksugahar/Radia tools/eqnedit64/web/equation-editor.js
 *
 * 数式エディタ（Web版） --- 3Dで学ぶ電磁気学 トップページ
 *
 * 研究室のネイティブ数式エディタ（Eqnedit64）の操作思想を、ページに載る
 * 最小構成へ写したもの。TeXソースを唯一の正とし、描画はサイト共通の
 * MathJax に任せる。パレットは「マウスで全記号に届く」規則の縮約版。
 *
 * 操作:
 *   - パレットのボタンでテンプレート・記号を挿入（Ctrl+Z で取り消せる）
 *   - 範囲選択して分数などを押すと、選択が最初の空欄に包まれる
 *   - Tab / Shift+Tab で空欄 {} の間を移動
 *   - Enter は数式の行区切り \\（aligned / matrix / cases の行を増やす）、
 *     Shift+Enter はソースだけの改行
 *   - 「$$付きでコピー」「equation付きでコピー」で持ち出し
 */
(function () {
  "use strict";

  /* デプロイごとに上げる。ボタン行の右端に出て、開きっぱなしのタブが
   * 古い版を動かし続けていないかを一目で判別できる（.exe の
   * タイトルバー・ビルドスタンプと同じ教訓）。 */
  var BUILD = "3.0.15 (2026-09-05b)";

  var PALETTES = [
    {
      label: "構造",
      items: [
        ["a/b", "\\frac{}{}", "分数"],
        ["dfrac", "\\dfrac{}{}", "常に大きい分数（入れ子・文中でも縮まない）"],
        ["√", "\\sqrt{}", "平方根"],
        ["ⁿ√", "\\sqrt[]{}", "n乗根"],
        ["x²", "^{}", "上付き"],
        ["xᵢ", "_{}", "下付き"],
        ["Σ", "\\sum_{}^{}", "総和"],
        ["∏", "\\prod_{}^{}", "総乗"],
        ["∫", "\\int_{}^{}", "積分"],
        ["∬", "\\iint_{} ", "二重積分（面積分）"],
        ["∭", "\\iiint_{} ", "三重積分（体積積分）"],
        ["∮", "\\oint_{}", "周回積分"],
        ["lim", "\\lim_{ \\to }", "極限"],
        ["d/dx", "\\frac{\\mathrm{d} {}}{\\mathrm{d} {}}", "常微分"],
        ["∂/∂x", "\\frac{\\partial {}}{\\partial {}}", "偏微分"],
        /* 行列と場合分けの各セルは空欄 {} にする。Tab で順に埋められ、
         * 挿入時にネイティブ版と同じ改行・字下げで整形して表示する。 */
        ["行列", "\\begin{pmatrix} {} & {} \\\\ {} & {} \\end{pmatrix}", "2×2行列（丸括弧）"],
        ["∷", "\\begin{matrix} {} & {} \\\\ {} & {} \\end{matrix}", "行列（括弧なし）"],
        ["[∷]", "\\begin{bmatrix} {} & {} \\\\ {} & {} \\end{bmatrix}", "行列（角括弧）"],
        ["|∷|", "\\begin{vmatrix} {} & {} \\\\ {} & {} \\end{vmatrix}", "行列式"],
        ["場合", "\\begin{cases} {} & {} \\\\ {} & {} \\end{cases}", "場合分け"],
        ["整列", "\\begin{aligned} {} &= {} \\\\ {} &= {} \\end{aligned}", "複数行（Enter でも行が増える）"],
        ["↑□", "\\overset{}{}", "上側の注記（記号の上に小さく）"],
        ["↓□", "\\underset{}{}", "下側の注記（記号の下に小さく）"]
      ]
    },
    {
      /* すべて \left…\right 対。中央の {} が空欄になり、カーソルは
       * 括弧の内側に着地する（選択して押すと選択が包まれる）。 */
      label: "括弧",
      items: [
        ["( )", "\\left( {} \\right)", "丸括弧"],
        ["[ ]", "\\left[ {} \\right]", "角括弧"],
        ["{ }", "\\left\\{ {} \\right\\}", "波括弧"],
        ["| |", "\\left| {} \\right|", "絶対値"],
        ["‖ ‖", "\\left\\| {} \\right\\|", "ノルム"],
        ["⟨ ⟩", "\\left\\langle {} \\right\\rangle", "山括弧（内積）"],
        ["⌊ ⌋", "\\left\\lfloor {} \\right\\rfloor", "床関数"],
        ["⌈ ⌉", "\\left\\lceil {} \\right\\rceil", "天井関数"],
        ["⟨ | ⟩", "\\left\\langle {} \\middle| {} \\right\\rangle", "ブラケット（\\middle| で仕切る）"]
      ]
    },
    {
      label: "装飾",
      items: [
        ["x̂", "\\hat{}", "ハット"],
        ["x̄", "\\bar{}", "バー"],
        ["x⃗", "\\vec{}", "ベクトル"],
        ["ẋ", "\\dot{}", "ドット"],
        ["ẍ", "\\ddot{}", "二重ドット"],
        ["x⃛", "\\dddot{}", "三重ドット"],
        ["ã", "\\tilde{}", "チルダ"],
        ["x′", "'", "プライム"],
        ["x̸", "\\cancel{}", "打ち消し線"],
        ["¯□", "\\overline{}", "上線（伸縮）"],
        ["_□", "\\underline{}", "下線（伸縮）"],
        ["⏞□", "\\overbrace{}", "上の水平波括弧"],
        ["⏟□", "\\underbrace{}", "下の水平波括弧"],
        ["□→", "\\overrightarrow{}", "上の右向き伸縮矢印"],
        ["←□", "\\overleftarrow{}", "上の左向き伸縮矢印"],
        ["□↔", "\\overleftrightarrow{}", "上の両向き伸縮矢印"],
        ["ℝ", "\\mathbb{}", "黒板太字"],
        ["ℰ", "\\mathcal{}", "カリグラフィー体（起電力ℰなど）"],
        ["sf", "\\mathsf{}", "サンセリフ"],
        ["tt", "\\mathtt{}", "等幅"],
        ["fr", "\\mathfrak{}", "フラクトゥール"],
        ["Bα", "\\bm{}", "記号太字"],
        ["math", "\\mathnormal{}", "標準の数式書体へ戻す"],
        ["abc", "\\text{}", "テキスト（空白が使える）"]
      ]
    },
    {
      /* 関数名は立体で組む。rot/div/grad は和書の流儀（curl ではなく）。 */
      label: "関数",
      items: [
        ["sin", "\\sin ", ""], ["cos", "\\cos ", ""], ["tan", "\\tan ", ""],
        ["log", "\\log ", ""], ["ln", "\\ln ", ""], ["exp", "\\exp ", ""],
        ["rot", "\\operatorname{rot} ", "回転（和書流儀）"],
        ["curl", "\\operatorname{curl} ", "回転（洋書流儀）"],
        ["div", "\\operatorname{div} ", "発散"],
        ["grad", "\\operatorname{grad} ", "勾配"],
        ["arg min", "\\operatorname*{arg\\,min}_{} ", "最小値を与える変数"],
        ["arg max", "\\operatorname*{arg\\,max}_{} ", "最大値を与える変数"],
        ["∇⋅", "\\nabla\\cdot ", "発散（∇表記）"],
        ["∇×", "\\nabla\\times ", "回転（∇表記）"],
        ["Re", "\\operatorname{Re} ", "実部"],
        ["Im", "\\operatorname{Im} ", "虚部"]
      ]
    },
    {
      label: "関係",
      items: [
        ["±", "\\pm ", ""], ["×", "\\times ", ""], ["⋅", "\\cdot ", ""],
        ["∝", "\\propto ", "比例"], ["≈", "\\approx ", ""],
        ["≃", "\\simeq ", "ほぼ等しい"], ["∼", "\\sim ", "同程度"],
        ["≡", "\\equiv ", ""], ["≤", "\\leq ", ""], ["≥", "\\geq ", ""],
        ["≠", "\\neq ", ""], ["≪", "\\ll ", "十分小さい"],
        ["≫", "\\gg ", "十分大きい"], ["⊥", "\\perp ", "垂直"],
        ["∥", "\\parallel ", "平行"], ["∠", "\\angle ", "角・偏角（フェーザ）"],
        ["∓", "\\mp ", "マイナスプラス"], ["÷", "\\div ", "除算"],
        ["∗", "\\ast ", "アスタリスク"], ["•", "\\bullet ", "黒丸"],
        ["⊖", "\\ominus ", "丸マイナス"], ["⊙", "\\odot ", "丸ドット"],
        ["⋄", "\\diamond ", "ダイヤ演算子"], ["▽", "\\bigtriangledown ", "下向き三角"],
        ["∖", "\\setminus ", "差集合"], ["\\", "\\backslash ", "バックスラッシュ"],
        ["†", "\\dagger ", "ダガー（随伴）"], ["‡", "\\ddagger ", "二重ダガー"],
        ["≺", "\\prec ", "順序が先"], ["≻", "\\succ ", "順序が後"],
        ["≅", "\\cong ", "合同"], ["∣", "\\mid ", "割り切る・条件付き"],
        ["‖", "\\Vert ", "二重縦線"], ["⊢", "\\vdash ", "導出される"],
        ["⊨", "\\models ", "充足する"], ["◁", "\\triangleleft ", "左向き三角"],
        ["▷", "\\triangleright ", "右向き三角"], ["⊑", "\\sqsubseteq ", "角付き部分集合"],
        ["⊒", "\\sqsupseteq ", "角付き上位集合"], ["⌢", "\\frown ", "フラウン"],
        ["⌣", "\\smile ", "スマイル"]
      ]
    },
    {
      label: "矢印・集合",
      items: [
        ["→", "\\to ", ""], ["⇒", "\\Rightarrow ", ""],
        ["⇔", "\\Leftrightarrow ", "同値"], ["↦", "\\mapsto ", "写像"],
        ["∈", "\\in ", ""], ["∉", "\\notin ", ""], ["⊂", "\\subset ", ""],
        ["∪", "\\cup ", ""], ["∩", "\\cap ", ""], ["∀", "\\forall ", ""],
        ["∃", "\\exists ", ""], ["∅", "\\emptyset ", ""],
        ["∴", "\\therefore ", "ゆえに"], ["∵", "\\because ", "なぜならば"],
        ["←", "\\leftarrow ", ""], ["↔", "\\leftrightarrow ", ""],
        ["↑", "\\uparrow ", ""], ["↓", "\\downarrow ", ""], ["↕", "\\updownarrow ", ""],
        ["⇐", "\\Leftarrow ", ""], ["⇑", "\\Uparrow ", ""],
        ["⟹", "\\Longrightarrow ", "長い⇒"], ["⟸", "\\Longleftarrow ", "長い⇐"],
        ["⟺", "\\Longleftrightarrow ", "長い⇔"],
        ["↪", "\\hookrightarrow ", "右フック矢印"], ["↩", "\\hookleftarrow ", "左フック矢印"],
        ["↗", "\\nearrow ", ""], ["↘", "\\searrow ", ""], ["↙", "\\swarrow ", ""], ["↖", "\\nwarrow ", ""],
        ["⇀", "\\rightharpoonup ", "右上ハープーン"], ["⇁", "\\rightharpoondown ", "右下ハープーン"],
        ["↼", "\\leftharpoonup ", "左上ハープーン"], ["↽", "\\leftharpoondown ", "左下ハープーン"],
        ["∋", "\\ni ", "要素として含む"], ["⊃", "\\supset ", ""],
        ["⊆", "\\subseteq ", ""], ["⊇", "\\supseteq ", ""],
        ["⊄", "\\not\\subset ", "部分集合でない"],
        ["⊔", "\\sqcup ", "角付き和"], ["⊓", "\\sqcap ", "角付き共通部分"],
        ["⋃", "\\bigcup_{}^{}", "大きい和集合"], ["⋂", "\\bigcap_{}^{}", "大きい共通部分"],
        ["∐", "\\coprod_{}^{}", "余積"],
        ["¬", "\\neg ", "否定"], ["∨", "\\vee ", "論理和"],
        ["∄", "\\nexists ", "存在しない"], ["⊤", "\\top ", "真"]
      ]
    },
    {
      /* 第9章（微分形式・Hodge star・pullback・幾何代数）向け。
       * すべて標準 LaTeX + amsmath/amssymb の範囲で組める。 */
      label: "微分幾何",
      items: [
        ["∧", "\\wedge ", "ウェッジ積"],
        ["★", "\\star ", "Hodge star作用素"],
        ["dω", "\\mathrm{d} ", "外微分（立体のd）"],
        ["ι", "\\iota_{} ", "内部積（縮約）ι_X"],
        ["ℒ", "\\mathcal{L}_{} ", "Lie微分"],
        ["f^*", "^{*} ", "引き戻し（pullback）"],
        ["f_*", "_{*} ", "押し出し（pushforward）"],
        ["♭", "^{\\flat} ", "フラット（添字を下げる）"],
        ["♯", "^{\\sharp} ", "シャープ（添字を上げる）"],
        ["⊗", "\\otimes ", "テンソル積"],
        ["⊕", "\\oplus ", "直和"]
      ]
    },
    {
      label: "その他",
      items: [
        ["∂", "\\partial ", ""], ["∇", "\\nabla ", ""], ["∞", "\\infty ", ""],
        ["ℏ", "\\hbar ", ""], ["°", "^{\\circ} ", "度"],
        ["⋯", "\\cdots ", "横の点"], ["⋮", "\\vdots ", "縦の点"],
        ["⋱", "\\ddots ", "斜めの点"], ["…", "\\ldots ", "下寄りの点"],
        ["ℜ", "\\Re ", "実部（フラクトゥール）"], ["ℑ", "\\Im ", "虚部（フラクトゥール）"],
        ["ℓ", "\\ell ", "筆記体の l"], ["℧", "\\mho ", "モー（コンダクタンス）"],
        ["′", "\\prime ", "プライム記号"],
        ["⌊", "\\lfloor ", ""], ["⌋", "\\rfloor ", ""],
        ["⌈", "\\lceil ", ""], ["⌉", "\\rceil ", ""]
      ]
    },
    {
      /* TeX の空白命令。文字幅は TeX 標準（\, = 1/6 em など）。 */
      label: "空白",
      items: [
        ["\\,", "\\, ", "細い空き（1/6 em）"],
        ["\\:", "\\: ", "中位の空き（2/9 em）"],
        ["\\;", "\\; ", "広い空き（5/18 em）"],
        ["\\!", "\\! ", "負の空き（詰める）"],
        ["quad", "\\quad ", "1 em の空き"],
        ["qquad", "\\qquad ", "2 em の空き"],
        ["~", "~", "改行しない空白"],
        ["\\␣", "\\ ", "通常幅の空白（\\ の後に空白）"]
      ]
    },
    {
      label: "ギリシャ",
      items: [
        ["α", "\\alpha ", ""], ["β", "\\beta ", ""], ["γ", "\\gamma ", ""],
        ["δ", "\\delta ", ""], ["ε", "\\varepsilon ", ""], ["ζ", "\\zeta ", ""],
        ["η", "\\eta ", ""], ["θ", "\\theta ", ""], ["κ", "\\kappa ", ""],
        ["λ", "\\lambda ", ""], ["μ", "\\mu ", ""], ["ν", "\\nu ", ""],
        ["ξ", "\\xi ", ""], ["π", "\\pi ", ""], ["ρ", "\\rho ", ""],
        ["σ", "\\sigma ", ""], ["τ", "\\tau ", ""],
        ["ϕ", "\\phi ", "ファイ（スカラーポテンシャル）"],
        ["φ", "\\varphi ", ""], ["χ", "\\chi ", ""], ["ψ", "\\psi ", ""],
        ["ω", "\\omega ", ""],
        ["ϵ", "\\epsilon ", "イプシロン（\\varepsilon とは別字形）"],
        ["ϑ", "\\vartheta ", "シータ（異体）"], ["ϰ", "\\varkappa ", "カッパ（異体）"],
        ["ϖ", "\\varpi ", "パイ（異体）"], ["ϱ", "\\varrho ", "ロー（異体）"],
        ["ς", "\\varsigma ", "シグマ（語末形）"], ["υ", "\\upsilon ", "ウプシロン"]
      ]
    },
    {
      label: "ギリシャ大",
      items: [
        ["Γ", "\\Gamma ", ""], ["Δ", "\\Delta ", ""], ["Θ", "\\Theta ", ""],
        ["Λ", "\\Lambda ", ""], ["Φ", "\\Phi ", "磁束Φなど"],
        ["Ψ", "\\Psi ", ""], ["Ω", "\\Omega ", ""],
        ["Ξ", "\\Xi ", ""], ["Π", "\\Pi ", ""], ["Σ", "\\Sigma ", ""],
        ["Υ", "\\Upsilon ", ""]
      ]
    }
  ];

  /* 行を用途別のタブに束ねる。ラベルは従来のまま残すので、タブを
   * 開いた後は以前と同じ位置関係で記号を探せる。 */
  var PALETTE_TABS = [
    { id: "basic", label: "基本", palettes: ["構造", "括弧", "装飾"] },
    { id: "analysis", label: "解析", palettes: ["関数", "関係"] },
    { id: "sets", label: "集合・記号", palettes: ["矢印・集合", "その他", "空白"] },
    { id: "geometry", label: "幾何", palettes: ["微分幾何"] },
    { id: "greek", label: "ギリシャ", palettes: ["ギリシャ", "ギリシャ大"] }
  ];

  /* Category-independent math alphabets.  These stay beside the tabs so a
   * user never has to remember which subject palette owns a writing style. */
  var MATH_ALPHABETS = [
    { label: "R x", css: "roman", snippet: "\\mathrm{}", name: "立体" },
    { label: "I x", css: "italic", snippet: "\\mathit{}", name: "変数（斜体）" },
    { label: "B x", css: "vector", snippet: "\\mathbf{}", name: "ベクトル" }
  ];


  var CSS = [
    ".eqed { --eqed-source-font: ui-monospace, 'Cascadia Mono', 'Yu Gothic UI', Meiryo, monospace; border: 1px solid #d9d9d6; border-radius: 10px; padding: 14px 16px; background: #fcfcfb; }",
    ".eqed-palette-head { display: flex; flex-wrap: wrap; gap: 6px 12px; align-items: flex-start; padding: 0 0 7px; }",
    ".eqed-tabs { display: flex; flex: 1 1 auto; gap: 4px; max-width: 100%; overflow-x: auto; scrollbar-width: thin; }",
    ".eqed-tab { flex: 0 0 auto; font: inherit; font-size: 0.86rem; padding: 6px 12px; border: 1px solid #cfcfcb; border-radius: 7px 7px 3px 3px; background: #f7f7f5; color: inherit; cursor: pointer; }",
    ".eqed-tab:hover { background: #eef3f8; border-color: #9db8d2; }",
    ".eqed-tab[aria-selected='true'] { background: #fff; border-color: #6f98bd; box-shadow: inset 0 -2px #3977ad; font-weight: 600; }",
    ".eqed-tab:focus-visible { outline: 3px solid rgba(57,119,173,0.3); outline-offset: 1px; }",
    ".eqed-stylebar { display: flex; flex: 0 0 auto; gap: 4px; margin-left: auto; }",
    ".eqed-style-key { min-width: 3.3em; font-family: 'Times New Roman', serif; font-size: 0.95rem; line-height: 1; padding: 6px 8px; border: 1px solid #b9c6d1; border-radius: 6px; background: #fff; color: inherit; cursor: pointer; }",
    ".eqed-style-key:hover { background: #eef3f8; border-color: #6f98bd; }",
    ".eqed-style-key--italic { font-style: italic; }",
    ".eqed-style-key--vector { font-weight: 700; }",
    ".eqed-tab-panel[hidden] { display: none; }",
    ".eqed-row { display: flex; flex-wrap: wrap; gap: 4px; align-items: center; margin-bottom: 6px; }",
    ".eqed-row-label { font-size: 0.78rem; color: var(--muted, #717170); min-width: 3.6em; }",
    ".eqed-key { font: inherit; font-size: 0.95rem; line-height: 1; padding: 5px 8px; border: 1px solid #cfcfcb; border-radius: 6px; background: #fff; cursor: pointer; min-width: 2.1em; }",
    ".eqed-key:hover { background: #eef3f8; border-color: #9db8d2; }",
    ".eqed-recent { display: block; width: fit-content; max-width: 100%; box-sizing: border-box; margin: 4px 0 2px; padding: 3px 8px; border: 1px solid #8fb7dc; border-radius: 5px; background: #dceeff; color: #163f63; font-family: var(--eqed-source-font); font-size: 0.88rem; white-space: pre-wrap; overflow-wrap: anywhere; }",
    ".eqed-recent[hidden] { display: none; }",
    ".eqed-recent--flash { animation: eqed-recent-flash 700ms ease-out; }",
    "@keyframes eqed-recent-flash { from { background: #8fc9ff; } to { background: #dceeff; } }",
    "@media (prefers-reduced-motion: reduce) { .eqed-recent--flash { animation: none; } }",
    ".eqed-source { width: 100%; box-sizing: border-box; font-family: var(--eqed-source-font); font-size: 0.95rem; padding: 8px 10px; border: 1px solid #cfcfcb; border-radius: 6px; margin-top: 4px; }",
    /* The native status bar's fifth part: the Tab / Enter / & hint. */
    ".eqed-hint { font-size: 0.76rem; color: var(--muted, #717170); margin: 2px 0 4px; overflow-wrap: anywhere; }",
    /* Why the equation cannot be shown, in place of the preview. */
    ".eqed-problem { margin: 0; padding: 8px 10px; border-left: 3px solid #c2793a; background: #fdf3e7; color: #6a3c11; font-size: 0.9rem; overflow-wrap: anywhere; }",
    /* Match the native canvas: white ground, sunken frame, top-left anchor. */
    ".eqed-preview { min-height: 5.5em; padding: 12px 14px; overflow-x: auto; background: #fff; border: 1px solid #c9c9c5; border-radius: 4px; box-shadow: inset 1px 1px 3px rgba(0,0,0,0.07); margin: 8px 0; }",
    ".eqed-preview mjx-container[display='true'] { text-align: left !important; margin: 0 !important; }",
    ".eqed-empty { color: var(--muted, #717170); font-size: 0.9rem; margin: 0.4em 0; }",
    ".eqed-actions { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; }",
    ".eqed-actions button { font: inherit; font-size: 0.88rem; padding: 6px 12px; border: 1px solid #cfcfcb; border-radius: 6px; background: #fff; cursor: pointer; }",
    ".eqed-actions button:hover { background: #eef3f8; border-color: #9db8d2; }",
    ".eqed-status { font-size: 0.85rem; color: var(--muted, #717170); }",
    ".eqed-build { margin-left: auto; align-self: center; font-size: 0.72rem; color: var(--muted, #a5a5a2); }"
  ];

  function installStyles() {
    if (document.getElementById("eqed-style")) return;
    var style = document.createElement("style");
    style.id = "eqed-style";
    style.textContent = CSS.join("");
    document.head.appendChild(style);
  }

  /* MathJax exposes \boldsymbol but not the bm package's shorter \bm name.
   * Keep \bm in editable/saved TeX and translate only at the renderer
   * boundary, so preview, SVG/PNG and Office MathML agree. */
  function leftAlignUnanchoredAligned(tex) {
    var begin = "\\begin{aligned}";
    var end = "\\end{aligned}";
    var output = "";
    var cursor = 0;

    function alignBody(body) {
      var rowStarts = [0];
      var braceDepth = 0;
      var environmentDepth = 0;
      var hasTopLevelTab = false;
      for (var i = 0; i < body.length; i += 1) {
        if (body.slice(i, i + 7) === "\\begin{") {
          environmentDepth += 1;
          i = body.indexOf("}", i + 7);
          if (i < 0) break;
          continue;
        }
        if (body.slice(i, i + 5) === "\\end{") {
          environmentDepth = Math.max(0, environmentDepth - 1);
          i = body.indexOf("}", i + 5);
          if (i < 0) break;
          continue;
        }
        if (body.charAt(i) === "{") { braceDepth += 1; continue; }
        if (body.charAt(i) === "}") {
          braceDepth = Math.max(0, braceDepth - 1);
          continue;
        }
        if (body.charAt(i) === "\\" && body.charAt(i + 1) === "&") {
          i += 1;
          continue;
        }
        if (braceDepth === 0 && environmentDepth === 0 &&
            body.charAt(i) === "&") {
          hasTopLevelTab = true;
        }
        if (braceDepth === 0 && environmentDepth === 0 &&
            body.charAt(i) === "\\" && body.charAt(i + 1) === "\\") {
          rowStarts.push(i + 2);
          i += 1;
        }
      }
      if (hasTopLevelTab) return body;
      for (var r = rowStarts.length - 1; r >= 0; r -= 1) {
        var at = rowStarts[r];
        while (at < body.length && /\s/.test(body.charAt(at))) at += 1;
        body = body.slice(0, at) + "& " + body.slice(at);
      }
      return body;
    }

    while (cursor < tex.length) {
      var start = tex.indexOf(begin, cursor);
      if (start < 0) { output += tex.slice(cursor); break; }
      output += tex.slice(cursor, start) + begin;
      var bodyStart = start + begin.length;
      var scan = bodyStart;
      var depth = 1;
      while (depth > 0) {
        var nextBegin = tex.indexOf(begin, scan);
        var nextEnd = tex.indexOf(end, scan);
        if (nextEnd < 0) {
          output += tex.slice(bodyStart);
          return output;
        }
        if (nextBegin >= 0 && nextBegin < nextEnd) {
          depth += 1;
          scan = nextBegin + begin.length;
        } else {
          depth -= 1;
          if (depth === 0) {
            var body = tex.slice(bodyStart, nextEnd);
            output += alignBody(leftAlignUnanchoredAligned(body)) + end;
            cursor = nextEnd + end.length;
          } else {
            scan = nextEnd + end.length;
          }
        }
      }
    }
    return output;
  }

  function mathJaxTex(tex) {
    return leftAlignUnanchoredAligned(
      tex.replace(/\\bm(?=\s*\{)/g, "\\boldsymbol"));
  }

  function splitUnanchoredOfficeRows(math) {
    function elements(node) {
      return Array.prototype.filter.call(
        node.children, function () { return true; });
    }

    function containsOnly(node, child) {
      return Array.prototype.every.call(node.childNodes, function (candidate) {
        return candidate === child ||
          (candidate.nodeType === 3 && candidate.nodeValue.trim() === "");
      });
    }

    function syntheticCells(table) {
      var alignment = (table.getAttribute("columnalign") || "")
        .trim().split(/\s+/);
      if (alignment.indexOf("right") < 0 ||
          alignment.indexOf("left") < 0) return null;
      var rows = elements(table).filter(function (row) {
        return row.localName === "mtr";
      });
      var cellsByRow = rows.map(function (row) {
        return elements(row).filter(function (cell) {
          return cell.localName === "mtd";
        });
      });
      if (!rows.length || !cellsByRow.every(function (cells) {
        return cells.length === 2 && cells[0].textContent.trim() === "";
      })) return null;
      return cellsByRow;
    }

    function onlyNestedSyntheticTable(cell) {
      var node = cell;
      while (true) {
        var children = elements(node);
        if (children.length !== 1 || !containsOnly(node, children[0])) {
          return null;
        }
        var child = children[0];
        if (child.localName === "mtable") {
          return syntheticCells(child) ? child : null;
        }
        if (child.localName !== "mrow" && child.localName !== "mstyle") {
          return null;
        }
        node = child;
      }
    }

    function collectLeafCells(table, output) {
      var cellsByRow = syntheticCells(table);
      if (!cellsByRow) return false;
      return cellsByRow.every(function (cells) {
        var content = cells[1];
        var nested = onlyNestedSyntheticTable(content);
        if (nested) return collectLeafCells(nested, output);
        if (content.querySelector("mtable")) return false;
        output.push(content);
        return true;
      });
    }

    var tables = Array.prototype.filter.call(
      math.querySelectorAll("mtable"), function (table) {
        for (var parent = table.parentNode; parent && parent !== math;
             parent = parent.parentNode) {
          if (parent.localName === "mtable") return false;
        }
        return true;
      });
    if (tables.length !== 1) return null;
    var table = tables[0];
    var only = table;
    while (only.parentNode !== math) {
      var parent = only.parentNode;
      if (!parent || (parent.localName !== "mrow" &&
          parent.localName !== "mstyle") ||
          elements(parent).length !== 1 ||
          elements(parent)[0] !== only || !containsOnly(parent, only)) {
        return null;
      }
      only = parent;
    }

    var leafCells = [];
    if (!collectLeafCells(table, leafCells) || leafCells.length < 2) return null;

    var serializer = new window.XMLSerializer();
    return leafCells.map(function (cell) {
      var rowMath = math.cloneNode(false);
      Array.prototype.forEach.call(cell.childNodes, function (child) {
        rowMath.appendChild(child.cloneNode(true));
      });
      return serializer.serializeToString(rowMath);
    });
  }

  /* Officeへ渡すMathMLの製品境界。MathJax固有属性や一時mstyleを落とし、
   * native版と同じ inline / 18 pt / large-operator 属性へ正規化する。
   * TeXの構造決定はMathJaxに任せるため、積分はmsubsup、総和は
   * munderoverとなり、native emitterも同じ規則を試験する。 */
  /* MathJax loads \boldsymbol and \cancel on demand, and until the package
   * arrives the synchronous tex2mml the Office copy needs throws
   * "MathJax retry".  Pull both in once at startup so the first copy of a
   * `\bm` or `\cancel` equation cannot fail on a package that is still in
   * flight.  Fire and forget: the copy path stays synchronous, which is what
   * keeps the user gesture alive for the clipboard write. */
  function warmAutoloadedMacros() {
    if (!window.MathJax || typeof window.MathJax.tex2mmlPromise !== "function") return;
    window.MathJax.tex2mmlPromise("\\boldsymbol{x}+\\cancel{x}", { display: false })
      .catch(function () { /* the sync path reports its own failure */ });
  }

  function officeMathMl(tex) {
    var raw = window.MathJax.tex2mml(
      "\\displaystyle " + mathJaxTex(tex), { display: false });
    var doc = new window.DOMParser().parseFromString(raw, "application/xml");
    if (doc.querySelector("parsererror")) throw new Error("invalid MathML");
    var math = doc.documentElement;
    math.setAttribute("display", "inline");
    math.setAttribute("mathsize", "18pt");

    Array.prototype.forEach.call(math.querySelectorAll("*"), function (node) {
      Array.prototype.slice.call(node.attributes).forEach(function (attribute) {
        if (attribute.name.indexOf("data-") === 0) {
          node.removeAttribute(attribute.name);
        }
      });
      if (node.getAttribute("stretchy") === "false") {
        node.removeAttribute("stretchy");
      }
    });
    Array.prototype.forEach.call(math.querySelectorAll("mstyle"), function (node) {
      while (node.firstChild) node.parentNode.insertBefore(node.firstChild, node);
      node.parentNode.removeChild(node);
    });
    Array.prototype.forEach.call(math.querySelectorAll("mo"), function (node) {
      if ("∑∏∐⋃⋂∫∬∭∮∯∰".indexOf(node.textContent) !== -1) {
        node.setAttribute("largeop", "true");
        node.setAttribute("movablelimits", "true");
      }
      /* MathJax emits U+2015 for both \overline and \underline.  Office
       * imports that as an ordinary glyph, while the native emitter's
       * macron/underscore forms become editable m:acc and m:bar objects. */
      var parent = node.parentNode && node.parentNode.localName;
      if ((parent === "mover" || parent === "munder") &&
          node.textContent === "―") {
        node.textContent = parent === "mover" ? "¯" : "_";
        node.removeAttribute("accent");
        node.setAttribute("stretchy", "true");
      }
    });
    var officeRows = splitUnanchoredOfficeRows(math);
    if (officeRows) {
      return officeRows.join(
        '<span style="font-size:18pt">&#160;</span><br>');
    }
    return new window.XMLSerializer().serializeToString(math);
  }

  /* ClipboardItemのtext/htmlはChromiumによってHTML文書として包み直される。
   * 同期copyイベントなら選択範囲と同じCF_HTML断片として渡せるため、
   * PowerPointの通常Ctrl+Vで18pt・左揃えを保持できる。 */
  function writeOfficeClipboard(html, tex) {
    var handled = false;
    function onCopy(event) {
      if (!event.clipboardData) return;
      event.preventDefault();
      event.clipboardData.setData("text/html", html);
      event.clipboardData.setData("text/plain", tex);
      handled = true;
    }
    document.addEventListener("copy", onCopy, true);
    var executed = false;
    try {
      executed = document.execCommand("copy");
    } catch (error) {
      executed = false;
    }
    document.removeEventListener("copy", onCopy, true);
    if (executed && handled) return Promise.resolve("copy-event");
    if (window.ClipboardItem && navigator.clipboard && navigator.clipboard.write) {
      return navigator.clipboard.write([new window.ClipboardItem({
        "text/html": new Blob([html], { type: "text/html" }),
        "text/plain": new Blob([tex], { type: "text/plain" })
      })]).then(function () { return "clipboard-item"; });
    }
    return Promise.reject(new Error("clipboard API missing"));
  }

  /* TeX → SVG。サイト共通の MathJax は CHTML 出力なので、初回押下時に
   * SVG 出力部品を動的に読み込み、グリフをパスとして埋め込む変換器を
   * 組み立てる（fontCache: "none"）。外部参照の無い SVG なので、後段の
   * canvas は汚染されずクリップボードへ渡せる。 */
  var svgConvert = null;
  function getSvgConverter() {
    if (svgConvert) return Promise.resolve(svgConvert);
    if (!window.MathJax || !window.MathJax.loader ||
        typeof window.MathJax.loader.load !== "function") {
      return Promise.reject(new Error("MathJax loader missing"));
    }
    return window.MathJax.loader.load("output/svg").then(function () {
      var SVG = window.MathJax._.output.svg_ts.SVG;
      var TeXFont = window.MathJax._.output.svg.fonts.tex_ts.TeXFont;
      var mathjax = window.MathJax._.mathjax.mathjax;
      var jax = new SVG({ fontCache: "none", font: new TeXFont({}) });
      var doc = mathjax.document("", {
        InputJax: window.MathJax.startup.input[0],
        OutputJax: jax
      });
      svgConvert = function (tex) {
        var node = doc.convert(mathJaxTex(tex), { display: true });
        return node.querySelector("svg");
      };
      return svgConvert;
    });
  }

  /* TeX → 単体で開ける SVG テキスト。寸法は ex → px に直して書く
   * （PowerPoint などの取り込み側は ex 単位を解釈しない）。 */
  function texToSvgText(tex) {
    return getSvgConverter().then(function (convert) {
      var svg = convert(tex);
      if (!svg) throw new Error("no svg");
      svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      svg.setAttribute("width", Math.max(1, Math.ceil(parseFloat(svg.getAttribute("width")) * 8)) + "px");
      svg.setAttribute("height", Math.max(1, Math.ceil(parseFloat(svg.getAttribute("height")) * 8)) + "px");
      return '<?xml version="1.0" encoding="UTF-8"?>\n' + svg.outerHTML;
    });
  }

  /* TeX → PNG Blob。1ex ≈ 8px（本文16px基準）に SCALE を掛けた解像度で
   * 描く。背景は透明 -- TeXclip と同じく、貼り付け先の地色に馴染む。 */
  var PNG_SCALE = 4;
  function texToPngBlob(tex) {
    return getSvgConverter().then(function (convert) {
      var svg = convert(tex);
      if (!svg) throw new Error("no svg");
      svg.setAttribute("xmlns", "http://www.w3.org/2000/svg");
      var w = Math.max(1, Math.ceil(parseFloat(svg.getAttribute("width")) * 8 * PNG_SCALE));
      var h = Math.max(1, Math.ceil(parseFloat(svg.getAttribute("height")) * 8 * PNG_SCALE));
      var url = URL.createObjectURL(new Blob([svg.outerHTML], { type: "image/svg+xml" }));
      return new Promise(function (resolve, reject) {
        var img = new Image();
        img.onload = function () {
          URL.revokeObjectURL(url);
          var canvas = document.createElement("canvas");
          canvas.width = w;
          canvas.height = h;
          canvas.getContext("2d").drawImage(img, 0, 0, w, h);
          canvas.toBlob(function (blob) {
            if (blob) resolve(blob);
            else reject(new Error("toBlob failed"));
          }, "image/png");
        };
        img.onerror = function () {
          URL.revokeObjectURL(url);
          reject(new Error("svg image load failed"));
        };
        img.src = url;
      });
    });
  }

  function el(tag, className, text) {
    var node = document.createElement(tag);
    if (className) node.className = className;
    if (text) node.textContent = text;
    return node;
  }

  /* 空欄 {} を探して中へ。無ければ null。 */
  function nextHole(text, from, backwards) {
    if (backwards) {
      /* A caret inside a hole sits at from with "{" at from - 1, so the
       * previous hole can start no later than from - 2.  Starting from
       * text.length - 3 (the old bound) skipped a hole at the very end. */
      for (var i = from - 2; i >= 0; i--) {
        if (text.charAt(i) === "{" && text.charAt(i + 1) === "}") return i + 1;
      }
      return null;
    }
    for (var j = from; j + 1 < text.length; j++) {
      if (text.charAt(j) === "{" && text.charAt(j + 1) === "}") return j + 1;
    }
    return null;
  }

  /* Native GUI_SPEC §3.5 source layout, applied to inserted templates: a
   * line break after \begin{...}, before \end{...}, and after a row break
   * \\, with two spaces of indent per environment depth.  Only whitespace
   * that TeX ignores is touched, so the rendered equation is unchanged.
   * Text the user typed is never re-flowed under the caret. */
  function prettyTex(raw) {
    var out = "";
    var depth = 0;
    var lineStart = true;
    var i = 0;
    function newline() {
      out = out.replace(/[ \t]+$/, "") + "\n" + new Array(depth + 1).join("  ");
      lineStart = true;
    }
    while (i < raw.length) {
      var c = raw.charAt(i);
      if (c === "\n" || c === "\r" || (lineStart && (c === " " || c === "\t"))) {
        i += 1;
        continue;
      }
      var close = -1;
      if (raw.slice(i, i + 5) === "\\end{" && (close = raw.indexOf("}", i)) >= 0) {
        depth = Math.max(0, depth - 1);
        if (!lineStart) newline();
        out += raw.slice(i, close + 1);
        lineStart = false;
        i = close + 1;
        continue;
      }
      if (raw.slice(i, i + 7) === "\\begin{" && (close = raw.indexOf("}", i)) >= 0) {
        if (raw.charAt(close + 1) === "{") {
          /* \begin{array}{cc}: the column spec stays on the \begin line. */
          var spec = raw.indexOf("}", close + 1);
          if (spec >= 0) close = spec;
        }
        /* The rule breaks after \begin, not before it, so a nested
         * environment stays on the line that introduces it. */
        out += raw.slice(i, close + 1);
        i = close + 1;
        depth += 1;
        newline();
        continue;
      }
      if (c === "\\" && raw.charAt(i + 1) === "\\") {
        out += "\\\\";
        i += 2;
        newline();
        continue;
      }
      if (c === "\\" && i + 1 < raw.length) {
        /* Keep control symbols such as \{ and \} intact. */
        out += raw.slice(i, i + 2);
        lineStart = false;
        i += 2;
        continue;
      }
      out += c;
      lineStart = false;
      i += 1;
    }
    return out.replace(/[ \t]+$/, "");
  }

  /* Row environments: the ones whose \\ means "next row" and whose & means
   * "next column".  Enter adds rows to these. */
  var ROW_ENVIRONMENT = /^(aligned|alignedat|gathered|split|cases|dcases|rcases|array|(?:p|b|B|v|V|small)?matrix)$/;

  function environmentToken(text, from) {
    /* End of a \begin{name} token at from, including an array column spec. */
    var close = text.indexOf("}", from);
    if (close < 0) return from;
    if (text.charAt(close + 1) === "{" && /^\\begin\{array\}/.test(text.slice(from))) {
      var spec = text.indexOf("}", close + 1);
      if (spec >= 0) close = spec;
    }
    return close + 1;
  }

  /* Environments open at pos, outermost first.  A position inside a
   * \begin token counts as the start of that body; inside an \end token it
   * counts as the end of that body. */
  function openEnvironments(text, pos) {
    var token = /\\(begin|end)\{[^}]*\}/g;
    var stack = [];
    var match;
    while ((match = token.exec(text)) !== null) {
      var from = match.index;
      if (from >= pos) break;
      if (match[1] === "begin") {
        stack.push({
          name: match[0].slice(7, -1),
          begin: from,
          body: environmentToken(text, from)
        });
      } else if (pos < from + match[0].length) {
        break;
      } else if (stack.length) {
        stack.pop();
      }
    }
    return stack;
  }

  /* Index of the \end token that closes env (text.length when unclosed). */
  function environmentEnd(text, env) {
    var token = /\\(begin|end)\{[^}]*\}/g;
    token.lastIndex = env.body;
    var level = 0;
    var match;
    while ((match = token.exec(text)) !== null) {
      if (match[1] === "begin") level += 1;
      else if (level === 0) return match.index;
      else level -= 1;
    }
    return text.length;
  }

  /* The single row environment that makes up the whole source, if any:
   * Enter outside it (before \begin or after \end) still adds a row to it
   * instead of nesting a second aligned. */
  function soleRowEnvironment(text) {
    var begin = text.indexOf("\\begin{");
    if (begin < 0 || text.slice(0, begin).trim() !== "") return null;
    var close = text.indexOf("}", begin);
    if (close < 0) return null;
    var env = { name: text.slice(begin + 7, close), begin: begin, body: environmentToken(text, begin) };
    if (!ROW_ENVIRONMENT.test(env.name)) return null;
    var end = environmentEnd(text, env);
    var after = text.indexOf("}", end);
    if (end >= text.length || after < 0 || text.slice(after + 1).trim() !== "") return null;
    return env;
  }

  /* The row of body [bodyStart, bodyEnd) that contains pos: its bounds, the
   * column (top-level & count before pos) and the brace/environment depth
   * at pos.  A pos inside a \\ token is moved in front of it. */
  function rowContext(text, bodyStart, bodyEnd, pos) {
    var depth = 0;
    var depthAtPos = null;
    var rowStart = bodyStart;
    var rowEnd = bodyEnd;
    var column = 0;
    var i = bodyStart;
    while (i < bodyEnd) {
      if (depthAtPos === null && i >= pos) depthAtPos = depth;
      var c = text.charAt(i);
      if (c === "\\") {
        if (text.charAt(i + 1) === "\\") {
          if (i < pos && pos < i + 2) pos = i;
          if (depth === 0) {
            if (i + 2 <= pos) {
              rowStart = i + 2;
              column = 0;
            } else {
              rowEnd = i;
              break;
            }
          }
          i += 2;
          continue;
        }
        if (text.slice(i, i + 7) === "\\begin{") {
          depth += 1;
          i = environmentToken(text, i);
          continue;
        }
        if (text.slice(i, i + 5) === "\\end{") {
          depth = Math.max(0, depth - 1);
          var close = text.indexOf("}", i);
          i = close >= 0 ? close + 1 : i + 5;
          continue;
        }
        i += 2;                        /* \{ \} \& and the first letter of a control word */
        continue;
      }
      if (c === "{") depth += 1;
      else if (c === "}") depth = Math.max(0, depth - 1);
      else if (c === "&" && depth === 0 && i < pos) column += 1;
      i += 1;
    }
    if (depthAtPos === null) depthAtPos = depth;
    var contentStart = rowStart;
    while (contentStart < rowEnd && /\s/.test(text.charAt(contentStart))) contentStart += 1;
    var contentEnd = rowEnd;
    while (contentEnd > contentStart && /\s/.test(text.charAt(contentEnd - 1))) contentEnd -= 1;
    return {
      pos: pos,
      rowStart: rowStart,
      contentStart: contentStart,
      contentEnd: contentEnd,
      column: column,
      depth: depthAtPos
    };
  }

  function wrapInAligned(text, pos) {
    var row = rowContext(text, 0, text.length, pos);
    var first;
    var second;
    if (row.depth > 0) {
      /* Enter inside \frac{}{} or ^{} at top level: keep the expression whole
       * as row one and open an empty row two. */
      first = text.trim();
      second = "";
    } else {
      first = text.slice(0, row.pos).trim();
      second = text.slice(row.pos).trim();
    }
    var open = "\\begin{aligned}\n  ";
    var between = " \\\\\n  ";
    var value = open + (first || "{}") + between + (second || "{}") + "\n\\end{aligned}";
    return {
      value: value,
      caret: open.length + (first || "{}").length + between.length + (second ? 0 : 1)
    };
  }

  /* Structural Enter (native GUI_SPEC §3.5), as a pure function like
   * composeInsertion: the selection is deleted, then
   *   - in a row environment, the row is split at the caret: the rest of the
   *     row moves to a new row that keeps the column (leading "& "), and an
   *     empty row gets a {} hole;
   *   - inside \frac{}{}, ^{}, \text{} or a nested environment the group is
   *     never split: an empty {} row is opened after the current row;
   *   - at top level the whole source is wrapped in aligned;
   *   - before \begin or after \end of a source that is one row environment,
   *     an empty first or last row is added instead of nesting a second one. */
  function composeRowBreak(value, selStart, selEnd) {
    selStart = Math.max(0, Math.min(selStart, value.length));
    selEnd = Math.max(selStart, Math.min(selEnd, value.length));
    var text = value.slice(0, selStart) + value.slice(selEnd);
    var pos = selStart;
    var stack = openEnvironments(text, pos);
    var env = null;
    var depth = 0;
    for (var s = stack.length - 1; s >= 0; s--) {
      if (ROW_ENVIRONMENT.test(stack[s].name)) {
        env = stack[s];
        depth = s + 1;
        break;
      }
    }
    var firstRow = false;
    if (!env) {
      env = soleRowEnvironment(text);
      if (env) {
        depth = 1;
        firstRow = pos <= env.begin;
        pos = firstRow ? env.body : environmentEnd(text, env);
      }
    }
    if (!env) return wrapInAligned(text, pos);
    var bodyEnd = environmentEnd(text, env);
    pos = Math.max(env.body, Math.min(pos, bodyEnd));
    var indent = new Array(depth + 1).join("  ");
    var row = rowContext(text, env.body, bodyEnd, pos);
    var head;
    var ins;
    if (firstRow) {
      head = text.slice(0, row.contentStart);
      ins = "{} \\\\\n" + indent;
      return { value: head + ins + text.slice(row.contentStart), caret: head.length + 1 };
    }
    if (row.depth > 0) {
      head = text.slice(0, row.contentEnd);
      ins = " \\\\\n" + indent + "{}";
      return { value: head + ins + text.slice(row.contentEnd), caret: head.length + ins.length - 1 };
    }
    /* Split inside the row's own content: the indentation in front of it and
     * the newline that keeps \end on its own line are layout, not material to
     * move into the new row. */
    var at = Math.max(row.contentStart, Math.min(row.pos, row.contentEnd));
    var left = at;
    while (left > row.contentStart && /\s/.test(text.charAt(left - 1))) left -= 1;
    var right = at;
    while (right < row.contentEnd && /\s/.test(text.charAt(right))) right += 1;
    var prefix = "";
    for (var k = 0; k < row.column; k++) prefix += "& ";
    var hole = right >= row.contentEnd ? "{}" : "";
    head = text.slice(0, left);
    ins = " \\\\\n" + indent + prefix + hole;
    return {
      value: head + ins + text.slice(right),
      caret: head.length + ins.length - (hole ? 1 : 0)
    };
  }

  /* Replace input[start, end) with text through the browser's editing
   * command so the change stays in the undo history (Ctrl+Z).  Assigning
   * input.value discards that history, which made every palette click
   * irreversible.  setRangeText is the standards path when execCommand is
   * unavailable; it does not register in the history on every engine. */
  function replaceRangeUndoable(input, start, end, text) {
    var old = input.value;
    var expected = old.slice(0, start) + text + old.slice(end);
    input.focus();
    input.setSelectionRange(start, end);
    var done = false;
    try {
      done = document.execCommand(text ? "insertText" : "delete", false, text);
    } catch (error) {
      done = false;
    }
    if (done && input.value === expected) return;
    if (input.value === old) {
      input.setRangeText(text, start, end, "end");
    } else {
      input.value = expected;
    }
    input.dispatchEvent(new Event("input", { bubbles: true }));
  }

  /* Apply a {value, caret} edit as the smallest replacement that turns the
   * current value into the new one, then place the caret. */
  function applyEdit(input, edit) {
    var old = input.value;
    var value = edit.value;
    var limit = Math.min(old.length, value.length);
    var head = 0;
    while (head < limit && old.charAt(head) === value.charAt(head)) head += 1;
    var tail = 0;
    while (tail < limit - head &&
           old.charAt(old.length - 1 - tail) === value.charAt(value.length - 1 - tail)) {
      tail += 1;
    }
    if (head < old.length - tail || head < value.length - tail) {
      replaceRangeUndoable(input, head, old.length - tail, value.slice(head, value.length - tail));
    } else {
      input.focus();
    }
    input.setSelectionRange(edit.caret, edit.caret);
  }

  /* An unclosed group is the beginner's most frequent mistake, and it is the
   * one MathJax handles worst for a learner: it refuses the whole expression,
   * so the preview is left holding raw TeX with no hint of what went wrong.
   * Name the problem, and say where it is, before handing the source over.
   * The native editor cannot reach this state - its parser keeps brace
   * nesting balanced by construction - so this is a Web-only obligation. */
  function braceProblem(tex) {
    var depth = 0;
    for (var i = 0; i < tex.length; i++) {
      var c = tex.charAt(i);
      if (c === "\\") { i += 1; continue; }   /* \{ and \} are symbols */
      if (c === "{") depth += 1;
      else if (c === "}") {
        depth -= 1;
        if (depth < 0) {
          return (i + 1) + " 文字目の } に対応する { がありません。";
        }
      }
    }
    if (depth > 0) {
      return "閉じていない { が " + depth + " 個あります。";
    }
    return null;
  }

  /* TeX reads ' as a superscript, so a prime placed straight after another
   * superscript is a double exponent: `a^{2}'` fails to convert with
   * "Prime causes double exponent: use braces to clarify".  The structural
   * editor cannot produce that because it attaches the prime to the base; a
   * source pane can, because the palette inserts at the caret.  Detect the
   * case so the insertion can carry MathJax's own remedy, an empty group. */
  function endsWithSuperscript(text) {
    var i = text.length - 1;
    while (i >= 0 && /\s/.test(text.charAt(i))) i -= 1;
    if (i < 0) return false;
    if (text.charAt(i) === "'") return false;      /* a'' is legal */
    var before = i;
    if (text.charAt(i) === "}" && text.charAt(i - 1) !== "\\") {
      var depth = 0;
      for (; i >= 0; i--) {
        var c = text.charAt(i);
        if (text.charAt(i - 1) === "\\") continue;  /* \{ and \} are symbols */
        if (c === "}") depth += 1;
        else if (c === "{") {
          depth -= 1;
          if (depth === 0) break;
        }
      }
      if (i <= 0) return false;
      before = i;
    }
    var j = before - 1;
    while (j >= 0 && /\s/.test(text.charAt(j))) j -= 1;
    return text.charAt(j) === "^";
  }

  /* Pure insertion contract, shared by the live textarea and the Node CI
   * test.  Keeping caret arithmetic out of DOM event code makes selection
   * wrapping a tested result instead of a browser-manual assumption. */
  function composeInsertion(value, start, end, snippet) {
    start = Math.max(0, Math.min(start, value.length));
    end = Math.max(start, Math.min(end, value.length));
    var before = value.slice(0, start);
    var selected = value.slice(start, end);
    var after = value.slice(end);
    if (/^'+$/.test(snippet)) {
      /* A prime is its own superscript, so it never takes a hole and never
       * wraps a selection.  Clarify with an empty group when the caret sits
       * right after another superscript, which is exactly what MathJax asks
       * for; the rendered result is unchanged where no group is needed. */
      var prime = (endsWithSuperscript(before) ? "{}" : "") + snippet;
      return { value: before + prime + after, caret: before.length + prime.length };
    }
    var body = snippet;
    var caret;
    var firstHole = body.indexOf("{}");
    if (selected && firstHole >= 0) {
      body = body.slice(0, firstHole + 1) + selected + body.slice(firstHole + 1);
      var second = nextHole(body, firstHole + 1 + selected.length + 1, false);
      caret = second !== null ? second : firstHole + 1 + selected.length + 1;
    } else if (firstHole >= 0) {
      caret = firstHole + 1;
    } else {
      caret = body.length;
    }
    return {
      value: before + body + after,
      caret: before.length + caret
    };
  }

  function init(root) {
    installStyles();
    var input = root.querySelector(".eqed-source");
    var preview = root.querySelector(".eqed-preview");
    var status = root.querySelector(".eqed-status");
    /* It now carries the TeX of the cell under the cursor on purpose. */
    status.setAttribute("data-tex-literal-ok", "true");
    var paletteHost = root.querySelector(".eqed-palettes");
    var recent = el("output", "eqed-recent");
    recent.hidden = true;
    recent.setAttribute("aria-live", "polite");
    /* Shows TeX source on purpose (the learning surface); the homepage QA
     * must not read it as a MathJax leak. */
    recent.setAttribute("data-tex-literal-ok", "true");
    input.parentNode.insertBefore(recent, input);
    var hint = el("p", "eqed-hint",
      "Tab / Shift+Tab: 次・前の空欄へ　" +
      "Enter: 数式の改行 \\\\（aligned・行列・場合分けの行が増える）　" +
      "Shift+Enter: ソースだけ改行　Ctrl+Z: 取り消し");
    hint.setAttribute("data-tex-literal-ok", "true");
    input.parentNode.insertBefore(hint, input.nextSibling);
    var timer = null;

    /* Replace the preview with the reason the equation cannot be shown.  The
     * source stays untouched in the pane above, so the reader can compare the
     * message against what they typed. */
    function showRenderProblem(reason) {
      preview.innerHTML = "";
      var note = el("p", "eqed-problem");
      note.appendChild(el("strong", "", "数式を表示できません"));
      note.appendChild(document.createTextNode(" — " + reason));
      preview.appendChild(note);
    }

    /* MathJax rejects more than unbalanced braces, and when it does the page
     * is left showing raw TeX.  Ask it for the reason and put that on screen
     * instead, unless the source has moved on in the meantime. */
    function reportUnrendered(source) {
      if (input.value.trim() !== source) return;
      if (preview.querySelector("mjx-container")) return;
      if (!window.MathJax || typeof window.MathJax.tex2mmlPromise !== "function") return;
      window.MathJax.tex2mmlPromise(mathJaxTex(source), { display: true }).then(
        function () { /* converted after all: leave the preview alone */ },
        function (error) {
          if (input.value.trim() !== source) return;
          showRenderProblem(String((error && error.message) || error));
        }
      );
    }

    function render() {
      var tex = input.value.trim();
      preview.innerHTML = "";
      if (!tex) {
        preview.appendChild(el("p", "eqed-empty", "ここに数式が表示されます"));
        return;
      }
      var problem = braceProblem(tex);
      if (problem) {
        showRenderProblem(problem);
        return;
      }
      var box = el("div", "");
      box.textContent = "\\[" + mathJaxTex(tex) + "\\]";
      preview.appendChild(box);
      if (window.MathJax && typeof window.MathJax.typesetClear === "function") {
        window.MathJax.typesetClear([preview]);
      }
      if (window.SugaharaMath && typeof window.SugaharaMath.typeset === "function") {
        window.SugaharaMath.typeset(preview);
      }
      /* The shared typeset chain swallows its own failures, so the only
       * reliable signal is the absence of a rendered container afterwards. */
      window.setTimeout(function () { reportUnrendered(tex); }, 600);
    }

    function scheduleRender() {
      if (timer) window.clearTimeout(timer);
      timer = window.setTimeout(render, 250);
    }

    function say(message) {
      status.textContent = message;
      if (message) {
        window.setTimeout(function () {
          if (status.textContent === message) status.textContent = "";
        }, 4000);
      }
    }

    /* The native status bar explains the cell under the cursor before it is
     * inserted.  A title tooltip cannot do that on a touch screen or for a
     * keyboard user, and a two-letter face such as "sf" or "tt" says nothing
     * on its own, so mirror the native behaviour in the status line. */
    function showKeyHelp(text) {
      status.textContent = text;
    }

    function clearKeyHelp(text) {
      if (status.textContent === text) status.textContent = "";
    }

    function showRecentInsertion(snippet) {
      recent.textContent = snippet.replace(/\s+$/, "");
      recent.hidden = !recent.textContent;
      recent.classList.remove("eqed-recent--flash");
      /* Restart the highlight when the same palette item is pressed twice. */
      void recent.offsetWidth;
      recent.classList.add("eqed-recent--flash");
    }

    function clearRecentInsertion() {
      recent.hidden = true;
      recent.classList.remove("eqed-recent--flash");
    }

    /* 挿入。選択があり、断片が空欄を持つなら選択を最初の空欄へ包む。
     * 環境テンプレートはネイティブ版と同じ改行・字下げで入れる。 */
    function insert(snippet) {
      if (snippet.indexOf("\\begin{") >= 0) snippet = prettyTex(snippet);
      var edit = composeInsertion(
        input.value, input.selectionStart, input.selectionEnd, snippet);
      applyEdit(input, edit);
      showRecentInsertion(snippet);
      scheduleRender();
    }

    var idPrefix = "eqed-palette";
    var paletteHead = el("div", "eqed-palette-head");
    var tabList = el("div", "eqed-tabs");
    var styleBar = el("div", "eqed-stylebar");
    var panelHost = el("div", "eqed-tab-panels");
    var tabButtons = [];
    var tabPanels = [];
    tabList.setAttribute("role", "tablist");
    tabList.setAttribute("aria-label", "数式記号の分類");
    styleBar.setAttribute("role", "toolbar");
    styleBar.setAttribute("aria-label", "数式の文字スタイル");
    paletteHost.setAttribute("aria-label", "数式記号パレット");
    paletteHead.appendChild(tabList);
    paletteHead.appendChild(styleBar);
    paletteHost.appendChild(paletteHead);
    paletteHost.appendChild(panelHost);

    MATH_ALPHABETS.forEach(function (alphabet) {
      var button = el("button",
        "eqed-style-key eqed-style-key--" + alphabet.css, alphabet.label);
      button.type = "button";
      button.setAttribute("data-tex-literal-ok", "true");
      button.setAttribute("aria-label", alphabet.name + " " + alphabet.snippet);
      button.title = alphabet.snippet;
      button.addEventListener("click", function () { insert(alphabet.snippet); });
      styleBar.appendChild(button);
    });

    function activatePaletteTab(index, moveFocus) {
      tabButtons.forEach(function (button, buttonIndex) {
        var selected = buttonIndex === index;
        button.setAttribute("aria-selected", selected ? "true" : "false");
        button.tabIndex = selected ? 0 : -1;
        tabPanels[buttonIndex].hidden = !selected;
      });
      if (moveFocus) tabButtons[index].focus();
    }

    PALETTE_TABS.forEach(function (group, groupIndex) {
      var tabId = idPrefix + "-tab-" + group.id;
      var panelId = idPrefix + "-panel-" + group.id;
      var tab = el("button", "eqed-tab", group.label);
      tab.type = "button";
      tab.id = tabId;
      tab.setAttribute("role", "tab");
      tab.setAttribute("aria-controls", panelId);
      tab.addEventListener("click", function () {
        activatePaletteTab(groupIndex, false);
      });
      tab.addEventListener("keydown", function (event) {
        var target = null;
        if (event.key === "ArrowRight") {
          target = (groupIndex + 1) % PALETTE_TABS.length;
        } else if (event.key === "ArrowLeft") {
          target = (groupIndex + PALETTE_TABS.length - 1) %
            PALETTE_TABS.length;
        } else if (event.key === "Home") {
          target = 0;
        } else if (event.key === "End") {
          target = PALETTE_TABS.length - 1;
        }
        if (target === null) return;
        event.preventDefault();
        activatePaletteTab(target, true);
      });
      tabList.appendChild(tab);
      tabButtons.push(tab);

      var panel = el("div", "eqed-tab-panel");
      panel.id = panelId;
      panel.setAttribute("role", "tabpanel");
      panel.setAttribute("aria-labelledby", tabId);
      PALETTES.forEach(function (palette) {
        if (group.palettes.indexOf(palette.label) < 0) return;
        var row = el("div", "eqed-row");
        row.setAttribute("role", "toolbar");
        row.setAttribute("aria-label", palette.label);
        row.appendChild(el("span", "eqed-row-label", palette.label));
        palette.items.forEach(function (item) {
          var button = el("button", "eqed-key", item[0]);
          var description = item[2]
            ? item[2] + "  " + item[1].trim()
            : item[1].trim();
          button.type = "button";
          button.setAttribute("data-tex-literal-ok", "true");
          /* title only: it becomes the accessible DESCRIPTION while the
           * visible face stays the accessible NAME.  An aria-label here would
           * replace that name, so "arg min" would no longer address the key
           * for a screen reader, for voice control, or for the homepage QA. */
          button.title = description;
          button.addEventListener("mouseenter", function () { showKeyHelp(description); });
          button.addEventListener("focus", function () { showKeyHelp(description); });
          button.addEventListener("mouseleave", function () { clearKeyHelp(description); });
          button.addEventListener("blur", function () { clearKeyHelp(description); });
          button.addEventListener("click", function () { insert(item[1]); });
          row.appendChild(button);
        });
        panel.appendChild(row);
      });
      panelHost.appendChild(panel);
      tabPanels.push(panel);
    });
    activatePaletteTab(0, false);

    input.addEventListener("input", function () {
      clearRecentInsertion();
      scheduleRender();
    });
    input.addEventListener("keydown", function (event) {
      if (event.key !== "Tab") return;
      var pos = event.shiftKey
        ? nextHole(input.value, input.selectionStart, true)
        : nextHole(input.value, input.selectionEnd, false);
      if (pos === null) return;
      event.preventDefault();
      input.setSelectionRange(pos, pos);
    });
    /* Enter は数式の行区切り（ネイティブ GUI_SPEC §3.5 の構造的 Enter）。
     * Shift+Enter はブラウザ標準のソース改行のまま。日本語入力の確定
     * Enter（isComposing / keyCode 229）は IME に渡す。 */
    input.addEventListener("keydown", function (event) {
      if (event.key !== "Enter" || event.shiftKey || event.ctrlKey ||
          event.altKey || event.metaKey) return;
      if (event.isComposing || event.keyCode === 229) return;
      event.preventDefault();
      applyEdit(input, composeRowBreak(
        input.value, input.selectionStart, input.selectionEnd));
      scheduleRender();
    });

    /* Office へのコピー。text/html に MathML、text/plain に TeX の2形式
     * だけを積む（動作実績のある構成）。Word / PowerPoint は HTML の
     * MathML をネイティブ数式として取り込む。
     *
     * 教訓（2026-08-24）: ここに image/png を同居させてはいけない。
     * PowerPoint が HTML より画像を優先して掴むことがあり、透過背景の
     * PNG は暗いスライドで見えない。「貼り付けたのに何も見えない」に
     * なる。画像が欲しいときは専用の「PNGでコピー」を使う。 */
    root.querySelector(".eqed-copy-office").addEventListener("click", function () {
      var tex = input.value.trim();
      if (!tex) { say("数式が空です"); return; }
      if (!window.MathJax || typeof window.MathJax.tex2mml !== "function") {
        say("MathMLの生成にはページの再読み込みが必要です");
        return;
      }
      var mml;
      try {
        /* inline の MathML として渡す。display(block) だと Office は
         * 「数式段落」を作り、段落の左寄せを無視して中央に置く。
         * 数式の後ろの NBSP が段落を「文中数式入りテキスト」に格下げし、
         * 左寄せが効く（2026-08-24 に実機で実測: block はボックス中央
         * x=220、この形は左端 x=8）。前置きのゼロ幅スペースは Word で
         * 豆腐（□）に見えたため、後置きの NBSP（どのフォントにもある
         * 不可視の空白）に置き換えた。\displaystyle は Office では
         * 落ちるが害はなく、正しく解釈する他アプリでは表示を保つ。 */
        mml = officeMathMl(tex);
      } catch (error) {
        say("この数式はMathMLに変換できませんでした");
        return;
      }
      /* native版と同じinline 18 pt MathMLを同期CF_HTMLで渡す。
       * PowerPointにOMMLを直接渡すと、MathMLからのOffice変換と
       * 総和記号の大きさや上下限の位置が異なるため使わない。 */
      /* The trailing inline sentinel keeps the equation out of a centred
       * display paragraph. Give that otherwise invisible run the same 18 pt
       * size as the equation so the caret and next insertion agree. */
      var html = mml + '<span style="font-size:18pt">&#160;</span>';
      writeOfficeClipboard(html, tex).then(
        function () { say("18 pt・左揃えのPowerPoint数式をコピーしました"); },
        function () { say("コピーできませんでした（ブラウザの権限を確認してください）"); }
      );
    });
    root.querySelector(".eqed-copy-display").addEventListener("click", function () {
      navigator.clipboard.writeText("$$\n" + input.value.trim() + "\n$$").then(
        function () { say("$$付きでコピーしました"); },
        function () { say("コピーできませんでした（手動で選択してください）"); }
      );
    });
    root.querySelector(".eqed-copy-equation").addEventListener("click", function () {
      var tex = input.value.trim();
      if (!tex) { say("数式が空です"); return; }
      navigator.clipboard.writeText(
        "\\begin{equation}\n" + tex + "\n\\end{equation}"
      ).then(
        function () { say("equation環境付きでコピーしました"); },
        function () { say("コピーできませんでした（手動で選択してください）"); }
      );
    });
    /* 画像（PNG）でコピー。MathML を受け付けない Google スライド等へは
     * 画像で貼る。ClipboardItem へ Promise を渡し、クリック（ユーザー
     * 操作）と同期でクリップボード書き込みを開始する（Safari の制約）。 */
    root.querySelector(".eqed-copy-png").addEventListener("click", function () {
      var tex = input.value.trim();
      if (!tex) { say("数式が空です"); return; }
      if (!window.ClipboardItem || !navigator.clipboard || !navigator.clipboard.write) {
        say("このブラウザは画像コピーに対応していません");
        return;
      }
      say("PNGを作成中…");
      var item;
      try {
        item = new window.ClipboardItem({ "image/png": texToPngBlob(tex) });
      } catch (error) {
        /* Promise を受け取れない古い ClipboardItem 実装への退避。 */
        texToPngBlob(tex).then(function (blob) {
          return navigator.clipboard.write([
            new window.ClipboardItem({ "image/png": blob })
          ]);
        }).then(
          function () { say("PNGをコピーしました"); },
          function () { say("PNGコピーできませんでした"); }
        );
        return;
      }
      navigator.clipboard.write([item]).then(
        function () { say("PNGをコピーしました。スライドにそのまま貼れます"); },
        function () { say("PNGコピーできませんでした（ブラウザの権限を確認してください）"); }
      );
    });
    /* SVG で保存。ベクトルのままなので PowerPoint（挿入→画像）や
     * Inkscape で劣化なく使える。SVG をクリップボード経由で受け取る
     * アプリは無いため、コピーではなくファイル保存にする。 */
    root.querySelector(".eqed-save-svg").addEventListener("click", function () {
      var tex = input.value.trim();
      if (!tex) { say("数式が空です"); return; }
      texToSvgText(tex).then(function (text) {
        var url = URL.createObjectURL(new Blob([text], { type: "image/svg+xml" }));
        var a = document.createElement("a");
        a.href = url;
        a.download = "equation.svg";
        document.body.appendChild(a);
        a.click();
        a.remove();
        window.setTimeout(function () { URL.revokeObjectURL(url); }, 1000);
        say("equation.svg を保存しました");
      }, function () {
        say("SVGを作成できませんでした");
      });
    });
    root.querySelector(".eqed-clear").addEventListener("click", function () {
      input.value = "";
      clearRecentInsertion();
      input.focus();
      render();
    });

    root.querySelector(".eqed-actions")
        .appendChild(el("span", "eqed-build", "build " + BUILD));

    render();
    if (window.MathJax && window.MathJax.startup && window.MathJax.startup.promise) {
      window.MathJax.startup.promise.then(warmAutoloadedMacros).catch(function () {});
    } else {
      warmAutoloadedMacros();
    }
  }

  function ready() {
    var root = document.querySelector("[data-equation-editor]");
    if (root) init(root);
  }
  if (typeof module !== "undefined" && module.exports) {
    module.exports = {
      composeInsertion: composeInsertion,
      composeRowBreak: composeRowBreak,
      braceProblem: braceProblem,
      prettyTex: prettyTex,
      nextHole: nextHole,
      mathJaxTex: mathJaxTex,
      mathAlphabets: MATH_ALPHABETS.map(function (item) {
        return { label: item.label, snippet: item.snippet };
      })
    };
  }
  if (typeof document === "undefined") return;
  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", ready);
  } else {
    ready();
  }
})();
