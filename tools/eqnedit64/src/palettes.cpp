/* The eighteen palettes.  Categories follow Eqnedt32's own toolbar, which
 * groups by what a reader is looking for rather than by internal node type.
 *
 * Rules this table has to keep (tests/test_palettes.py checks all of them):
 *   - every command in the symbol table appears exactly once;
 *   - every template kind appears exactly once;
 *   - no command appears on two palettes.
 */
#include "palettes.h"

namespace eqnedit {
namespace {

using I = PaletteItem;

I sym(const char* c, const char* f, const char* l) {
    return I{std::string("symbol.") + c, f, l};
}
I tpl(const char* c, const char* f, const char* l) {
    return I{std::string("template.") + c, f, l};
}
I raw(const char* c, const char* f, const char* l) {
    return I{std::string("latex.") + c, f, l};
}

std::vector<Palette> build() {
    std::vector<Palette> p;

    /* ---------------- symbol palettes ---------------- */

    p.push_back({u8"関係演算子", u8"≤≠", 3, {
        sym("\\leq", u8"≤", u8"以下"),
        sym("\\geq", u8"≥", u8"以上"),
        sym("\\ll", u8"≪", u8"はるかに小さい"),
        sym("\\gg", u8"≫", u8"はるかに大きい"),
        sym("\\prec", u8"≺", u8"順序が先"),
        sym("\\succ", u8"≻", u8"順序が後"),
        sym("\\triangleleft", u8"◁", u8"左向き三角"),
        sym("\\triangleright", u8"▷", u8"右向き三角"),
        sym("\\neq", u8"≠", u8"等しくない"),
        sym("\\equiv", u8"≡", u8"恒等・合同"),
        sym("\\approx", u8"≈", u8"ほぼ等しい"),
        sym("\\cong", u8"≅", u8"合同"),
        sym("\\sim", u8"∼", u8"相似"),
        sym("\\simeq", u8"≃", u8"漸近的に等しい"),
        sym("\\propto", u8"∝", u8"比例"),
        sym("\\sqsubseteq", u8"⊑", u8"角付き部分集合"),
        sym("\\sqsupseteq", u8"⊒", u8"角付き上位集合"),
        sym("\\parallel", u8"∥", u8"平行"),
        sym("\\Vert", u8"‖", u8"二重縦線（ノルム）"),
        sym("\\mid", u8"∣", u8"割り切る"),
        sym("\\vdash", u8"⊢", u8"導出される"),
        sym("\\models", u8"⊨", u8"充足する"),
        sym("\\frown", u8"⌢", u8"フラウン（記号）"),
        sym("\\smile", u8"⌣", u8"スマイル（記号）"),
    }});

    p.push_back({u8"空白と点", u8"⋯", 3, {
        raw("\\!", u8"−", u8"負の空き"),
        raw("\\,", u8"␣", u8"細い空き（1/6 em）"),
        raw("\\:", u8"␣␣", u8"中位の空き（2/9 em）"),
        raw("\\;", u8"␣␣␣", u8"広い空き（5/18 em）"),
        raw("\\quad", u8"1em", u8"1 em の空き"),
        raw("\\qquad", u8"2em", u8"2 em の空き"),
        sym("\\ldots", u8"…", u8"下寄りの省略記号"),
        sym("\\cdots", u8"⋯", u8"中央の省略記号"),
        sym("\\vdots", u8"⋮", u8"縦の省略記号"),
        sym("\\ddots", u8"⋱", u8"斜めの省略記号"),
        sym("\\therefore", u8"∴", u8"ゆえに"),
        sym("\\because", u8"∵", u8"なぜならば"),
    }});

    p.push_back({u8"装飾", u8"â", 3, {
        tpl("hat", u8"x̂", u8"ハット"),
        tpl("tilde", u8"x̃", u8"チルダ"),
        tpl("bar", u8"x̄", u8"バー（1文字）"),
        tpl("vec", u8"x⃗", u8"ベクトル矢印"),
        tpl("dot", u8"ẋ", u8"ドット"),
        tpl("ddot", u8"ẍ", u8"二重ドット"),
        tpl("dddot", u8"x⃛", u8"三重ドット"),
        tpl("prime", u8"x′", u8"プライム"),
        tpl("dprime", u8"x″", u8"二重プライム"),
        tpl("tprime", u8"x‴", u8"三重プライム"),
        tpl("strike", u8"x̸", u8"打ち消し線"),
        tpl("frown", u8"x⌢", u8"フラウン"),
        tpl("smile", u8"x⌣", u8"スマイル"),
    }});

    p.push_back({u8"演算子", u8"±×", 3, {
        sym("\\pm", u8"±", u8"プラスマイナス"),
        sym("\\mp", u8"∓", u8"マイナスプラス"),
        sym("\\times", u8"×", u8"乗算"),
        sym("\\div", u8"÷", u8"除算"),
        sym("\\ast", u8"∗", u8"アスタリスク"),
        sym("\\cdot", u8"⋅", u8"中黒（積）"),
        sym("\\star", u8"⋆", u8"星"),
        sym("\\bullet", u8"•", u8"黒丸"),
        sym("\\otimes", u8"⊗", u8"テンソル積"),
        sym("\\oplus", u8"⊕", u8"直和"),
        sym("\\ominus", u8"⊖", u8"丸マイナス"),
        sym("\\odot", u8"⊙", u8"丸ドット"),
        sym("\\diamond", u8"⋄", u8"ダイヤ演算子"),
        sym("\\bigtriangledown", u8"▽", u8"下向き三角"),
        sym("\\dagger", u8"†", u8"ダガー"),
        sym("\\ddagger", u8"‡", u8"二重ダガー"),
        sym("\\dag", u8"†", u8"ダガー（別名）"),
        sym("\\setminus", u8"∖", u8"差集合"),
        sym("\\backslash", u8"\\", u8"バックスラッシュ"),
        sym("\\surd", u8"√", u8"根号記号"),
    }});

    p.push_back({u8"矢印", u8"→⇒", 3, {
        sym("\\rightarrow", u8"→", u8"右矢印"),
        sym("\\leftarrow", u8"←", u8"左矢印"),
        sym("\\leftrightarrow", u8"↔", u8"両向き矢印"),
        sym("\\uparrow", u8"↑", u8"上矢印"),
        sym("\\downarrow", u8"↓", u8"下矢印"),
        sym("\\updownarrow", u8"↕", u8"上下矢印"),
        sym("\\Rightarrow", u8"⇒", u8"二重右矢印"),
        sym("\\Leftarrow", u8"⇐", u8"二重左矢印"),
        sym("\\Leftrightarrow", u8"⇔", u8"二重両向き矢印"),
        sym("\\Uparrow", u8"⇑", u8"二重上矢印"),
        sym("\\to", u8"→", u8"右矢印（\\to）"),
        sym("\\mapsto", u8"↦", u8"写像"),
        sym("\\Longrightarrow", u8"⟹", u8"長い二重右矢印"),
        sym("\\Longleftarrow", u8"⟸", u8"長い二重左矢印"),
        sym("\\Longleftrightarrow", u8"⟺", u8"長い二重両向き矢印"),
        sym("\\hookrightarrow", u8"↪", u8"右フック矢印"),
        sym("\\hookleftarrow", u8"↩", u8"左フック矢印"),
        sym("\\nearrow", u8"↗", u8"右上矢印"),
        sym("\\searrow", u8"↘", u8"右下矢印"),
        sym("\\swarrow", u8"↙", u8"左下矢印"),
        sym("\\nwarrow", u8"↖", u8"左上矢印"),
        sym("\\rightharpoonup", u8"⇀", u8"右上ハープーン"),
        sym("\\rightharpoondown", u8"⇁", u8"右下ハープーン"),
        sym("\\leftharpoonup", u8"↼", u8"左上ハープーン"),
        sym("\\leftharpoondown", u8"↽", u8"左下ハープーン"),
    }});

    p.push_back({u8"論理記号", u8"∀∧", 2, {
        sym("\\forall", u8"∀", u8"すべての"),
        sym("\\exists", u8"∃", u8"存在する"),
        sym("\\nexists", u8"∄", u8"存在しない"),
        sym("\\neg", u8"¬", u8"否定"),
        sym("\\wedge", u8"∧", u8"論理積"),
        sym("\\vee", u8"∨", u8"論理和"),
        sym("\\top", u8"⊤", u8"真"),
        sym("\\perp", u8"⊥", u8"垂直・偽"),
        sym("\\not", u8"̸", u8"打ち消し"),
    }});

    p.push_back({u8"集合記号", u8"∈∪", 3, {
        sym("\\in", u8"∈", u8"要素である"),
        sym("\\notin", u8"∉", u8"要素でない"),
        sym("\\ni", u8"∋", u8"要素として含む"),
        sym("\\subset", u8"⊂", u8"部分集合"),
        sym("\\supset", u8"⊃", u8"上位集合"),
        sym("\\not\\subset", u8"⊄", u8"部分集合でない"),
        sym("\\subseteq", u8"⊆", u8"部分集合または等しい"),
        sym("\\supseteq", u8"⊇", u8"上位集合または等しい"),
        sym("\\emptyset", u8"∅", u8"空集合"),
        sym("\\cup", u8"∪", u8"和集合"),
        sym("\\cap", u8"∩", u8"共通部分"),
        sym("\\sqcup", u8"⊔", u8"角付き和"),
        sym("\\sqcap", u8"⊓", u8"角付き共通部分"),
        sym("\\bigcup", u8"⋃", u8"大きい和集合"),
        sym("\\bigcap", u8"⋂", u8"大きい共通部分"),
    }});

    /* Ordered by what this lab actually reaches for, not by Eqnedt32's
     * arrangement.  The first rows are the working symbols of
     * electromagnetics and analysis; the curiosities keep their place at the
     * end so pasted TeX still parses, but they no longer sit above the
     * things people want.  \hbar was missing outright, which is a strange
     * gap for a Planck constant. */
    p.push_back({u8"その他の記号", u8"∂∞", 3, {
        sym("\\partial", u8"∂", u8"偏微分"),
        sym("\\nabla", u8"∇", u8"ナブラ"),
        sym("\\infty", u8"∞", u8"無限大"),
        sym("\\Re", u8"ℜ", u8"実部"),
        sym("\\Im", u8"ℑ", u8"虚部"),
        sym("\\hbar", u8"ℏ", u8"換算プランク定数"),
        sym("\\angle", u8"∠", u8"角"),
        sym("\\circ", u8"∘", u8"合成・度（^\\circ）"),
        sym("\\mho", u8"℧", u8"モー（コンダクタンス）"),
        sym("\\ell", u8"ℓ", u8"筆記体の l"),
        sym("\\prime", u8"′", u8"プライム記号"),
        sym("\\coprod", u8"∐", u8"余積記号"),
        /* The bare delimiter glyphs.  The 括弧 palette inserts matched,
         * stretching pairs; these are the single characters, for the times a
         * half-open interval or a lone bracket is what is wanted. */
        sym("\\langle", u8"⟨", u8"左山括弧（単体）"),
        sym("\\rangle", u8"⟩", u8"右山括弧（単体）"),
        sym("\\lfloor", u8"⌊", u8"左床記号（単体）"),
        sym("\\rfloor", u8"⌋", u8"右床記号（単体）"),
        sym("\\lceil", u8"⌈", u8"左天井記号（単体）"),
        sym("\\rceil", u8"⌉", u8"右天井記号（単体）"),
        /* Rarely wanted here, kept so pasted TeX still reads. */
    }});

    p.push_back({u8"ギリシャ小文字", u8"αβ", 4, {
        sym("\\alpha", u8"α", u8"アルファ"),
        sym("\\beta", u8"β", u8"ベータ"),
        sym("\\gamma", u8"γ", u8"ガンマ"),
        sym("\\delta", u8"δ", u8"デルタ"),
        sym("\\epsilon", u8"ϵ", u8"イプシロン"),
        sym("\\varepsilon", u8"ε", u8"イプシロン（異体）"),
        sym("\\zeta", u8"ζ", u8"ゼータ"),
        sym("\\eta", u8"η", u8"エータ"),
        sym("\\theta", u8"θ", u8"シータ"),
        sym("\\vartheta", u8"ϑ", u8"シータ（異体）"),
        sym("\\iota", u8"ι", u8"イオタ"),
        sym("\\kappa", u8"κ", u8"カッパ"),
        sym("\\varkappa", u8"ϰ", u8"カッパ（異体）"),
        sym("\\lambda", u8"λ", u8"ラムダ"),
        sym("\\mu", u8"μ", u8"ミュー"),
        sym("\\nu", u8"ν", u8"ニュー"),
        sym("\\xi", u8"ξ", u8"クサイ"),
        sym("\\pi", u8"π", u8"パイ"),
        sym("\\varpi", u8"ϖ", u8"パイ（異体）"),
        sym("\\rho", u8"ρ", u8"ロー"),
        sym("\\varrho", u8"ϱ", u8"ロー（異体）"),
        sym("\\sigma", u8"σ", u8"シグマ"),
        sym("\\varsigma", u8"ς", u8"シグマ（語末形）"),
        sym("\\tau", u8"τ", u8"タウ"),
        sym("\\upsilon", u8"υ", u8"ウプシロン"),
        sym("\\phi", u8"ϕ", u8"ファイ"),
        sym("\\varphi", u8"φ", u8"ファイ（異体）"),
        sym("\\chi", u8"χ", u8"カイ"),
        sym("\\psi", u8"ψ", u8"プサイ"),
        sym("\\omega", u8"ω", u8"オメガ"),
    }});

    p.push_back({u8"ギリシャ大文字", u8"ΓΩ", 4, {
        sym("\\Gamma", u8"Γ", u8"ガンマ"),
        sym("\\Delta", u8"Δ", u8"デルタ"),
        sym("\\Theta", u8"Θ", u8"シータ"),
        sym("\\Lambda", u8"Λ", u8"ラムダ"),
        sym("\\Xi", u8"Ξ", u8"クサイ"),
        sym("\\Pi", u8"Π", u8"パイ"),
        sym("\\Sigma", u8"Σ", u8"シグマ"),
        sym("\\Upsilon", u8"Υ", u8"ウプシロン"),
        sym("\\Phi", u8"Φ", u8"ファイ"),
        sym("\\Psi", u8"Ψ", u8"プサイ"),
        sym("\\Omega", u8"Ω", u8"オメガ"),
    }});

    /* ---------------- template palettes ---------------- */

    p.push_back({u8"括弧", u8"(□)", 3, {
        tpl("paren", u8"(▯)", u8"丸括弧"),
        tpl("bracket", u8"[▯]", u8"角括弧"),
        tpl("brace", u8"{▯}", u8"波括弧"),
        tpl("angle", u8"⟨▯⟩", u8"山括弧"),
        tpl("abs", u8"|▯|", u8"絶対値"),
        tpl("norm", u8"‖▯‖", u8"ノルム"),
        tpl("floor", u8"⌊▯⌋", u8"床関数"),
        tpl("ceil", u8"⌈▯⌉", u8"天井関数"),
        tpl("dirac", u8"⟨▯|▯⟩", u8"ブラケット"),
    }});

    p.push_back({u8"分数と根号", u8"½√", 2, {
        tpl("frac", u8"▯/▯", u8"分数"),
        tpl("slashfrac", u8"▯⁄▯", u8"スラッシュ分数"),
        tpl("sqrt", u8"√▯", u8"平方根"),
        tpl("nthroot", u8"ⁿ√▯", u8"n 乗根"),
    }});

    p.push_back({u8"上下付き", u8"x²", 3, {
        tpl("sup", u8"▯▫", u8"上付き"),
        tpl("sub", u8"▯₅", u8"下付き"),
        tpl("subsup", u8"▯₅▫", u8"上下付き"),
        tpl("over", u8"▫̅▯", u8"上側の注記"),
        tpl("under", u8"▯̲▫", u8"下側の注記"),
        tpl("lim", u8"lim", u8"極限（条件は下）"),
    }});

    p.push_back({u8"総和", u8"∑", 2, {
        tpl("sum", u8"∑▯", u8"総和（上下限つき）"),
        tpl("prod", u8"∏▯", u8"総乗（上下限つき）"),
    }});

    p.push_back({u8"積分", u8"∫", 3, {
        tpl("int", u8"∫", u8"積分"),
        tpl("iint", u8"∬", u8"二重積分"),
        tpl("iiint", u8"∭", u8"三重積分"),
        tpl("oint", u8"∮", u8"周回積分"),
        /* No \oiint or \oiiint: one needs esint and the other exists only
         * in packages that replace the whole math font, so either would give
         * the reader TeX their paper cannot typeset.  Pasted ones still
         * display; the editor simply never originates them. */
    }});

    p.push_back({u8"上線と下線", u8"¯_", 2, {
        tpl("overline", u8"▯̅", u8"上線（伸縮）"),
        tpl("underline", u8"▯̲", u8"下線（伸縮）"),
        tpl("overbrace", u8"⏞▯", u8"上の水平波括弧"),
        tpl("underbrace", u8"⏟▯", u8"下の水平波括弧"),
        tpl("overrightarrow", u8"▯⃗", u8"上の右向き伸縮矢印"),
        tpl("overleftarrow", u8"▯⃖", u8"上の左向き伸縮矢印"),
        tpl("overleftrightarrow", u8"▯⃡", u8"上の両向き伸縮矢印"),
    }});

    p.push_back({u8"総乗と集合演算", u8"∏∐", 2, {
        tpl("coprod", u8"∐▯", u8"余積（上下限つき）"),
        tpl("bigcup", u8"⋃▯", u8"大きい和集合（上下限つき）"),
        tpl("bigcap", u8"⋂▯", u8"大きい共通部分（上下限つき）"),
    }});

    p.push_back({u8"行列", u8"⊞", 4, {
        tpl("matrix1x2", u8"1×2", u8"1×2"),
        tpl("matrix1x3", u8"1×3", u8"1×3"),
        tpl("matrix2x1", u8"2×1", u8"2×1"),
        tpl("matrix2x2", u8"2×2", u8"2×2"),
        tpl("matrix2x3", u8"2×3", u8"2×3"),
        tpl("matrix3x1", u8"3×1", u8"3×1"),
        tpl("matrix3x2", u8"3×2", u8"3×2"),
        tpl("matrix3x3", u8"3×3", u8"3×3"),
        tpl("matrix4x4", u8"4×4", u8"4×4"),
        tpl("matrix5x5", u8"5×5", u8"5×5"),
        tpl("matrix6x6", u8"6×6", u8"6×6"),
        I{"matrix.add_row", u8"＋行", u8"行を下へ追加"},
        I{"matrix.remove_row", u8"－行", u8"現在行を削除"},
        I{"matrix.add_column", u8"＋列", u8"列を右へ追加"},
        I{"matrix.remove_column", u8"－列", u8"現在列を削除"},
        tpl("cases", u8"{▯", u8"場合分け"),
    }});

    return p;
}

}  // namespace

const std::vector<Palette>& palettes() {
    static const std::vector<Palette> kPalettes = build();
    return kPalettes;
}

int symbol_palette_count() { return 10; }

const std::vector<PaletteCategory>& palette_categories() {
    /* Keep the five concepts and their order aligned with the web editor.
     * The native catalogue is larger, so each tab contains its existing
     * compact popup palettes rather than duplicating individual commands. */
    static const std::vector<PaletteCategory> kCategories = {
        /* Follow the web editor's learning order: structures first, then
         * brackets, then decoration.  Matrix is a structure, not geometry. */
        {u8"基本", {11, 12, 17, 10, 2}},
        {u8"解析", {0, 3, 7, 13, 14}},
        {u8"集合・記号", {4, 6, 5, 16, 1}},
        {u8"幾何", {15}},
        {u8"ギリシャ", {8, 9}},
    };
    return kCategories;
}

}  // namespace eqnedit
