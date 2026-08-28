#include "math_symbols.h"

#include "equation_types.h"

#include <utility>

namespace eqnedit {
namespace {

using Entry = std::pair<const char*, int>;
constexpr Entry kSymbols[] = {
    {"\\backslash", 0x005C},
    {"\\Delta", 0x0394}, {"\\Gamma", 0x0393}, {"\\Lambda", 0x039B},
    {"\\Leftarrow", 0x21D0}, {"\\Leftrightarrow", 0x21D4},
    {"\\Longleftarrow", 0x27F8}, {"\\Longleftrightarrow", 0x27FA},
    {"\\Longrightarrow", 0x27F9}, {"\\Omega", 0x03A9},
    {"\\Phi", 0x03A6}, {"\\Pi", 0x03A0}, {"\\Psi", 0x03A8},
    {"\\Rightarrow", 0x21D2}, {"\\Sigma", 0x03A3},
    {"\\Theta", 0x0398}, {"\\Uparrow", 0x21D1},
    {"\\Upsilon", 0x03A5}, {"\\Xi", 0x039E},
    /* No \aleph, \wp, \imath -- dropped with the card suits at the user's
     * request; nothing in this lab's notation uses them. */
    {"\\alpha", 0x03B1},
    {"\\approx", 0x2248}, {"\\ast", 0x2217},
    {"\\because", 0x2235}, {"\\beta", 0x03B2},
    {"\\bigcap", 0x22C2}, {"\\bigcup", 0x22C3},
    {"\\bigtriangledown", 0x25BD}, {"\\bullet", 0x2022},
    /* \cdots is the centred ellipsis U+22EF; U+2026 is \ldots on the
     * baseline.  Sharing one code point made \ldots come back as \cdots. */
    {"\\cap", 0x2229}, {"\\cdot", 0x22C5}, {"\\cdots", 0x22EF},
    {"\\chi", 0x03C7},
    {"\\cong", 0x2245}, {"\\coprod", 0x2210}, {"\\cup", 0x222A},
    {"\\dag", 0x2020}, {"\\dagger", 0x2020},
    {"\\ddagger", 0x2021}, {"\\ddots", 0x22F1},
    {"\\delta", 0x03B4}, {"\\diamond", 0x22C4},
    {"\\div", 0x00F7},
    {"\\downarrow", 0x2193}, {"\\ell", 0x2113},
    /* LaTeX's \epsilon is the lunate symbol and \varepsilon the ordinary
     * letter; the same inversion holds for \phi/\varphi.  Six var- commands
     * shared a code point with their plain form, so the canvas drew the
     * wrong glyph and Office received it. */
    {"\\emptyset", 0x2205}, {"\\epsilon", 0x03F5},
    {"\\equiv", 0x2261}, {"\\eta", 0x03B7}, {"\\exists", 0x2203},
    {"\\forall", 0x2200}, {"\\frown", 0x2322},
    {"\\gamma", 0x03B3}, {"\\geq", 0x2265}, {"\\gg", 0x226B},
    /* No card suits.  They are not mathematics, the user asked them gone,
     * and a pasted \heartsuit degrading to its letters is an acceptable
     * price for a palette that holds only working symbols. */
    {"\\hookrightarrow", 0x21AA},
    {"\\hookleftarrow", 0x21A9},
    {"\\in", 0x2208}, {"\\infty", 0x221E}, {"\\iota", 0x03B9},
    {"\\kappa", 0x03BA}, {"\\lambda", 0x03BB},
    {"\\langle", 0x27E8}, {"\\lceil", 0x2308},
    {"\\ldots", 0x2026}, {"\\leftarrow", 0x2190},
    {"\\leftharpoondown", 0x21BD}, {"\\leftharpoonup", 0x21BC},
    {"\\leftrightarrow", 0x2194}, {"\\leq", 0x2264},
    {"\\lfloor", 0x230A}, {"\\ll", 0x226A}, {"\\mapsto", 0x21A6},
    {"\\mid", 0x2223}, {"\\models", 0x22A8}, {"\\mp", 0x2213},
    {"\\mu", 0x03BC}, {"\\nabla", 0x2207}, {"\\nearrow", 0x2197},
    {"\\neg", 0x00AC}, {"\\neq", 0x2260}, {"\\nexists", 0x2204},
    {"\\ni", 0x220B}, {"\\not", 0x0338}, {"\\notin", 0x2209},
    /* LaTeX has no \nsubset, in any package -- it was emitted for years and
     * would not compile.  The way to spell U+2284 is \not applied to a
     * relation, and the parser folds the pair back into one character. */
    {"\\not\\subset", 0x2284},
    {"\\nu", 0x03BD}, {"\\nwarrow", 0x2196}, {"\\odot", 0x2299},
    {"\\omega", 0x03C9}, {"\\ominus", 0x2296},
    {"\\oplus", 0x2295}, {"\\otimes", 0x2297},
    /* The double bar of a norm.  Distinct from \parallel (a relation) and
     * from a single | (absolute value): without it, \|x\| outside
     * \left/\right read back as |x|, quietly changing a norm into an
     * absolute value. */
    {"\\Vert", 0x2016},
    /* Physics and engineering symbols the lab actually reaches for.  \angle
     * and \circ could already be WRITTEN by the emitter but not read back --
     * the same one-way asymmetry \overbrace had.  \hbar was in neither, which
     * is a strange omission for a Planck constant. */
    {"\\hbar", 0x210F}, {"\\angle", 0x2220}, {"\\circ", 0x2218},
    {"\\mho", 0x2127},
    {"\\parallel", 0x2225}, {"\\partial", 0x2202},
    {"\\perp", 0x22A5}, {"\\phi", 0x03D5}, {"\\pi", 0x03C0},
    {"\\pm", 0x00B1}, {"\\prec", 0x227A}, {"\\prime", 0x2032},
    {"\\propto", 0x221D}, {"\\psi", 0x03C8}, {"\\rangle", 0x27E9},
    /* Fraktur real and imaginary parts.  These are glyphs, not operator
     * names, so they must not be treated as \sin-style function words. */
    {"\\Re", 0x211C}, {"\\Im", 0x2111},
    {"\\rceil", 0x2309}, {"\\rfloor", 0x230B}, {"\\rho", 0x03C1},
    {"\\rightarrow", 0x2192}, {"\\rightharpoondown", 0x21C1},
    {"\\rightharpoonup", 0x21C0}, {"\\searrow", 0x2198},
    {"\\setminus", 0x2216}, {"\\sigma", 0x03C3}, {"\\sim", 0x223C},
    {"\\simeq", 0x2243}, {"\\smile", 0x2323},
    {"\\sqcap", 0x2293},
    {"\\sqcup", 0x2294}, {"\\sqsubseteq", 0x2291},
    {"\\sqsupseteq", 0x2292}, {"\\star", 0x22C6},
    {"\\subset", 0x2282}, {"\\subseteq", 0x2286},
    {"\\succ", 0x227B}, {"\\supset", 0x2283},
    {"\\supseteq", 0x2287}, {"\\surd", 0x221A},
    {"\\swarrow", 0x2199}, {"\\tau", 0x03C4},
    {"\\therefore", 0x2234}, {"\\theta", 0x03B8},
    {"\\times", 0x00D7}, {"\\to", 0x2192}, {"\\top", 0x22A4},
    {"\\triangleleft", 0x25C1}, {"\\triangleright", 0x25B7},
    {"\\uparrow", 0x2191}, {"\\updownarrow", 0x2195},
    {"\\upsilon", 0x03C5}, {"\\varepsilon", 0x03B5},
    {"\\varkappa", 0x03F0}, {"\\varphi", 0x03C6},
    {"\\varpi", 0x03D6}, {"\\varrho", 0x03F1},
    {"\\varsigma", 0x03C2}, {"\\vartheta", 0x03D1},
    {"\\vdash", 0x22A2}, {"\\vdots", 0x22EE}, {"\\vee", 0x2228},
    {"\\wedge", 0x2227}, {"\\xi", 0x03BE},
    {"\\zeta", 0x03B6},
};

}  // namespace

int latex_symbol_codepoint(const std::string& command) {
    for (const auto& item : kSymbols)
        if (command == item.first) return item.second;
    return -1;
}

bool is_latex_symbol_codepoint(uint32_t codepoint) {
    for (const auto& item : kSymbols)
        if (uint32_t(item.second) == codepoint) return true;
    return false;
}

int typeface_for_code(uint32_t codepoint) {
    if (codepoint >= 0x0391 && codepoint <= 0x03A9) return TF_UCGREEK;
    if (codepoint >= 0x03B1 && codepoint <= 0x03C9) return TF_LCGREEK;
    /* Variant lowercase Greek: theta, phi, pi, kappa, rho, epsilon symbols.
     * They sit outside the Greek block, so a range test alone left them
     * upright while every other lowercase Greek letter was italic. */
    switch (codepoint) {
    case 0x03D1: case 0x03D5: case 0x03D6:
    case 0x03F0: case 0x03F1: case 0x03F5:
        return TF_LCGREEK;
    default:
        return TF_SYMBOL;
    }
}

std::vector<std::string> latex_symbol_commands() {
    std::vector<std::string> commands;
    commands.reserve(sizeof(kSymbols) / sizeof(kSymbols[0]));
    for (const auto& item : kSymbols) commands.emplace_back(item.first);
    return commands;
}

}  // namespace eqnedit
