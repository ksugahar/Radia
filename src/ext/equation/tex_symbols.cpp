/* SPDX-License-Identifier: BSD-2-Clause */
#include "tex_symbols.h"

#include <cstdint>
#include <cstdlib>
#include <cstring>

typedef struct { const char *cmd; uint16_t code; } CmdEntry;

static const CmdEntry LATEX_TO_UNICODE[] = {
    {"\\Delta", 0x0394}, {"\\Gamma", 0x0393}, {"\\Lambda", 0x039B},
    {"\\Leftarrow", 0x21D0}, {"\\Leftrightarrow", 0x21D4},
    {"\\Longleftarrow", 0x27F8}, {"\\Longleftrightarrow", 0x27FA},
    {"\\Longrightarrow", 0x27F9},
    {"\\Omega", 0x03A9}, {"\\Phi", 0x03A6}, {"\\Pi", 0x03A0},
    {"\\Psi", 0x03A8}, {"\\Rightarrow", 0x21D2},
    {"\\Sigma", 0x03A3}, {"\\Theta", 0x0398},
    {"\\Uparrow", 0x21D1}, {"\\Upsilon", 0x03A5}, {"\\Xi", 0x039E},
    {"\\aleph", 0x2135}, {"\\alpha", 0x03B1}, {"\\approx", 0x2248},
    {"\\ast", 0x2217},
    {"\\because", 0x2235}, {"\\beta", 0x03B2}, {"\\bigcap", 0x22C2},
    {"\\bigcup", 0x22C3}, {"\\bigtriangledown", 0x25BD},
    {"\\bullet", 0x2022},
    {"\\cap", 0x2229}, {"\\cdot", 0x22C5}, {"\\cdots", 0x2026},  /* EQNEDT32 uses U+2026 (HORIZONTAL ELLIPSIS), not U+22EF */
    {"\\chi", 0x03C7}, {"\\clubsuit", 0x2663}, {"\\cong", 0x2245},
    {"\\coprod", 0x2210}, {"\\cup", 0x222A},
    {"\\dag", 0x2020}, {"\\dagger", 0x2020}, {"\\ddagger", 0x2021},
    {"\\ddots", 0x22F1},
    {"\\delta", 0x03B4}, {"\\diamond", 0x22C4},
    {"\\diamondsuit", 0x2666}, {"\\div", 0x00F7},
    {"\\downarrow", 0x2193},
    {"\\ell", 0x2113}, {"\\emptyset", 0x2205},
    {"\\epsilon", 0x03B5}, {"\\equiv", 0x2261}, {"\\eta", 0x03B7},
    {"\\exists", 0x2203},
    {"\\forall", 0x2200}, {"\\frown", 0x2322},
    {"\\gamma", 0x03B3}, {"\\geq", 0x2265}, {"\\gg", 0x226B},
    {"\\heartsuit", 0x2665},
    /* Sorted: bsearch walks past anything out of order, and these two were
     * swapped, so \hookleftarrow was in the table and unreachable. */
    {"\\hookleftarrow", 0x21A9}, {"\\hookrightarrow", 0x21AA},
    {"\\imath", 0x0131}, {"\\in", 0x2208}, {"\\infty", 0x221E},
    {"\\iota", 0x03B9},
    {"\\kappa", 0x03BA},
    {"\\lambda", 0x03BB}, {"\\langle", 0x2329},
    {"\\lceil", 0x2308}, {"\\ldots", 0x2026},
    {"\\leftarrow", 0x2190}, {"\\leftharpoondown", 0x21BD},
    {"\\leftharpoonup", 0x21BC}, {"\\leftrightarrow", 0x2194},
    {"\\leq", 0x2264}, {"\\lfloor", 0x230A}, {"\\ll", 0x226A},
    {"\\mapsto", 0x21A6}, {"\\mid", 0x2223}, {"\\models", 0x22A8},
    {"\\mp", 0x2213}, {"\\mu", 0x03BC},
    {"\\nabla", 0x2207}, {"\\nearrow", 0x2197}, {"\\neg", 0x00AC},
    {"\\neq", 0x2260}, {"\\nexists", 0x2204}, {"\\ni", 0x220B},
    {"\\not", 0x0338}, {"\\notin", 0x2209}, {"\\nu", 0x03BD},
    {"\\nwarrow", 0x2196},
    {"\\odot", 0x2299}, {"\\omega", 0x03C9}, {"\\ominus", 0x2296},
    {"\\oplus", 0x2295}, {"\\otimes", 0x2297},
    {"\\parallel", 0x2225}, {"\\partial", 0x2202},
    {"\\perp", 0x22A5}, {"\\phi", 0x03C6}, {"\\pi", 0x03C0},
    {"\\pm", 0x00B1}, {"\\prec", 0x227A}, {"\\prime", 0x2032},
    {"\\propto", 0x221D}, {"\\psi", 0x03C8},
    {"\\rangle", 0x232A}, {"\\rceil", 0x2309},
    {"\\rfloor", 0x230B}, {"\\rho", 0x03C1},
    {"\\rightarrow", 0x2192}, {"\\rightharpoondown", 0x21C1},
    {"\\rightharpoonup", 0x21C0},
    {"\\searrow", 0x2198}, {"\\setminus", 0x2216},
    {"\\sigma", 0x03C3}, {"\\sim", 0x223C}, {"\\simeq", 0x2243},
    {"\\smile", 0x2323}, {"\\spadesuit", 0x2660},
    {"\\sqcap", 0x2293}, {"\\sqcup", 0x2294},
    {"\\sqsubseteq", 0x2291}, {"\\sqsupseteq", 0x2292},
    {"\\star", 0x22C6}, {"\\subset", 0x2282},
    {"\\subseteq", 0x2286}, {"\\succ", 0x227B},
    {"\\supset", 0x2283}, {"\\supseteq", 0x2287},
    {"\\surd", 0x221A}, {"\\swarrow", 0x2199},
    {"\\tau", 0x03C4}, {"\\therefore", 0x2234},
    {"\\theta", 0x03B8}, {"\\times", 0x00D7}, {"\\to", 0x2192},
    {"\\top", 0x22A4}, {"\\triangleleft", 0x25C1},
    {"\\triangleright", 0x25B7},
    {"\\uparrow", 0x2191}, {"\\updownarrow", 0x2195},
    {"\\upsilon", 0x03C5},
    {"\\varepsilon", 0x03B5},
    /* \var* Greek: U+03D1/03D5/03D6/03F0/03F1 are "SYMBOL" codepoints absent in EQNEDT32 fonts.
     * Map each to its base-letter equivalent so the glyph appears. */
    {"\\varkappa", 0x03BA},  /* ϰ U+03F0 → κ U+03BA */
    {"\\varphi",   0x03D5},  /* ϕ U+03D5 (EQNEDT32 renders correctly) */
    {"\\varpi",    0x03C0},  /* ϖ U+03D6 → π U+03C0 */
    {"\\varrho",   0x03C1},  /* ϱ U+03F1 → ρ U+03C1 */
    {"\\varsigma", 0x03C2},  /* ς U+03C2 (standard Greek block, fine as-is) */
    {"\\vartheta", 0x03B8},  /* ϑ U+03D1 → θ U+03B8 */
    {"\\vdash", 0x22A2}, {"\\vdots", 0x22EE}, {"\\vee", 0x2228},
    {"\\wedge", 0x2227}, {"\\wp", 0x2118},
    {"\\xi", 0x03BE},
    {"\\zeta", 0x03B6},
};
#define LATEX_TO_UNICODE_N (sizeof(LATEX_TO_UNICODE)/sizeof(LATEX_TO_UNICODE[0]))

static int cmp_cmd(const void *a, const void *b) {
    return strcmp(((const CmdEntry *)a)->cmd, ((const CmdEntry *)b)->cmd);
}

static int lookup_latex_unicode(const char *cmd) {
    CmdEntry key = { cmd, 0 };
    CmdEntry *found = (CmdEntry *)bsearch(&key, LATEX_TO_UNICODE,
        LATEX_TO_UNICODE_N, sizeof(CmdEntry), cmp_cmd);
    return found ? (int)found->code : -1;
}

int tex_command_to_unicode(const char *cmd) {
    return lookup_latex_unicode(cmd);
}

int tex_command_count(void) {
    return (int)LATEX_TO_UNICODE_N;
}

const char *tex_command_name(int index) {
    if (index < 0 || index >= (int)LATEX_TO_UNICODE_N) return 0;
    return LATEX_TO_UNICODE[index].cmd;
}

int tex_command_code_at(int index) {
    if (index < 0 || index >= (int)LATEX_TO_UNICODE_N) return -1;
    return (int)LATEX_TO_UNICODE[index].code;
}
