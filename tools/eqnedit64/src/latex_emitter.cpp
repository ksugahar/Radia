/* Structural equation tree -> normalized TeX text. */
#include "latex_emitter.h"
#include <algorithm>
#include <cctype>
#include <cstring>
#include <cstdio>

namespace eqnedit {

/* ============================================================
 * Unicode -> TeX symbol lookup tables
 * ============================================================ */
struct MapEntry { uint16_t key; const char* val; };

static const MapEntry SYMBOL_MAP[] = {
    {0x0021, "!"}, {0x0022, "\""}, {0x0023, "\\#"}, {0x0024, "\\$"},
    {0x0025, "\\%"},
    {0x0026, "\\&"}, {0x0028, "("}, {0x0029, ")"}, {0x002A, "*"},
    {0x002B, "+"}, {0x002C, ","}, {0x002D, "-"}, {0x002E, "."},
    {0x002F, "/"}, {0x003A, ":"}, {0x003B, ";"}, {0x003C, "<"},
    {0x003D, " = "}, {0x003E, ">"}, {0x003F, "?"}, {0x0040, "@"},
    {0x005B, "["}, {0x005C, "\\backslash "}, {0x005D, "]"},
    {0x005E, " \\wedge "}, {0x005F, "\\_"},
    {0x007B, "\\{"}, {0x007C, "|"}, {0x007D, "\\}"},
    /* An ASCII tilde typed as a glyph is not the U+223C relation below. */
    {0x007E, "\\textasciitilde{}"},
    {0x00AC, " \\neg "}, {0x00B0, "^{\\circ }"}, {0x00B1, " \\pm "},
    {0x00B7, " \\cdot "}, {0x00D7, " \\times "}, {0x00F7, " \\div "},
};
#define SYMBOL_MAP_N (sizeof(SYMBOL_MAP)/sizeof(SYMBOL_MAP[0]))

/* Sorted by key: map_lookup() binary-searches this table.  Every entry the
 * parser's symbol table can produce must appear here, otherwise emitChar()
 * falls through to a raw UTF-8 glyph and the saved .tex stops compiling.
 * tests/test_symbols.py checks that correspondence for the whole table. */
static const MapEntry UNICODE_MAP[] = {
    {0x0338, "\\not "}, {0x0393, "\\Gamma "},
    {0x0394, "\\Delta "}, {0x0398, "\\Theta "}, {0x039B, "\\Lambda "},
    {0x039E, "\\Xi "}, {0x03A0, "\\Pi "}, {0x03A3, "\\Sigma "},
    {0x03A5, "\\Upsilon "}, {0x03A6, "\\Phi "}, {0x03A8, "\\Psi "},
    {0x03A9, "\\Omega "}, {0x03B1, "\\alpha "}, {0x03B2, "\\beta "},
    {0x03B3, "\\gamma "}, {0x03B4, "\\delta "}, {0x03B5, "\\varepsilon "},
    {0x03B6, "\\zeta "}, {0x03B7, "\\eta "}, {0x03B8, "\\theta "},
    {0x03B9, "\\iota "}, {0x03BA, "\\kappa "}, {0x03BB, "\\lambda "},
    {0x03BC, "\\mu "}, {0x03BD, "\\nu "}, {0x03BE, "\\xi "},
    {0x03C0, "\\pi "}, {0x03C1, "\\rho "}, {0x03C2, "\\varsigma "},
    {0x03C3, "\\sigma "}, {0x03C4, "\\tau "}, {0x03C5, "\\upsilon "},
    {0x03C6, "\\varphi "}, {0x03C7, "\\chi "}, {0x03C8, "\\psi "},
    {0x03C9, "\\omega "}, {0x03D1, "\\vartheta "}, {0x03D5, "\\phi "},
    {0x03D6, "\\varpi "}, {0x03F0, "\\varkappa "}, {0x03F1, "\\varrho "},
    {0x03F5, "\\epsilon "}, {0x2016, "\\Vert "}, {0x2019, "'"},
    {0x2020, "\\dagger "}, {0x2021, "\\ddagger "}, {0x2022, "\\bullet "},
    {0x2026, "\\ldots "}, {0x2032, "'"}, {0x2033, "''"}, {0x2034, "'''"},
    {0x2102, "\\mathbb{C}"}, {0x210D, "\\mathbb{H}"},
    {0x210F, "\\hbar "}, {0x2111, "\\Im "},
    {0x2113, "\\ell "}, {0x2115, "\\mathbb{N}"},
    {0x2119, "\\mathbb{P}"}, {0x211A, "\\mathbb{Q}"}, {0x211C, "\\Re "},
    {0x211D, "\\mathbb{R}"}, {0x2124, "\\mathbb{Z}"},
    {0x2127, "\\mho "},
    {0x2190, "\\leftarrow "}, {0x2191, "\\uparrow "},
    {0x2192, "\\rightarrow "}, {0x2193, "\\downarrow "},
    {0x2194, "\\leftrightarrow "}, {0x2195, "\\updownarrow "},
    {0x2196, "\\nwarrow "}, {0x2197, "\\nearrow "}, {0x2198, "\\searrow "},
    {0x2199, "\\swarrow "}, {0x21A6, " \\mapsto "},
    {0x21A9, " \\hookleftarrow "}, {0x21AA, " \\hookrightarrow "},
    {0x21BC, " \\leftharpoonup "}, {0x21BD, " \\leftharpoondown "},
    {0x21C0, " \\rightharpoonup "}, {0x21C1, " \\rightharpoondown "},
    {0x21D0, "\\Leftarrow "}, {0x21D1, "\\Uparrow "},
    {0x21D2, "\\Rightarrow "}, {0x21D3, "\\Downarrow "},
    {0x21D4, "\\Leftrightarrow "}, {0x2200, " \\forall "},
    {0x2202, "\\partial "}, {0x2203, " \\exists "}, {0x2204, " \\nexists "},
    {0x2205, "\\emptyset "}, {0x2207, "\\nabla "}, {0x2208, " \\in "},
    {0x2209, " \\notin "}, {0x220B, " \\ni "}, {0x2210, " \\coprod "},
    {0x2211, "\\sum "}, {0x2212, " - "}, {0x2213, " \\mp "}, {0x2215, "/"},
    {0x2216, " \\setminus "}, {0x2217, "\\ast "}, {0x2218, " \\circ "},
    {0x221A, "\\surd "}, {0x221D, " \\propto "}, {0x221E, "\\infty "},
    {0x2220, "\\angle "}, {0x2223, " \\mid "}, {0x2225, "\\parallel "},
    {0x2227, " \\wedge "}, {0x2228, " \\vee "}, {0x2229, " \\cap "},
    {0x222A, " \\cup "}, {0x222B, "\\int "}, {0x222C, "\\iint "},
    {0x222D, "\\iiint "}, {0x222E, "\\oint "}, {0x222F, "\\oiint "},
    {0x2230, "\\oiiint "}, {0x2234, " \\therefore "},
    {0x2235, " \\because "}, {0x223C, " \\sim "}, {0x2243, " \\simeq "},
    {0x2245, " \\cong "}, {0x2248, " \\approx "}, {0x2260, " \\neq "},
    {0x2261, " \\equiv "}, {0x2264, " \\leq "}, {0x2265, " \\geq "},
    {0x226A, " \\ll "}, {0x226B, " \\gg "}, {0x227A, " \\prec "},
    {0x227B, " \\succ "}, {0x2282, " \\subset "}, {0x2283, " \\supset "},
    {0x2284, " \\not\\subset "}, {0x2286, " \\subseteq "},
    {0x2287, " \\supseteq "}, {0x2291, " \\sqsubseteq "},
    {0x2292, " \\sqsupseteq "}, {0x2293, " \\sqcap "},
    {0x2294, " \\sqcup "}, {0x2295, " \\oplus "}, {0x2296, " \\ominus "},
    {0x2297, " \\otimes "}, {0x2299, " \\odot "}, {0x22A2, " \\vdash "},
    {0x22A4, " \\top "}, {0x22A5, " \\perp "}, {0x22A8, " \\models "},
    {0x22C2, " \\bigcap "}, {0x22C3, " \\bigcup "}, {0x22C4, " \\diamond "},
    {0x22C5, " \\cdot "}, {0x22C6, " \\star "}, {0x22EE, "\\vdots "},
    {0x22EF, "\\cdots "}, {0x22F0, "\\iddots "}, {0x22F1, "\\ddots "},
    {0x2308, "\\lceil "}, {0x2309, "\\rceil "}, {0x230A, "\\lfloor "},
    {0x230B, "\\rfloor "}, {0x2322, " \\frown "}, {0x2323, " \\smile "},
    {0x2329, "\\langle "}, {0x232A, "\\rangle "},
    {0x25B7, " \\triangleright "}, {0x25BD, "\\bigtriangledown "},
    {0x25C1, " \\triangleleft "},
    {0x27E8, "\\langle "},
    {0x27E9, "\\rangle "}, {0x27F8, " \\Longleftarrow "},
    {0x27F9, " \\Longrightarrow "}, {0x27FA, " \\Longleftrightarrow "},
};
#define UNICODE_MAP_N (sizeof(UNICODE_MAP)/sizeof(UNICODE_MAP[0]))

static const char* map_lookup(const MapEntry* map, int n, uint16_t key) {
    int lo = 0, hi = n - 1;
    while (lo <= hi) {
        int mid = (lo + hi) / 2;
        if (map[mid].key == key) return map[mid].val;
        if (map[mid].key < key) lo = mid + 1;
        else hi = mid - 1;
    }
    return nullptr;
}

/* Embellishment format table */
struct EmbellFmt { const char* prefix; const char* suffix; };
static const EmbellFmt EMBELL_MAP[] = {
    /* 0 */ {nullptr, nullptr}, /* 1 */ {nullptr, nullptr},
    /* EM_DOT    2 */ {"\\dot{", "}"},
    /* EM_DDOT   3 */ {"\\ddot{", "}"},
    /* EM_TDOT   4 */ {"\\dddot{", "}"},
    /* EM_PRIME  5 */ {nullptr, "'"},
    /* EM_DPRIME 6 */ {nullptr, "''"},
    /* EM_BPRIME 7 */ {nullptr, "'''"},
    /* EM_TILDE  8 */ {"\\tilde{", "}"},
    /* EM_HAT    9 */ {"\\hat{", "}"},
    /* EM_NOT   10 */ {"\\not{", "}"},
    /* EM_RARROW 11 */ {"\\vec{", "}"},
    /* EM_LARROW 12 */ {"\\overleftarrow{", "}"},
    /* EM_BARROW 13 */ {"\\overleftrightarrow{", "}"},
    /* The spacing here matches what emitting \rightharpoonup as a symbol
     * produces, so loading and saving does not rewrite the file. */
    /* EM_R1ARROW 14 */ {"\\overset{ \\rightharpoonup }{", "}"},
    /* EM_L1ARROW 15 */ {"\\overset{ \\leftharpoonup }{", "}"},
    /* EM_MBAR  16 */ {"\\cancel{", "}"},
    /* EM_OBAR  17 */ {"\\bar{", "}"},
    /* EM_TPRIME 18 */ {nullptr, "'''"},
    /* EM_FROWN 19 */ {"\\overset{ \\frown }{", "}"},
    /* EM_SMILE 20 */ {"\\overset{ \\smile }{", "}"},
};
#define EMBELL_MAP_N (sizeof(EMBELL_MAP)/sizeof(EMBELL_MAP[0]))

/* ============================================================
 * Emit NodeList → string
 * ============================================================ */
std::string LaTeXEmitter::emitNodes(const NodeList& nodes) {
    std::string s;
    emitSequence(nodes, 0, nodes.size(), s);
    return s;
}

/* ============================================================
 * Main entry point
 * ============================================================ */
std::string LaTeXEmitter::emit(const LineNode& root) {
    std::string raw;
    emitLine(root, raw);
    return postProcess(raw);
}

static std::string utf8_of(uint32_t cp) {
    std::string s;
    if (cp < 0x80) s += char(cp);
    else if (cp < 0x800) {
        s += char(0xC0 | (cp >> 6)); s += char(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
        s += char(0xE0 | (cp >> 12)); s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xF0 | (cp >> 18)); s += char(0x80 | ((cp >> 12) & 0x3F));
        s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F));
    }
    return s;
}

static void append_text_char(std::string& out, uint32_t cp) {
    switch (cp) {
    case '\\': out += "\\textbackslash{}"; break;
    case '{':  out += "\\{"; break;
    case '}':  out += "\\}"; break;
    case '%':  out += "\\%"; break;
    case '#':  out += "\\#"; break;
    case '$':  out += "\\$"; break;
    case '_':  out += "\\_"; break;
    case '&':  out += "\\&"; break;
    case '^':  out += "\\textasciicircum{}"; break;
    case '~':  out += "\\textasciitilde{}"; break;
    default:   if (cp >= 0x20) out += utf8_of(cp); break;
    }
}

/* A math alphabet changes glyph style, not parsing mode.  Keep operators and
 * TeX symbol names mathematical inside \mathrm/\mathit; append_text_char()
 * is deliberately different because \text treats the same bytes as prose. */
static void append_math_char(std::string& out, uint32_t cp) {
    if (cp < 0x80) {
        switch (cp) {
        case '\\': out += "\\backslash "; break;
        case '{':  out += "\\{"; break;
        case '}':  out += "\\}"; break;
        case '%':  out += "\\%"; break;
        case '#':  out += "\\#"; break;
        case '$':  out += "\\$"; break;
        case '_':  out += "\\_"; break;
        case '&':  out += "\\&"; break;
        case '^':  out += "\\wedge "; break;
        case '~':  out += "\\sim "; break;
        default:   if (cp >= 0x20) out += char(cp); break;
        }
        return;
    }
    const char* symbol = cp <= 0xFFFF
        ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(cp)) : nullptr;
    if (!symbol) {
        out += utf8_of(cp);
        return;
    }
    /* Symbol-table entries include inter-atom spaces for the top-level
     * emitter (" \\leq ").  Inside a math alphabet wrapper those spaces
     * must not become stored content.  Retain only one terminator after a
     * control word so the next letter cannot turn \alpha b into \alphab. */
    std::string token(symbol);
    const size_t first = token.find_first_not_of(' ');
    if (first == std::string::npos) return;
    const size_t last = token.find_last_not_of(' ');
    token = token.substr(first, last - first + 1);
    out += token;
    if (!token.empty() &&
        std::isalpha(static_cast<unsigned char>(token.back())))
        out += ' ';
}

static const char* math_alphabet_command(int typeface) {
    switch (typeface) {
    case TF_ROMAN: return "\\mathrm";
    case TF_MATH_ITALIC: return "\\mathit";
    case TF_MATH_SANS: return "\\mathsf";
    case TF_MATH_MONO: return "\\mathtt";
    case TF_MATH_SCRIPT: return "\\mathcal";
    case TF_MATH_DOUBLE: return "\\mathbb";
    case TF_MATH_FRAKTUR: return "\\mathfrak";
    case TF_BOLD_SYMBOL: return "\\bm";
    default: return nullptr;
    }
}

std::string LaTeXEmitter::emit_range(const NodeList& nodes,
                                     size_t first, size_t last) {
    first = std::min(first, nodes.size());
    last = std::min(last, nodes.size());
    if (last < first) std::swap(first, last);
    std::string raw;
    emitSequence(nodes, first, last, raw);
    return postProcess(raw);
}

/* ============================================================
 * LINE emission (simplified — no 7-pass system yet)
 * ============================================================ */
void LaTeXEmitter::emitLine(const LineNode& line, std::string& out) {
    if (line.isNull) return;
    emitSequence(line.children, 0, line.children.size(), out);
}

void LaTeXEmitter::emitSequence(const NodeList& nodes, size_t first,
                                size_t last, std::string& out) {
    first = std::min(first, nodes.size());
    last = std::min(last, nodes.size());
    /* Output phase with text/function/math-alphabet grouping.
     * Consecutive TF_TEXT chars → \text{...}
     * Consecutive TF_FUNCTION chars → \sin, \cos, \operatorname{...}
     * Consecutive explicit Roman/italic chars → one readable wrapper. */
    std::string textBuf;
    std::string funcBuf;
    std::string alphabetBuf;
    int alphabetTypeface = -1;

    auto flushText = [&]() {
        if (!textBuf.empty()) {
            out += "\\text{"; out += textBuf; out += "}";
            textBuf.clear();
        }
    };
    auto flushFunc = [&]() {
        if (funcBuf.empty()) return;
        /* Known LaTeX function names */
        static const char* KNOWN_FUNCS[] = {
            "arccos","arcsin","arctan","arg","cos","cosh","cot","coth",
            "csc","deg","det","dim","exp","gcd","hom","inf",
            "ker","lg","lim","liminf","limsup","ln","log","max",
            "min","Pr","sec","sin","sinh","sup","tan","tanh",NULL
        };
        static const char* KNOWN_OPS[] = {
            "rot","curl","div","grad","sgn","tr","diag","mod",
            "Re","Im","Res","const",NULL
        };
        bool emitted = false;
        for (int i = 0; KNOWN_FUNCS[i]; ++i) {
            if (funcBuf == KNOWN_FUNCS[i]) {
                out += "\\"; out += KNOWN_FUNCS[i]; out += " ";
                emitted = true; break;
            }
        }
        if (!emitted) {
            for (int i = 0; KNOWN_OPS[i]; ++i) {
                if (funcBuf == KNOWN_OPS[i]) {
                    out += "\\operatorname{"; out += KNOWN_OPS[i]; out += "}";
                    emitted = true; break;
                }
            }
        }
        /* Explicit Function style must survive save/reopen even for a
         * user-defined name.  Raw letters would silently become variables. */
        if (!emitted) {
            out += "\\operatorname{"; out += funcBuf; out += "}";
        }
        funcBuf.clear();
    };
    auto flushAlphabet = [&]() {
        if (alphabetBuf.empty()) return;
        out += math_alphabet_command(alphabetTypeface);
        out += '{';
        out += alphabetBuf;
        out += '}';
        alphabetBuf.clear();
        alphabetTypeface = -1;
    };

    for (size_t i = first; i < last; ++i) {
        const auto& child = nodes[i];
        if (!child) continue;
        if (child->tag() == Node::kChar) {
            auto* ch = static_cast<const CharNode*>(child.get());
            if (ch->typeface == TF_TEXT) {
                flushFunc();
                flushAlphabet();
                append_text_char(textBuf, ch->charCode);
                continue;
            }
            if (ch->typeface == TF_FUNCTION) {
                flushText();
                flushAlphabet();
                append_text_char(funcBuf, ch->charCode);
                continue;
            }
            if (math_alphabet_command(ch->typeface) && ch->embells.empty()) {
                flushText();
                flushFunc();
                if (alphabetTypeface != ch->typeface) flushAlphabet();
                alphabetTypeface = ch->typeface;
                append_math_char(alphabetBuf, ch->charCode);
                continue;
            }
        }
        /* Non-text/func node: flush buffers and emit */
        flushText();
        flushFunc();
        flushAlphabet();
        emitNode(child.get(), out);
    }
    flushText();
    flushFunc();
    flushAlphabet();

}

/* ============================================================
 * Node dispatch
 * ============================================================ */
void LaTeXEmitter::emitNode(const Node* node, std::string& out) {
    if (!node) return;
    switch (node->tag()) {
    case Node::kLine:
        emitLine(*static_cast<const LineNode*>(node), out);
        break;
    case Node::kChar:
        emitChar(*static_cast<const CharNode*>(node), out);
        break;
    case Node::kFence:
        emitFence(*static_cast<const FenceNode*>(node), out);
        break;
    case Node::kFrac:
        emitFrac(*static_cast<const FracNode*>(node), out);
        break;
    case Node::kSqrt:
        emitSqrt(*static_cast<const SqrtNode*>(node), out);
        break;
    case Node::kScript:
        emitScript(*static_cast<const ScriptNode*>(node), out);
        break;
    case Node::kIntegral:
        emitIntegral(*static_cast<const IntegralNode*>(node), out);
        break;
    case Node::kBigOp:
        emitBigOp(*static_cast<const BigOpNode*>(node), out);
        break;
    case Node::kDecoration:
        emitDecoration(*static_cast<const DecorationNode*>(node), out);
        break;
    case Node::kBraceDeco:
        emitBraceDeco(*static_cast<const BraceDecoNode*>(node), out);
        break;
    case Node::kDirac:
        emitDirac(*static_cast<const DiracNode*>(node), out);
        break;
    case Node::kLim:
        emitLim(*static_cast<const LimNode*>(node), out);
        break;
    case Node::kOverset:
        emitOverset(*static_cast<const OversetNode*>(node), out);
        break;
    case Node::kPile:
        emitPile(*static_cast<const PileNode*>(node), out);
        break;
    case Node::kMatrix:
        emitMatrix(*static_cast<const MatrixNode*>(node), out);
        break;
    case Node::kSize:
        /* SIZE nodes are consumed by passes, not emitted directly */
        break;
    case Node::kEmbell:
        emitEmbell(*static_cast<const EmbellNode*>(node), out);
        break;
    case Node::kFont:
        /* FONT records ignored in output */
        break;
    default:
        break;
    }
}

/* ============================================================
 * CHAR emission
 * ============================================================ */
void LaTeXEmitter::emitChar(const CharNode& ch, std::string& out) {
    uint32_t code = ch.charCode;
    int tf = ch.typeface;

    /* Apply embellishments (prefix) */
    for (auto it = ch.embells.rbegin(); it != ch.embells.rend(); ++it) {
        int et = *it;
        if (et >= 0 && et < (int)EMBELL_MAP_N && EMBELL_MAP[et].prefix)
            out += EMBELL_MAP[et].prefix;
    }

    /* Explicit spacing is written back as the command it came from.  A
     * trailing space keeps `\quad x` from running together as `\quadx`. */
    if (tf == TF_SPACE) {
        out += ch.latex;
        out += ' ';
        return;
    }

    /* A character parsed from a named TeX command keeps that command.  The
     * code-point tables are a lossy inverse: several commands share one code
     * point (\epsilon/\varepsilon, \dag/\dagger, \ldots/\cdots), so a
     * round trip through them silently rewrites the author's source.  Only
     * the three table-driven typefaces are shortcut here; \mathbf{\alpha}
     * still has to go through the TF_VECTOR branch to keep its wrapper. */
    const bool namedSymbol = !ch.latex.empty() &&
        (tf == TF_SYMBOL || tf == TF_LCGREEK || tf == TF_UCGREEK);

    if (namedSymbol) {
        /* The tables carry inter-atom spacing as well as the name: " \in "
         * leads with a space, "\alpha " does not.  Keep that spacing and take
         * only the name from the node, otherwise a relation loses the space
         * in front of it the second time the line is serialized. */
        const char* spaced = code <= 0xFFFF ? lookupSymbol(uint16_t(code))
                                            : nullptr;
        if (spaced && spaced[0] == ' ') out += ' ';
        out += ch.latex;
        out += ' ';
    }
    /* Symbol typeface */
    else if (tf == TF_SYMBOL) {
        const char* s = code <= 0xFFFF
                      ? map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, uint16_t(code))
                      : nullptr;
        if (s) { out += s; }
        else {
            s = code <= 0xFFFF
              ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
              : nullptr;
            if (s) out += s;
            else {
                out += utf8_of(code);
            }
        }
    }
    /* Greek lowercase */
    else if (tf == TF_LCGREEK) {
        const char* s = code <= 0xFFFF
                      ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
                      : nullptr;
        if (s) out += s;
        else { char c = (code < 128) ? (char)code : '?'; out += c; }
    }
    /* Greek uppercase */
    else if (tf == TF_UCGREEK) {
        const char* s = code <= 0xFFFF
                      ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
                      : nullptr;
        if (s) out += s;
        else { char c = (code < 128) ? (char)code : '?'; out += c; }
    }
    /* Vector (bold) */
    else if (tf == TF_VECTOR) {
        if ((code >= 'A' && code <= 'Z') || (code >= 'a' && code <= 'z')) {
            out += "\\mathbf{";
            out += (char)code;
            out += '}';
        } else if (code >= 0x20 && code < 0x7F) {
            out += (char)code;
        } else {
            out += "\\mathbf{";
            const char* s = code <= 0xFFFF
                          ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
                          : nullptr;
            if (s) out += s;
            else append_text_char(out, code);
            out += '}';
        }
    }
    /* Function */
    else if (tf == TF_FUNCTION) {
        append_text_char(out, code);
    }
    /* Text */
    else if (tf == TF_TEXT) {
        out += "\\text{";
        append_text_char(out, code);
        out += '}';
    }
    /* Explicit math alphabets.  Unembellished runs are grouped by
     * emitSequence(); this path handles a styled character carrying an
     * accent or other embellishment. */
    else if (const char* alphabet = math_alphabet_command(tf)) {
        out += alphabet;
        out += '{';
        append_math_char(out, code);
        out += '}';
    }
    /* Number */
    else if (tf == TF_NUMBER) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
    }
    /* Variable (italic) */
    else if (tf == TF_VARIABLE) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
        else {
            const char* s = code <= 0xFFFF
                          ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
                          : nullptr;
            if (s) out += s;
            else {
                out += "\\mathit{";
                append_text_char(out, code);
                out += '}';
            }
        }
    }
    /* Display (fence brackets) */
    else if (tf == (TF_DISPLAY & 0x7F)) {
        /* Display chars are handled by fence/pass system, not emitted directly */
    }
    /* Fallback */
    else {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
        else {
            const char* s = code <= 0xFFFF
                          ? map_lookup(UNICODE_MAP, UNICODE_MAP_N, uint16_t(code))
                          : nullptr;
            if (s) out += s;
        }
    }

    /* Apply embellishments (suffix) */
    for (auto et : ch.embells) {
        if (et >= 0 && et < (int)EMBELL_MAP_N && EMBELL_MAP[et].suffix)
            out += EMBELL_MAP[et].suffix;
    }
}

/* ============================================================
 * Template emission — typed nodes
 * ============================================================ */

void LaTeXEmitter::emitFence(const FenceNode& fence, std::string& out) {
    static const char* leftBrackets[] = {
        "\\langle ", "(", "\\{", "[", "|", "\\|",
        "\\lfloor ", "\\lceil ", "[", "]", "]", "[", "("
    };
    static const char* rightBrackets[] = {
        "\\rangle ", ")", "\\}", "]", "|", "\\|",
        "\\rfloor ", "\\rceil ", "[", "]", "[", ")", "]"
    };

    auto clamp_selector = [](int s) { return (s < 0 || s > 12) ? tmPAREN : s; };
    const int sel = clamp_selector(fence.selector);
    const int rsel = clamp_selector(fence.right_selector());

    /* \begin{cases} parses to a left brace wrapped around a cases-layout
     * matrix, and `cases` already draws its own brace.  Write that shape back
     * as the environment so both it and the column split it carries survive
     * save and reopen.  A plain matrix in a brace is left alone: it is not
     * the same thing, and rewriting it would left-align its cells. */
    if (sel == tmBRACE && fence.variation == 1 && fence.content.size() == 1 &&
        fence.content.front() &&
        fence.content.front()->tag() == Node::kMatrix &&
        static_cast<const MatrixNode&>(*fence.content.front()).layoutKind ==
            MatrixNode::kCasesLayout) {
        emitNode(fence.content.front().get(), out);
        return;
    }

    std::string content = emitNodes(fence.content);

    if (fence.variation == 0) {
        out += "\\left";
        out += leftBrackets[sel];
        out += " ";
        out += content;
        out += " \\right";
        out += rightBrackets[rsel];
    } else if (fence.variation == 1) {
        out += "\\left";
        out += leftBrackets[sel];
        out += " ";
        out += content;
        out += " \\right.";
    } else {
        out += "\\left. ";
        out += content;
        out += " \\right";
        out += rightBrackets[rsel];
    }
}

void LaTeXEmitter::emitFrac(const FracNode& frac, std::string& out) {
    std::string n = emitNodes(frac.numer);
    std::string d = emitNodes(frac.denom);
    if (frac.slashed) {
        out += "{}^{"; out += n; out += "}/{}_{"; out += d; out += "}";
    } else {
        out += (frac.display ? "\\dfrac{" : "\\frac{");
        out += n; out += "}{"; out += d; out += "}";
    }
}

void LaTeXEmitter::emitSqrt(const SqrtNode& sq, std::string& out) {
    std::string content = emitNodes(sq.content);
    std::string idx = sq.hasIndex ? emitNodes(sq.index) : std::string();
    if (!idx.empty()) {
        out += "\\sqrt[";
        /* TeX ends an optional argument at the first unbraced ']', so an
         * index containing one has to be wrapped.  Without this, the root
         * lost its index and took the wrong radicand when the file was read
         * back: \sqrt[{]}]{x} came back as \sqrt{]}x. */
        if (idx.find(']') != std::string::npos) {
            out += "{"; out += idx; out += "}";
        } else {
            out += idx;
        }
        out += "]{"; out += content; out += "}";
    } else {
        out += "\\sqrt{"; out += content; out += "}";
    }
}

/* One TeX atom: a single character, or one control word such as \lim.
 * Wrapping an atom in braces is harmless for spacing but not for meaning --
 * {\lim} is an ordinary symbol, so in display style its limits move from
 * under the operator to beside it. */
static bool is_single_tex_atom(const std::string& s) {
    auto letter = [](char c) {
        return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
    };
    if (s.empty()) return false;
    size_t i = 0;
    if (s[0] == '\\') {
        ++i;
        if (i < s.size() && !letter(s[i])) ++i;      /* \{  \%  \,  ... */
        else while (i < s.size() && letter(s[i])) ++i;
    } else {
        const unsigned char lead = static_cast<unsigned char>(s[0]);
        i = lead < 0x80             ? 1
          : (lead & 0xE0) == 0xC0   ? 2
          : (lead & 0xF0) == 0xE0   ? 3
          : (lead & 0xF8) == 0xF0   ? 4 : 1;
        if (i > s.size()) return false;
    }
    while (i < s.size() && s[i] == ' ') ++i;
    return i == s.size();
}

void LaTeXEmitter::emitScript(const ScriptNode& script, std::string& out) {
    std::string base = emitNodes(script.base);
    /* An empty script base still needs an explicit TeX atom.  Emitting only
     * `^{...}` or `_{...}` is invalid at the start of a group and required a
     * second load/save cycle to disappear.  `{}` preserves the structural
     * placeholder and makes serialization reach a fixed point immediately.
     *
     * Whether the rest needs braces is decided from the emitted text, not
     * from the node shape.  One character node can serialize to several
     * tokens -- a text character becomes `\text{...}`, a vector character
     * `\mathbf{x}` -- and those came out unbraced while the same content
     * re-read from a file came out braced, so the file changed every other
     * time it was opened.  Judging the string covers both, and it is also
     * what keeps `\lim` an operator rather than an ordinary atom. */
    /* `\lim` is a single atom but it greedily swallows a following subscript
     * as its own condition, so a script written on it unbraced would rebind
     * to the lim on reparse -- {\lim}_{x} came back as \lim_{x}, a different
     * tree.  Brace such a base so the script stays where it was authored. */
    auto trimmed = [](const std::string& s) {
        size_t a = s.find_first_not_of(' '), b = s.find_last_not_of(' ');
        return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
    };
    const bool scriptGreedyBase = trimmed(base) == "\\lim";

    if (base.empty()) out += "{}";
    else if (scriptGreedyBase || !is_single_tex_atom(base)) {
        out += "{"; out += base; out += "}";
    }
    else out += base;
    if (script.hasSub) {
        std::string sub = emitNodes(script.sub);
        out += "_{"; out += sub; out += "}";
    }
    if (script.hasSup) {
        std::string sup = emitNodes(script.sup);
        out += "^{"; out += sup; out += "}";
    }
}

void LaTeXEmitter::emitIntegral(const IntegralNode& integ, std::string& out) {
    /* Symbol name based on selector */
    const char* sym = "\\int ";
    if (integ.selector == tmDINT) sym = "\\iint ";
    else if (integ.selector == tmTINT) sym = "\\iiint ";
    else if (integ.selector == tmSSINT) sym = "\\oint ";
    else if (integ.selector == tmDSINT) sym = "\\oiint ";
    else if (integ.selector == tmTSINT) sym = "\\oiiint ";

    out += sym;

    if (integ.hasLower || integ.hasUpper) {
        if (integ.hasLower) {
            std::string lo = emitNodes(integ.lower);
            out += "_{"; out += lo; out += "}";
        }
        if (integ.hasUpper) {
            std::string hi = emitNodes(integ.upper);
            out += "^{"; out += hi; out += "}";
        }
    }

    /* Body */
    std::string body = emitNodes(integ.body);
    out += body;
}

void LaTeXEmitter::emitBigOp(const BigOpNode& bigop, std::string& out) {
    const char* sym = "\\sum ";
    int sel = bigop.selector;
    if (sel == tmSUM || sel == tmISUM) sym = "\\sum ";
    else if (sel == tmPROD || sel == tmIPROD) sym = "\\prod ";
    else if (sel == tmCOPROD || sel == tmICOPROD) sym = "\\coprod ";
    else if (sel == tmUNION || sel == tmIUNION) sym = "\\bigcup ";
    else if (sel == tmINTER || sel == tmIINTER) sym = "\\bigcap ";

    out += sym;

    if (bigop.hasLower || bigop.hasUpper) {
        if (bigop.hasLower) {
            std::string lo = emitNodes(bigop.lower);
            out += "_{"; out += lo; out += "}";
        }
        if (bigop.hasUpper) {
            std::string hi = emitNodes(bigop.upper);
            out += "^{"; out += hi; out += "}";
        }
        out += " ";
    }

    std::string body = emitNodes(bigop.body);
    out += body;
}

void LaTeXEmitter::emitDecoration(const DecorationNode& deco, std::string& out) {
    std::string content = emitNodes(deco.content);
    switch (deco.selector) {
    case tmOBAR: out += "\\overline{"; out += content; out += "}"; break;
    case tmUBAR: out += "\\underline{"; out += content; out += "}"; break;
    case tmRARROW: out += "\\overrightarrow{"; out += content; out += "}"; break;
    case tmLARROW: out += "\\overleftarrow{"; out += content; out += "}"; break;
    case tmBARROW: out += "\\overleftrightarrow{"; out += content; out += "}"; break;
    default: out += content; break;
    }
}

void LaTeXEmitter::emitBraceDeco(const BraceDecoNode& bd, std::string& out) {
    std::string content = emitNodes(bd.content);
    std::string label = emitNodes(bd.label);
    const bool over = (bd.selector == tmUHBRACE);
    out += over ? "\\overbrace{" : "\\underbrace{";
    out += content;
    out += "}";
    /* An unlabelled brace must not emit an empty script: `\overbrace{x}^{}`
     * is not what the user wrote, and the parser would rebuild it with a
     * stray empty slot the caret can land in. */
    if (!label.empty()) {
        out += over ? "^{" : "_{";
        out += label;
        out += "}";
    }
}

void LaTeXEmitter::emitDirac(const DiracNode& dirac, std::string& out) {
    std::string bra = emitNodes(dirac.bra);
    std::string ket = emitNodes(dirac.ket);
    if (dirac.variation == 0) {
        out += "\\left\\langle "; out += bra;
        out += " \\middle| "; out += ket;
        out += " \\right\\rangle ";
    } else {
        out += "\\left\\langle "; out += bra; out += " \\right|";
    }
}

void LaTeXEmitter::emitLim(const LimNode& lim, std::string& out) {
    out += "\\lim";
    std::string content = emitNodes(lim.content);
    if (!content.empty()) {
        out += "_{"; out += content; out += "}";
    }
    out += " ";
}

void LaTeXEmitter::emitPile(const PileNode& pile, std::string& out) {
    const char* envName = "gathered";
    if (pile.halign == 1 || pile.kind == 1) envName = "aligned";
    if (pile.halign >= 20) envName = "matrix"; /* bmatrix/pmatrix handled by fence wrapper */

    out += "\\begin{"; out += envName; out += "}\n";
    for (size_t i = 0; i < pile.lines.size(); i++) {
        if (i > 0) out += " \\\\\n";
        emitNode(pile.lines[i].get(), out);
    }
    out += "\n\\end{"; out += envName; out += "}";
}

void LaTeXEmitter::emitMatrix(const MatrixNode& mat, std::string& out) {
    const char* env = mat.layoutKind == MatrixNode::kAlignedLayout ? "aligned"
                    : mat.layoutKind == MatrixNode::kCasesLayout   ? "cases"
                                                                   : "matrix";
    const int rows = std::max(0, mat.rows);
    out += "\\begin{"; out += env; out += "}\n";
    for (int r = 0; r < rows; r++) {
        for (int c = 0; c < mat.cols; c++) {
            if (c > 0) out += " & ";
            int idx = r * mat.cols + c;
            std::string cell;
            if (idx < (int)mat.elements.size())
                emitNode(mat.elements[idx].get(), cell);
            /* A final empty row in a one-column row container cannot be
             * distinguished from an optional final \\ in TeX.  Add one
             * canonical empty group so an Enter-created line survives save
             * and reopen.  Multi-column rows already carry `&` delimiters. */
            if (mat.cols == 1 && r == rows - 1 &&
                cell.find_first_not_of(" \t\r\n") == std::string::npos)
                out += "{}";
            else
                out += cell;
        }
        if (r < rows - 1) out += " \\\\\n";
    }
    out += "\n\\end{"; out += env; out += "}";
}

void LaTeXEmitter::emitEmbell(const EmbellNode& embell, std::string& out) {
    const std::string content = emitNodes(embell.content);
    const int type = static_cast<int>(embell.embellType);
    /* The prime family has no prefix -- it is written entirely as a suffix
     * (x', x'').  Requiring both halves silently dropped the mark, so a
     * prime template produced a bare x. */
    if (type >= 0 && type < static_cast<int>(EMBELL_MAP_N) &&
        EMBELL_MAP[type].suffix) {
        if (EMBELL_MAP[type].prefix) out += EMBELL_MAP[type].prefix;
        out += content;
        out += EMBELL_MAP[type].suffix;
    } else {
        out += content;
    }
}

void LaTeXEmitter::emitOverset(const OversetNode& stacked, std::string& out) {
    out += stacked.under ? "\\underset{" : "\\overset{";
    out += emitNodes(stacked.over);
    out += "}{";
    out += emitNodes(stacked.base);
    out += "}";
}

/* ============================================================
 * Post-processing
 * ============================================================ */
std::string LaTeXEmitter::postProcess(const std::string& raw) {
    std::string s = raw;
    /* Whitespace outside a math expression is semantically empty.  Some
     * legacy symbol spellings include padding on both sides; trimming only
     * the tail made their first serialization differ from the next parse. */
    auto outer_space = [](char c) {
        return c == ' ' || c == '\t' || c == '\r' || c == '\n';
    };
    size_t first = 0;
    while (first < s.size() && outer_space(s[first])) ++first;
    if (first) s.erase(0, first);
    while (!s.empty() && outer_space(s.back()))
        s.pop_back();
    return s;
}

/* ============================================================
 * Symbol lookup helpers (delegating to map_lookup)
 * ============================================================ */
const char* LaTeXEmitter::lookupSymbol(uint16_t code) const {
    const char* s = map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, code);
    if (s) return s;
    return map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
}

const char* LaTeXEmitter::lookupGreekLower(uint16_t code) const {
    return map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
}

const char* LaTeXEmitter::lookupGreekUpper(uint16_t code) const {
    return map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
}

std::string tree_to_latex(const LineNode& root) {
    LaTeXEmitter em;
    return em.emit(root);
}

} /* namespace eqnedit */
