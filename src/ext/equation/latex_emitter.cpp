/*
 * latex_emitter.cpp -- Node tree → LaTeX text conversion
 *
 * Emits normalized TeX from the structural editing tree.
 */
#include <algorithm>
#include "latex_emitter.h"
#include "math_layout.h"      /* is_cjk: the same rule that picks the face */
#include <cstring>
#include <cstdio>

namespace mtef {

/* ============================================================
 * Symbol lookup tables used by the normalized TeX emitter
 * ============================================================ */
struct MapEntry { uint16_t key; const char* val; };

static const MapEntry SYMBOL_MAP[] = {
    {0x0021, "!"}, {0x0022, "\""}, {0x0023, "#"}, {0x0025, "\\%"},
    {0x0026, "\\&"}, {0x0028, "("}, {0x0029, ")"}, {0x002A, "*"},
    {0x002B, "+"}, {0x002C, ","}, {0x002D, "-"}, {0x002E, "."},
    {0x002F, "/"}, {0x003A, ":"}, {0x003B, ";"}, {0x003C, "<"},
    {0x003D, " = "}, {0x003E, ">"}, {0x003F, "?"}, {0x0040, "@"},
    {0x005B, "["}, {0x005D, "]"}, {0x005E, " \\wedge "},
    {0x007B, "\\{"}, {0x007C, "|"}, {0x007D, "\\}"},
    {0x007E, " \\sim "},
    {0x00AC, " \\neg "}, {0x00B0, "^{\\circ }"}, {0x00B1, " \\pm "},
    {0x00B7, " \\cdot "}, {0x00D7, " \\times "}, {0x00F7, " \\div "},
};
#define SYMBOL_MAP_N (sizeof(SYMBOL_MAP)/sizeof(SYMBOL_MAP[0]))

static const MapEntry UNICODE_MAP[] = {
    {0x0393, "\\Gamma "}, {0x0394, "\\Delta "}, {0x0398, "\\Theta "},
    {0x039B, "\\Lambda "}, {0x039E, "\\Xi "}, {0x03A0, "\\Pi "},
    {0x03A3, "\\Sigma "}, {0x03A5, "\\Upsilon "}, {0x03A6, "\\Phi "},
    {0x03A8, "\\Psi "}, {0x03A9, "\\Omega "},
    {0x03B1, "\\alpha "}, {0x03B2, "\\beta "}, {0x03B3, "\\gamma "},
    {0x03B4, "\\delta "}, {0x03B5, "\\varepsilon "}, {0x03B6, "\\zeta "},
    {0x03B7, "\\eta "}, {0x03B8, "\\theta "}, {0x03B9, "\\iota "},
    {0x03BA, "\\kappa "}, {0x03BB, "\\lambda "}, {0x03BC, "\\mu "},
    {0x03BD, "\\nu "}, {0x03BE, "\\xi "}, {0x03C0, "\\pi "},
    {0x03C1, "\\rho "}, {0x03C2, "\\varsigma "}, {0x03C3, "\\sigma "},
    {0x03C4, "\\tau "}, {0x03C5, "\\upsilon "}, {0x03C6, "\\phi "},
    {0x03C7, "\\chi "}, {0x03C8, "\\psi "}, {0x03C9, "\\omega "},
    /* The explicit spaces, so they read back as what was typed.  The
     * table is binary-searched, so these sit in numeric order like
     * everything else -- put at the top they were simply never found. */
    {0x2003, "\\quad "}, {0x2005, "\\; "},
    {0x2006, "\\, "},
    {0x2019, "'"}, {0x2032, "'"}, {0x2033, "''"}, {0x2034, "'''"},
    {0x205F, "\\: "},
    {0x2102, "\\mathbb{C}"}, {0x210D, "\\mathbb{H}"},
    {0x2111, "\\Im "}, {0x2113, "\\ell "},
    {0x2115, "\\mathbb{N}"}, {0x2118, "\\wp "}, {0x2119, "\\mathbb{P}"},
    {0x211A, "\\mathbb{Q}"}, {0x211C, "\\Re "}, {0x211D, "\\mathbb{R}"},
    {0x2124, "\\mathbb{Z}"},
    {0x2135, "\\aleph "}, {0x2190, "\\leftarrow "},
    {0x2191, "\\uparrow "}, {0x2192, "\\rightarrow "},
    {0x2193, "\\downarrow "}, {0x2194, "\\leftrightarrow "},
    {0x2195, "\\updownarrow "}, {0x2196, "\\nwarrow "},
    {0x2197, "\\nearrow "}, {0x2198, "\\searrow "},
    {0x2199, "\\swarrow "},
    /* Harpoons.  Equation Editor draws its vector arrow as one, so these turn
     * up in every converted document; without a name they came out as the
     * bare character.  In order: the table is searched by halving. */
    {0x21BC, "\\leftharpoonup "}, {0x21BD, "\\leftharpoondown "},
    {0x21C0, "\\rightharpoonup "}, {0x21C1, "\\rightharpoondown "},
    /* 0x21CC is deliberately absent: it is the base character of
     * \xrightleftharpoons, and naming it here took that name away. */
    {0x21D0, "\\Leftarrow "},
    {0x21D2, "\\Rightarrow "}, {0x21D4, "\\Leftrightarrow "},
    {0x2200, " \\forall "}, {0x2202, "\\partial "},
    {0x2203, " \\exists "}, {0x2205, "\\emptyset "},
    {0x2207, "\\nabla "}, {0x2208, " \\in "},
    {0x2209, " \\notin "}, {0x220B, " \\ni "},
    {0x2210, " \\coprod "}, {0x2211, "\\sum "},
    {0x2212, " - "}, {0x2213, " \\mp "},
    {0x2215, "/"}, {0x2216, " \\setminus "},
    {0x2217, " * "}, {0x2218, " \\circ "},
    {0x221A, "\\sqrt{}"}, {0x221D, " \\propto "},
    {0x221E, "\\infty "}, {0x2220, "\\angle "},
    {0x2223, " \\mid "}, {0x2225, " \\| "},
    {0x2227, " \\wedge "}, {0x2228, " \\vee "},
    {0x2229, " \\cap "}, {0x222A, " \\cup "},
    {0x222B, "\\int "}, {0x222C, "\\iint "},
    {0x222D, "\\iiint "}, {0x222E, "\\oint "},
    {0x222F, "\\oiint "}, {0x2230, "\\oiiint "},
    {0x2234, " \\therefore "}, {0x2235, " \\because "},
    {0x223C, " \\sim "}, {0x2243, " \\simeq "},
    {0x2245, " \\cong "}, {0x2248, " \\approx "},
    {0x2260, " \\neq "}, {0x2261, " \\equiv "},
    {0x2264, " \\leq "}, {0x2265, " \\geq "},
    {0x226A, " \\ll "}, {0x226B, " \\gg "},
    {0x2282, " \\subset "}, {0x2283, " \\supset "},
    {0x2286, " \\subseteq "}, {0x2287, " \\supseteq "},
    {0x2295, " \\oplus "}, {0x2296, " \\ominus "},
    {0x2297, " \\otimes "}, {0x2299, " \\odot "},
    {0x22A5, " \\perp "}, {0x22C5, " \\cdot "},
    {0x22EE, "\\vdots "}, {0x22EF, "\\cdots "},
    {0x22F0, "\\iddots "}, {0x22F1, "\\ddots "},
    {0x2329, "\\langle "}, {0x232A, "\\rangle"},
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

/* Write a code point as UTF-8.
 *
 * The fallback for a character with no LaTeX command has to be the character
 * itself.  '?' loses it outright, and \symbol{N} neither round trips nor reads
 * -- both are the same defect that dropped every \vec and every kanji: an
 * output path that quietly discards what the author typed.
 *
 * Uppercase Greek is where this bites in practice.  Beta, Eta, Rho and the
 * rest look like Latin letters, so LaTeX has no \Beta to write them with, and
 * they were coming out as '?'. */
void emit_utf8(uint32_t cp, std::string& out) {
    if (cp < 0x80) { out += char(cp); return; }
    if (cp < 0x800) {
        out += char(0xC0 | (cp >> 6));
    } else {
        out += char(0xE0 | (cp >> 12));
        out += char(0x80 | ((cp >> 6) & 0x3F));
    }
    out += char(0x80 | (cp & 0x3F));
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
    /* EM_R1ARROW 14 */ {"\\overset{\\rightharpoonup}{", "}"},
    /* EM_L1ARROW 15 */ {"\\overset{\\leftharpoonup}{", "}"},
    /* EM_MBAR  16 */ {"\\cancel{", "}"},
    /* EM_OBAR  17 */ {"\\overline{", "}"},
    /* EM_TPRIME 18 */ {nullptr, "'''"},
    /* EM_FROWN 19 */ {"\\overset{\\frown}{", "}"},
    /* EM_SMILE 20 */ {"\\overset{\\smile}{", "}"},
};
#define EMBELL_MAP_N (sizeof(EMBELL_MAP)/sizeof(EMBELL_MAP[0]))

/* ============================================================
 * Constructor
 * ============================================================ */
LaTeXEmitter::LaTeXEmitter(int prodVer, bool run_passes)
    : prodVer_(prodVer), runPasses_(run_passes) {}

/* ============================================================
 * Emit NodeList → string
 * ============================================================ */
std::string LaTeXEmitter::emitNodes(const NodeList& nodes) {
    std::string s;
    for (auto& n : nodes) emitNode(n.get(), s);
    return s;
}

/* ============================================================
 * Main entry point
 * ============================================================ */
std::string LaTeXEmitter::emit(const LineNode& root) {
    std::string raw;
    depth_ = 0;
    emitLine(root, raw);
    return postProcess(raw);
}

/* ============================================================
 * LINE emission (simplified — no 7-pass system yet)
 * ============================================================ */
void LaTeXEmitter::emitLine(const LineNode& line, std::string& out) {
    if (line.isNull) return;

    /* Run pass pipeline on a mutable copy of children.
     * The pipeline resolves EQNEDT32 patterns:
     *   Pass 1: fence/decoration empty-slot → merge sibling content
     *   Pass 2: BigOp display data (SIZE_SUB path) → set displayLower/Upper
     *   Pass 2b: BigOp display data (no SIZE_SUB) → set displayLower/Upper */
    NodeList& mutableChildren = const_cast<NodeList&>(line.children);
    if (runPasses_) pipeline_.process(mutableChildren, depth_, prodVer_);

    depth_++;

    /* Output phase with text/function grouping.
     * Consecutive TF_TEXT chars → \text{...}
     * Consecutive TF_FUNCTION chars → \sin, \cos, \operatorname{...} */
    std::string textBuf;
    std::string funcBuf;
    int funcEmbell = -1;  /* embellishment on function group */

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
        /* Longest-prefix matching */
        std::string remaining = funcBuf;
        while (!remaining.empty()) {
            bool matched = false;
            /* Try known funcs */
            for (int i = 0; KNOWN_FUNCS[i]; i++) {
                size_t len = strlen(KNOWN_FUNCS[i]);
                if (remaining.substr(0, len) == KNOWN_FUNCS[i]) {
                    out += "\\"; out += KNOWN_FUNCS[i]; out += " ";
                    remaining = remaining.substr(len);
                    matched = true; break;
                }
            }
            if (matched) continue;
            /* Try known operators */
            for (int i = 0; KNOWN_OPS[i]; i++) {
                size_t len = strlen(KNOWN_OPS[i]);
                if (remaining.substr(0, len) == KNOWN_OPS[i]) {
                    out += "\\operatorname{"; out += KNOWN_OPS[i]; out += "}";
                    remaining = remaining.substr(len);
                    matched = true; break;
                }
            }
            if (matched) continue;
            /* No match: output single char */
            out += remaining[0];
            remaining = remaining.substr(1);
        }
        if (funcEmbell >= 0 && funcEmbell < (int)EMBELL_MAP_N) {
            /* Wrap the entire function output with embellishment */
        }
        funcBuf.clear();
        funcEmbell = -1;
    };

    for (auto& child : line.children) {
        if (!child) continue;
        if (child->tag() == Node::kChar) {
            auto* ch = static_cast<const CharNode*>(child.get());
            if (ch->typeface == TF_TEXT) {
                flushFunc();
                if (ch->charCode >= 0x20 && ch->charCode < 0x7F)
                    textBuf += (char)ch->charCode;
                continue;
            }
            if (ch->typeface == TF_FUNCTION) {
                flushText();
                if (ch->charCode >= 0x20 && ch->charCode < 0x7F)
                    funcBuf += (char)ch->charCode;
                continue;
            }
        }
        /* Non-text/func node: flush buffers and emit */
        flushText();
        flushFunc();
        emitNode(child.get(), out);
    }
    flushText();
    flushFunc();

    depth_--;
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
    case Node::kOverset: {
        const OversetNode& o = *static_cast<const OversetNode*>(node);
        const std::string over = emitNodes(o.over);
        /* An arrow labelled above AND below is one command with an
         * optional argument, not an \underset wrapped round an
         * \overset -- which is how it is BUILT, but not how it is
         * written, and the written form is what the author reads back. */
        if (o.under && o.base.size() == 1 && o.base[0] &&
            o.base[0]->tag() == Node::kOverset) {
            const OversetNode& in = static_cast<const OversetNode&>(*o.base[0]);
            if (!in.under && in.base.size() == 1 && in.base[0] &&
                in.base[0]->tag() == Node::kChar) {
                const CharNode& ch = static_cast<const CharNode&>(*in.base[0]);
                const unsigned cp = ch.charCode ? ch.charCode : unsigned(ch.ch);
                if (const char* name = arrow_command(cp)) {
                    out += name;
                    out += "[" + over + "]{" + emitNodes(in.over) + "}";
                    break;
                }
            }
        }
        /* An arrow with a label over it goes back out as \xrightarrow,
         * because that is what an author types; anything else is the
         * \overset it was built from. */
        int arrow = 0;
        if (!o.under && o.base.size() == 1 && o.base[0] &&
            o.base[0]->tag() == Node::kChar) {
            const CharNode& c = static_cast<const CharNode&>(*o.base[0]);
            const unsigned cp = c.charCode ? c.charCode : unsigned(c.ch);
            arrow = arrow_command(cp) ? int(cp) : 0;
        }
        if (arrow) {
            out += arrow_command(unsigned(arrow));
            out += "{" + over + "}";
            break;
        }
        out += o.under ? "\\underset{" : "\\overset{";
        out += over;
        out += "}{";
        out += emitNodes(o.base);
        out += "}";
        break;
    }
    case Node::kPhantom: {
        const PhantomNode& ph = *static_cast<const PhantomNode*>(node);
        out += ph.keepWidth ? (ph.keepHeight ? "\\phantom{" : "\\hphantom{")
                            : "\\vphantom{";
        out += emitNodes(ph.content);
        out += "}";
        break;
    }
    case Node::kBraceDeco:
        emitBraceDeco(*static_cast<const BraceDecoNode*>(node), out);
        break;
    case Node::kDirac:
        emitDirac(*static_cast<const DiracNode*>(node), out);
        break;
    case Node::kLim:
        emitLim(*static_cast<const LimNode*>(node), out);
        break;
    case Node::kPile:
        emitPile(*static_cast<const PileNode*>(node), out);
        break;
    case Node::kMatrix:
        emitMatrix(*static_cast<const MatrixNode*>(node), out);
        break;
    case Node::kSize: {
        /* A size marker is a SWITCH -- everything after it in the group is set
         * at that size -- and that is exactly what TeX's \scriptstyle is, so
         * the two map onto each other and the setting survives a save.
         *
         * Only three of Equation Editor's five sizes are written.  Its Symbol
         * and Sub-Symbol sizes exist to say how big a summation sign should
         * be, and that is now read from the font (displayOperatorMinHeight),
         * so a marker for them would say nothing the layout would act on. */
        const int t = static_cast<const SizeNode*>(node)->sizeType;
        /* Equation Editor writes one of these between the parts of every
         * display construct, so most of them repeat the size already in
         * force.  Writing them anyway put a dozen no-op \displaystyle into a
         * converted lecture file, in front of a person who then has to read
         * it. */
        if (t == sizeStyle_) break;
        if (t == SIZETYPE_SUB)       { out += "\\scriptstyle ";       sizeStyle_ = t; }
        else if (t == SIZETYPE_SUB2) { out += "\\scriptscriptstyle "; sizeStyle_ = t; }
        else if (t == SIZETYPE_FULL) { out += "\\displaystyle ";      sizeStyle_ = t; }
        break;
    }
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
    uint16_t code = ch.charCode;
    int tf = ch.typeface;

    /* The double-struck and script alphabets are recorded as a typeface on an
     * ordinary letter, so the command is read back off the typeface.  One per
     * letter rather than one per run: \mathbb{R}\mathbb{C} says the same
     * thing as \mathbb{RC} and needs no state carried between characters. */
    if (tf == TF_USER1 || tf == TF_USER2) {
        out += (tf == TF_USER1) ? "\\mathbb{" : "\\mathcal{";
        emit_utf8(uint32_t(code ? code : (unsigned char)ch.ch), out);
        out += "}";
        return;
    }

    /* Apply embellishments (prefix) */
    for (auto it = ch.embells.rbegin(); it != ch.embells.rend(); ++it) {
        int et = *it;
        if (et >= 0 && et < (int)EMBELL_MAP_N && EMBELL_MAP[et].prefix)
            out += EMBELL_MAP[et].prefix;
    }

    /* Symbol typeface */
    if (tf == TF_SYMBOL) {
        const char* s = map_lookup(SYMBOL_MAP, SYMBOL_MAP_N, code);
        if (s) { out += s; }
        else {
            s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
            if (s) out += s;
            else if (is_cjk(uint32_t(code))) {
                /* Japanese is written as itself.  \symbol{30913} is the right
                 * answer for an unmapped MATHS symbol and the wrong one for a
                 * character the author typed: it does not round trip, and no
                 * one can read it.  The same is_cjk that picks the face
                 * decides this, so one rule governs both. */
                emit_utf8(uint32_t(code), out);
            }
            else if (code >= 0x80) {
                /* Any character with no command is written as itself.
                 * \symbol{977} is not an escape hatch -- the parser does not
                 * read it back, so it loses the character just as surely as
                 * '?' did.  This was the last hole of the three: \vec, kanji,
                 * and now the Greek variants like theta-symbol that fall
                 * outside the plain alphabet ranges. */
                emit_utf8(uint32_t(code), out);
            }
            else {
                char buf[32];
                snprintf(buf, sizeof(buf), "\\symbol{%d}", code);
                out += buf;
            }
        }
    }
    /* Greek lowercase */
    else if (tf == TF_LCGREEK) {
        const char* s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
        if (s) out += s;
        else if (code < 128) out += char(code);
        else emit_utf8(uint32_t(code), out);   /* see emit_utf8: never '?' */
    }
    /* Greek uppercase */
    else if (tf == TF_UCGREEK) {
        const char* s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
        if (s) out += s;
        else if (code < 128) out += char(code);
        else emit_utf8(uint32_t(code), out);   /* see emit_utf8: never '?' */
    }
    /* Vector -- bold ITALIC, which is \bm.  Writing \mathbf here set every
     * vector upright, so the paste and the picture disagreed about the one
     * face the lab has a rule for. */
    else if (tf == TF_VECTOR || tf == TF_USER3) {
        const char* cmd = (tf == TF_VECTOR) ? "\\bm{" : "\\mathbf{";
        if ((code >= 'A' && code <= 'Z') || (code >= 'a' && code <= 'z') ||
            (code >= '0' && code <= '9')) {
            out += cmd;
            out += (char)code;
            out += '}';
        } else if (code >= 0x20 && code < 0x7F) {
            out += (char)code;
        } else {
            out += cmd;
            emit_utf8(uint32_t(code), out);
            out += '}';
        }
    }
    /* Function */
    else if (tf == TF_FUNCTION) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
    }
    /* Text */
    else if (tf == TF_TEXT) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
    }
    /* Number */
    else if (tf == TF_NUMBER) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
    }
    /* Variable (italic) */
    else if (tf == TF_VARIABLE) {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
        else {
            const char* s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
            if (s) out += s;
        }
    }
    /* Display (fence brackets) */
    else if (tf == TF_DISPLAY) {
        /* Display chars are handled by fence/pass system, not emitted directly */
    }
    /* Fallback */
    else {
        if (code >= 0x20 && code < 0x7F) out += (char)code;
        else {
            const char* s = map_lookup(UNICODE_MAP, UNICODE_MAP_N, code);
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

    int sel = fence.selector;
    if (sel < 0 || sel > 12) sel = 1; /* default to parens */

    /* Parentheses round a ruleless fraction and nothing else IS \\binom, and
     * writing it that way is what makes it read back.  The parts on their own
     * would go out as {n \\atop k}, and \\atop is infix -- the parser has no way
     * to take it, so the equation came back with the word "atop" in it. */
    if (sel == tmPAREN && fence.variation == 0 && fence.content.size() == 1) {
        const Node* only = fence.content.front().get();
        if (only && only->tag() == Node::kFrac &&
            !static_cast<const FracNode*>(only)->ruled) {
            const FracNode& f = *static_cast<const FracNode*>(only);
            out += "\\binom{"; out += emitNodes(f.numer);
            out += "}{";       out += emitNodes(f.denom); out += "}";
            return;
        }
    }

    std::string content = emitNodes(fence.content);
    if (fence.hasMiddle) {
        content += " \\middle";
        content += (fence.middle == 0x2016) ? "\\| " : "| ";
        content += emitNodes(fence.content2);
    }

    if (fence.variation == 0) {
        out += "\\left";
        out += leftBrackets[sel];
        out += " ";
        out += content;
        out += " \\right";
        out += rightBrackets[sel];
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
        out += rightBrackets[sel];
    }
}

/* The amsmath name for an extensible arrow, or nullptr. */
const char* LaTeXEmitter::arrow_command(unsigned cp) {
    switch (cp) {
        case 0x2192: return "\\xrightarrow";
        case 0x2190: return "\\xleftarrow";
        case 0x2194: return "\\xleftrightarrow";
        case 0x21D2: return "\\xRightarrow";
        case 0x21D0: return "\\xLeftarrow";
        case 0x21D4: return "\\xLeftrightarrow";
        case 0x21A6: return "\\xmapsto";
        case 0x21CC: return "\\xrightleftharpoons";
        default: return nullptr;
    }
}

void LaTeXEmitter::emitFrac(const FracNode& frac, std::string& out) {
    ++fracDepth_;
    std::string n = emitNodes(frac.numer);
    std::string d = emitNodes(frac.denom);
    --fracDepth_;
    if (!frac.ruled) {
        /* A ruleless fraction is a binomial; the parentheses round it are
         * the fence node outside, so \binom would double them.  {a \atop b}
         * is the ruleless fraction on its own, which is what this is. */
        out += "{"; out += n; out += " \\atop "; out += d; out += "}";
    } else if (frac.slashed) {
        out += "{}^{"; out += n; out += "}/{}_{"; out += d; out += "}";
    } else {
        /* \dfrac is the default, because the outermost fraction is DRAWN
         * at display size and a bare \frac pasted into running text is
         * not.  Inside another fraction the name is dropped again: LaTeX
         * already steps a nested \frac down, so the two rules agree level
         * by level and the paste matches the picture at each one. */
        const bool display = frac.styleOverride ? (frac.styleOverride > 0)
                                                : (fracDepth_ == 0);
        out += display    ? "\\dfrac{"
             : fracDepth_ ? "\\frac{"
                          : "\\tfrac{";
        out += n; out += "}{"; out += d; out += "}";
    }
}

/* A standalone embellishment: \vec{B}, \hat{n}, \bar{A}, \dot{x}.
 *
 * This case used to emit nothing, on the reasoning that an embellishment is
 * "usually attached to CHAR".  That is true of a tree read from MTEF, where
 * Equation Editor hangs the accent off the character -- and false of a tree
 * from the LaTeX parser, which builds a node of its own.  So every vector an
 * author typed was silently dropped on the way back out to LaTeX: the picture
 * and the Office paste were right and the saved file had lost the arrow.
 *
 * An assumption that held for one producer and not the other, which is why it
 * survived: the MTEF corpus tests never exercised this path. */
void LaTeXEmitter::emitEmbell(const EmbellNode& em, std::string& out) {
    const EmbellFmt* fmt =
        (em.embellType >= 0 && em.embellType < (int)EMBELL_MAP_N)
            ? &EMBELL_MAP[em.embellType] : nullptr;

    /* A prime is a suffix with no prefix (x'), an accent is both (\hat{x}),
     * and an unknown type still has to keep its content rather than lose it. */
    if (fmt && fmt->prefix) out += fmt->prefix;
    for (const auto& child : em.content) emitNode(child.get(), out);
    if (fmt && fmt->suffix) out += fmt->suffix;
}

void LaTeXEmitter::emitSqrt(const SqrtNode& sq, std::string& out) {
    std::string content = emitNodes(sq.content);
    if (sq.hasIndex && !sq.index.empty()) {
        std::string idx = emitNodes(sq.index);
        out += "\\sqrt["; out += idx; out += "]{"; out += content; out += "}";
    } else {
        out += "\\sqrt{"; out += content; out += "}";
    }
}

void LaTeXEmitter::emitScript(const ScriptNode& script, std::string& out) {
    std::string base = emitNodes(script.base);
    /* {}^{14}_{6}C -- the empty group is what carries the scripts, and it is
     * what makes them attach to the LEFT of the C rather than to whatever
     * happened to come before. */
    if (script.pre) {
        out += "{}";
        if (script.hasSup) { out += "^{"; out += emitNodes(script.sup); out += "}"; }
        if (script.hasSub) { out += "_{"; out += emitNodes(script.sub); out += "}"; }
        out += base;
        return;
    }
    out += base;
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

    /* Display data limits (populated by Pass 2 — from old pipeline) */
    if (integ.displayLower || integ.displayUpper) {
        if (integ.displayLower) {
            out += "\\limits_{";
            emitNode(integ.displayLower.get(), out);
            out += "}";
        }
        if (integ.displayUpper) {
            out += "^{";
            emitNode(integ.displayUpper.get(), out);
            out += "}";
        }
    }
    /* Slot-based limits (from MTEF binary variation bits) */
    else if (integ.hasLower || integ.hasUpper) {
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

    /* Display data limits */
    if (bigop.displayLower || bigop.displayUpper) {
        out += "\\limits";
        if (bigop.displayLower) {
            out += "_{";
            emitNode(bigop.displayLower.get(), out);
            out += "}";
        }
        if (bigop.displayUpper) {
            out += "^{";
            emitNode(bigop.displayUpper.get(), out);
            out += "}";
        }
        out += " ";
    }
    /* Slot-based limits */
    else if (bigop.hasLower || bigop.hasUpper) {
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
    if (bd.selector == tmUHBRACE) {
        out += "\\overbrace{"; out += content; out += "}^{"; out += label; out += "}";
    } else {
        out += "\\underbrace{"; out += content; out += "}_{"; out += label; out += "}";
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
    /* Flush left and flush right need no environment of their own: an
     * `aligned` row whose content is all in the SECOND column is flush left,
     * and all in the first is flush right.  Writing them that way means the
     * alignment survives being saved and read back, where a `gathered` would
     * come back centred. */
    const bool left  = (pile.halign == 2);
    const bool right = (pile.halign == 3);
    if (left || right) envName = "aligned";
    if (pile.halign >= 20) envName = "matrix"; /* bmatrix/pmatrix handled by fence wrapper */

    out += "\\begin{"; out += envName; out += "}\n";
    /* A real alignment holds its cells in row-major order, so the & goes
     * between them and the row break every ncols. */
    const int cols = (pile.halign == 1) ? std::max(1, pile.ncols) : 1;
    for (size_t i = 0; i < pile.lines.size(); i++) {
        if (i > 0) {
            if (int(i) % cols == 0) out += " \\\\\n";
            else                    out += " & ";
        }
        if (left) out += "&";
        emitNode(pile.lines[i].get(), out);
        if (right) out += "&";
    }
    out += "\n\\end{"; out += envName; out += "}";
}

void LaTeXEmitter::emitMatrix(const MatrixNode& mat, std::string& out) {
    out += "\\begin{matrix}\n";
    for (int r = 0; r < mat.rows; r++) {
        for (int c = 0; c < mat.cols; c++) {
            if (c > 0) out += " & ";
            int idx = r * mat.cols + c;
            if (idx < (int)mat.elements.size())
                emitNode(mat.elements[idx].get(), out);
        }
        if (r < mat.rows - 1) out += " \\\\\n";
    }
    out += "\n\\end{matrix}";
}

/* ============================================================
 * Post-processing
 * ============================================================ */
std::string LaTeXEmitter::postProcess(const std::string& raw) {
    std::string s = raw;
    /* Trim trailing whitespace */
    while (!s.empty() && (s.back() == ' ' || s.back() == '\n'))
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
    LaTeXEmitter em(3, /*run_passes=*/false);
    return em.emit(root);
}

} /* namespace mtef */
