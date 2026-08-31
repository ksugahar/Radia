/*
 * math_writer.cpp -- the shared tree walk behind every Office math output
 *
 * See math_writer.h for why this is one walk rather than one per format.
 */
#include "math_writer.h"
#include "line_pass.h"
#include "mtef_common.h"

#include <string>
#include <vector>

namespace mtef {
namespace {

/* Delimiters.  A null pair means the default parentheses. */
struct Fence { const char* beg; const char* end; };
Fence fence_chars(int selector) {
    switch (selector) {
        case tmPAREN: return {"(", ")"};
        case tmBRACK: return {"[", "]"};
        case tmBRACE: return {"{", "}"};
        case tmANGLE: return {"\xE2\x9F\xA8", "\xE2\x9F\xA9"};   /* U+27E8/9 */
        case tmBAR:   return {"|", "|"};
        case tmDBAR:  return {"\xE2\x80\x96", "\xE2\x80\x96"};   /* U+2016 */
        case tmFLOOR: return {"\xE2\x8C\x8A", "\xE2\x8C\x8B"};   /* U+230A/B */
        case tmCEIL:  return {"\xE2\x8C\x88", "\xE2\x8C\x89"};   /* U+2308/9 */
        default:      return {"(", ")"};
    }
}

uint32_t nary_char(int selector) {
    switch (selector) {
        case tmSINT:  return 0x222B;
        case tmDINT:  return 0x222C;
        case tmTINT:  return 0x222D;
        case tmSSINT: return 0x222E;
        case tmDSINT: return 0x222F;
        case tmTSINT: return 0x2230;
        case tmSUM: case tmISUM:       return 0x2211;
        case tmPROD: case tmIPROD:     return 0x220F;
        case tmCOPROD: case tmICOPROD: return 0x2210;
        case tmUNION: case tmIUNION:   return 0x22C3;
        case tmINTER: case tmIINTER:   return 0x22C2;
        default: return 0x2211;
    }
}

}  // namespace

/* MTEF embellishment code -> the combining accent, for markup that positions
 * the accent itself. */
uint32_t accent_char(int embellType) {
    switch (embellType) {
        case 2:  return 0x0307;   /* dot         */
        case 3:  return 0x0308;   /* double dot  */
        case 4:  return 0x20DB;   /* triple dot  */
        case 8:  return 0x0303;   /* tilde       */
        case 9:  return 0x0302;   /* hat         */
        case 11: return 0x20D7;   /* right arrow */
        case 12: return 0x20D6;   /* left arrow  */
        case 13: return 0x20E1;   /* both arrows */
        case 17: return 0x0304;   /* bar         */
        default: return 0;
    }
}

/* The same accents as SPACING characters, for a renderer that has to place
 * them itself.  Handing a combining character to TextOut gives it nothing to
 * combine with, so it lands as a stray mark or as nothing at all -- which is
 * how every vector came to draw as a bare letter.
 *
 * A bar returns 0: it is a rule the width of what it covers, not a glyph, and
 * a fixed-width macron over a wide expression looks wrong. */
uint32_t accent_drawing_char(int embellType) {
    /* The COMBINING forms, which is what TeX sets and what Word writes in an
     * OMML m:chr.  The spacing forms that were here are a different drawing:
     * a spacing circumflex is 8.30 pt where the combining one is 8.81, and a
     * spacing dot 7.88 where the combining one is 8.12 -- so every accented
     * letter came out a little short, and \vec, which had been an ordinary
     * rightwards arrow, came out two and a quarter points short.
     *
     * A combining mark carries no advance of its own; where it sits is the
     * font's MathTopAccentAttachment, which is read for exactly this. */
    switch (embellType) {
        case 2:  return 0x0307;   /* dot           */
        case 3:  return 0x0308;   /* double dot    */
        case 4:  return 0x20DB;   /* triple dot    */
        case 8:  return 0x0303;   /* tilde         */
        case 9:  return 0x0302;   /* hat           */
        case 11: return 0x20D7;   /* right arrow   */
        case 12: return 0x20D6;   /* left arrow    */
        case 13: return 0x20E1;   /* both arrows   */
        /* A bar over ONE letter is the combining macron, which is what TeX
         * sets for ar -- a fixed drawing, centred, narrower than the
         * letter.  A rule the width of what it covers is \overline, and
         * that is a different construct with its own node.  Drawing this
         * one as a rule too made it 0.024 pt tall and, more to the point,
         * as wide as whatever it sat on. */
        case 17: return 0x0304;   /* bar: combining macron */
        default: return 0;
    }
}

namespace {

/* MTEF marks style by typeface; every output marks it per run.
 * 0 upright, 1 italic, 2 bold italic. */
int run_style(int tf) {
    if (tf == TF_VECTOR) return 2;
    if (tf == TF_VARIABLE || tf == TF_LCGREEK || tf == TF_UCGREEK) return 1;
    return 0;
}

std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }

class Walker {
public:
    Walker(const MathSyntax& syn, bool run_passes)
        : syn_(syn), runPasses_(run_passes) {}

    std::string run(const LineNode& root, bool display) {
        return syn_.document(line(root), display);
    }

private:
    const MathSyntax& syn_;
    bool runPasses_;
    int depth_ = 0;
    PassPipeline pipeline_;

    /* ---- sequence assembly -------------------------------------------- */

    /* One LINE: run the passes, then assemble its children.  The passes move
     * sibling content into empty slots and drop the siblings, which is why the
     * child list is taken as mutable. */
    std::string line(const LineNode& ln) {
        if (ln.isNull) return std::string();
        NodeList& children = const_cast<NodeList&>(ln.children);
        if (runPasses_) pipeline_.process(children, depth_, 3);
        depth_++;
        std::string s = seq(children);
        depth_--;
        return s;
    }

    std::string seq(const NodeList& list, size_t start = 0) {
        std::vector<std::string> parts;
        for (size_t i = start; i < list.size(); ++i) {
            const Node* n = list[i].get();
            if (!n) continue;

            /* MTEF writes a script's base as the preceding sibling.  LaTeX does
             * not care (`x` then `_{c}` binds), but an empty base slot renders
             * as a placeholder box, so it is absorbed here. */
            if (n->tag() == Node::kScript) {
                const auto& s = static_cast<const ScriptNode&>(*n);
                if (s.base.empty() && !parts.empty()) {
                    std::string base = parts.back();
                    parts.pop_back();
                    parts.push_back(script(s, base));
                    continue;
                }
            }

            if (n->tag() == Node::kIntegral || n->tag() == Node::kBigOp) {
                /* Limits the template itself does not carry arrive as a
                 * base-less script right after the operator. */
                const ScriptNode* limits = nullptr;
                if (nary_limits_empty(*n) && i + 1 < list.size() &&
                    list[i + 1] && list[i + 1]->tag() == Node::kScript) {
                    const auto& s = static_cast<const ScriptNode&>(*list[i + 1]);
                    if (s.base.empty()) { limits = &s; ++i; }
                }
                /* The operand is normally in the operator's own body slot; when
                 * the passes could not fill it, it is the sibling run that
                 * follows and has to be absorbed. */
                std::string body = nary_body(*n);
                if (body.empty()) {
                    parts.push_back(nary(*n, seq(list, i + 1), limits));
                    break;
                }
                parts.push_back(nary(*n, body, limits));
                continue;
            }

            std::string part = node(*n);
            if (!part.empty()) parts.push_back(part);
        }
        std::string out;
        for (const auto& p : parts) out += p;
        return out;
    }

    std::string list_of(const NodeList& l) { return seq(l); }
    std::string slot(const NodeList& l) { return syn_.row(seq(l)); }

    /* ---- individual nodes ---------------------------------------------- */

    std::string script(const ScriptNode& s, const std::string& base) {
        if (!s.hasSub && !s.hasSup) return base;
        return syn_.script(syn_.row(base), slot(s.sub), slot(s.sup),
                           s.hasSub, s.hasSup, base_takes_limits(s.base));
    }

    /* Whether the thing being scripted is a function name that wants its
     * scripts UNDER and OVER: \lim, \max, \sup and their family.  Read off
     * the tree the same way the layout reads it, so a picture and a paste
     * cannot disagree about where the limit goes. */
    static bool base_takes_limits(const NodeList& base) {
        return mtef_name_takes_limits(function_name_of(base));
    }

    std::string nary_body(const Node& n) {
        if (n.tag() == Node::kIntegral)
            return list_of(static_cast<const IntegralNode&>(n).body);
        return list_of(static_cast<const BigOpNode&>(n).body);
    }

    static bool nary_limits_empty(const Node& n) {
        if (n.tag() == Node::kIntegral) {
            const auto& i = static_cast<const IntegralNode&>(n);
            return !i.displayLower && !i.displayUpper && !i.hasLower && !i.hasUpper;
        }
        const auto& b = static_cast<const BigOpNode&>(n);
        return !b.displayLower && !b.displayUpper && !b.hasLower && !b.hasUpper;
    }

    /* Limits the passes recovered win over the template's own slots, which
     * EQNEDT32 usually leaves empty -- and those slots count only when the
     * MTEF variation bits say they are present. */
    std::string nary(const Node& n, const std::string& body,
                     const ScriptNode* limits) {
        int selector; bool stacked;
        std::string lower, upper;
        if (n.tag() == Node::kIntegral) {
            const auto& i = static_cast<const IntegralNode&>(n);
            selector = i.selector; stacked = i.hasLimits;
            lower = i.displayLower ? line(*i.displayLower)
                  : i.hasLower    ? list_of(i.lower) : std::string();
            upper = i.displayUpper ? line(*i.displayUpper)
                  : i.hasUpper    ? list_of(i.upper) : std::string();
        } else {
            const auto& b = static_cast<const BigOpNode&>(n);
            selector = b.selector; stacked = b.hasLimits;
            lower = b.displayLower ? line(*b.displayLower)
                  : b.hasLower    ? list_of(b.lower) : std::string();
            upper = b.displayUpper ? line(*b.displayUpper)
                  : b.hasUpper    ? list_of(b.upper) : std::string();
        }
        if (limits) {
            if (limits->hasSub) lower = list_of(limits->sub);
            if (limits->hasSup) upper = list_of(limits->sup);
        }
        return syn_.nary(nary_char(selector), syn_.row(lower), syn_.row(upper),
                         syn_.row(body), stacked,
                         !lower.empty(), !upper.empty());
    }

    std::string node(const Node& n) {
        switch (n.tag()) {
            case Node::kLine:
                return line(static_cast<const LineNode&>(n));

            case Node::kSize:
            case Node::kFont:
                return std::string();   /* Office sizes the equation itself */

            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(n);
                uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));
                if (!cp) return std::string();
                /* OMML/RTF/MathML do not carry this tree's logical extended
                 * typeface codes.  Emit the corresponding Unicode math
                 * alphabet scalar, as the SVG renderer does. */
                if (c.typeface == TF_MATH_SANS)
                    cp = mtef_math_sans_of(cp) ? mtef_math_sans_of(cp) : cp;
                else if (c.typeface == TF_MATH_MONO)
                    cp = mtef_math_mono_of(cp) ? mtef_math_mono_of(cp) : cp;
                else if (c.typeface == TF_MATH_FRAKTUR)
                    cp = mtef_math_fraktur_of(cp) ? mtef_math_fraktur_of(cp) : cp;
                return syn_.run(utf8_of(cp), run_style(c.typeface));
            }
            case Node::kScript: {
                const auto& s = static_cast<const ScriptNode&>(n);
                return script(s, list_of(s.base));
            }
            case Node::kFrac: {
                const auto& f = static_cast<const FracNode&>(n);
                return syn_.fraction(slot(f.numer), slot(f.denom), f.slashed);
            }
            case Node::kSqrt: {
                const auto& s = static_cast<const SqrtNode&>(n);
                return syn_.radical(slot(s.content), slot(s.index), s.hasIndex);
            }
            case Node::kFence: {
                const auto& f = static_cast<const FenceNode&>(n);
                Fence ch = fence_chars(f.selector);
                /* variation 1 = left delimiter only, 2 = right only */
                const char* beg = (f.variation == 2) ? "" : ch.beg;
                const char* end = (f.variation == 1) ? "" : ch.end;
                return syn_.fence(slot(f.content), beg, end);
            }
            case Node::kIntegral:
            case Node::kBigOp:
                return nary(n, std::string(), nullptr);

            case Node::kMatrix: {
                const auto& m = static_cast<const MatrixNode&>(n);
                std::vector<std::string> cells;
                for (int r = 0; r < m.rows; ++r)
                    for (int c = 0; c < m.cols; ++c) {
                        size_t idx = size_t(r) * size_t(m.cols) + size_t(c);
                        std::string cell;
                        if (idx < m.elements.size() && m.elements[idx])
                            cell = node(*m.elements[idx]);
                        cells.push_back(syn_.row(cell));
                    }
                return syn_.matrix(m.rows, m.cols, cells);
            }
            case Node::kPile: {
                const auto& p = static_cast<const PileNode&>(n);
                std::vector<std::string> lines;
                for (const auto& ln : p.lines)
                    lines.push_back(syn_.row(ln ? node(*ln) : std::string()));
                return syn_.stack(lines);
            }
            case Node::kEmbell: {
                const auto& e = static_cast<const EmbellNode&>(n);
                return syn_.accent(accent_char(e.embellType), slot(e.content));
            }
            case Node::kDecoration: {
                const auto& d = static_cast<const DecorationNode&>(n);
                return syn_.bar(slot(d.content), d.selector == tmOBAR);
            }
            case Node::kBraceDeco:
                return list_of(static_cast<const BraceDecoNode&>(n).content);

            default:
                return std::string();
        }
    }
};

}  // namespace

std::string write_math(const LineNode& root, const MathSyntax& syntax,
                       bool display, bool run_passes) {
    Walker w(syntax, run_passes);
    return w.run(root, display);
}

}  // namespace mtef
