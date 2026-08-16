/*
 * mtef_omml.cpp -- node tree -> OMML
 *
 * Structural mapping only; Office performs the typesetting:
 *
 *   LineNode      -> the child sequence of <m:oMath>
 *   CharNode      -> <m:r><m:t>c</m:t></m:r>   (with m:rPr/m:sty for upright)
 *   ScriptNode    -> <m:sSub> / <m:sSup> / <m:sSubSup>
 *   FracNode      -> <m:f><m:num/><m:den/></m:f>
 *   SqrtNode      -> <m:rad>  (degHide for a square root)
 *   FenceNode     -> <m:d>    (begChr/endChr for anything but parentheses)
 *   BigOp/Integral-> <m:nary> (m:chr = the operator, subHide/supHide)
 *   MatrixNode    -> <m:m>
 *   PileNode      -> <m:eqArr>
 *
 * Two things make this more than a tree walk, both consequences of how
 * EQNEDT32 lays a LINE out:
 *
 *   1. The same LINE passes the LaTeX emitter runs must run here.  MTEF stores
 *      a fence or a big operator as an empty template followed by sibling
 *      content, and only the passes put that content back in the slot.
 *
 *   2. A script template carries no base: MTEF writes the base character as the
 *      *preceding sibling*.  LaTeX does not care (`x` then `_{c}` binds), but
 *      an OMML <m:e> left empty renders as an empty box in Word, so the base is
 *      absorbed from the part emitted just before it.
 */
#include "mtef_omml.h"
#include "mtef_parser.h"
#include "line_pass.h"
#include "tex_parser.h"

#include <sstream>
#include <string>
#include <vector>

namespace mtef {
namespace {

std::string utf8_of(uint32_t cp) {
    std::string s;
    if (cp < 0x80) {
        s += char(cp);
    } else if (cp < 0x800) {
        s += char(0xC0 | (cp >> 6));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xE0 | (cp >> 12));
        s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    }
    return s;
}

std::string esc(const std::string& s) {
    std::string o;
    for (char c : s) {
        switch (c) {
            case '&': o += "&amp;"; break;
            case '<': o += "&lt;"; break;
            case '>': o += "&gt;"; break;
            case '"': o += "&quot;"; break;
            default: o += c;
        }
    }
    return o;
}

/* Fence delimiters; OMML defaults to parentheses, so only other shapes need
 * begChr/endChr. */
struct Fence { const char* beg; const char* end; };
Fence fence_chars(int selector) {
    switch (selector) {
        case tmPAREN: return {nullptr, nullptr};
        case tmBRACK: return {"[", "]"};
        case tmBRACE: return {"{", "}"};
        case tmANGLE: return {"\xE2\x9F\xA8", "\xE2\x9F\xA9"};   /* U+27E8/9 */
        case tmBAR:   return {"|", "|"};
        case tmDBAR:  return {"\xE2\x80\x96", "\xE2\x80\x96"};   /* U+2016 */
        case tmFLOOR: return {"\xE2\x8C\x8A", "\xE2\x8C\x8B"};   /* U+230A/B */
        case tmCEIL:  return {"\xE2\x8C\x88", "\xE2\x8C\x89"};   /* U+2308/9 */
        default:      return {nullptr, nullptr};
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

/* MTEF embellishment code -> the combining accent OMML draws.  0 leaves the
 * OMML default (a circumflex). */
uint32_t accent_char(int embellType) {
    switch (embellType) {
        case 2:  return 0x0307;   /* dot            */
        case 3:  return 0x0308;   /* double dot     */
        case 4:  return 0x20DB;   /* triple dot     */
        case 8:  return 0x0303;   /* tilde          */
        case 9:  return 0x0302;   /* hat            */
        case 11: return 0x20D7;   /* right arrow    */
        case 12: return 0x20D6;   /* left arrow     */
        case 13: return 0x20E1;   /* both arrows    */
        case 17: return 0x0304;   /* bar            */
        default: return 0;
    }
}

/* MTEF marks style by typeface; OMML marks it per run.  Variables and Greek
 * letters are italic, everything else upright. */
bool upright_typeface(int tf) {
    return !(tf == TF_VARIABLE || tf == TF_LCGREEK || tf == TF_UCGREEK ||
             tf == TF_VECTOR);
}

class OmmlWriter {
public:
    explicit OmmlWriter(const OmmlOptions& opt) : opt_(opt) {}

    std::string run(const LineNode& root, int prodVer, bool run_passes) {
        run_passes_ = run_passes;
        prodVer_ = prodVer;
        depth_ = 0;
        std::ostringstream o;
        const char* ns = opt_.declare_namespace
            ? " xmlns:m=\"http://schemas.openxmlformats.org/officeDocument/2006/math\""
            : "";
        if (opt_.display) o << "<m:oMathPara" << ns << "><m:oMath>";
        else              o << "<m:oMath" << ns << ">";
        o << line_xml(root);
        if (opt_.display) o << "</m:oMath></m:oMathPara>";
        else              o << "</m:oMath>";
        return o.str();
    }

private:
    const OmmlOptions& opt_;
    int prodVer_ = 3;
    int depth_ = 0;
    bool run_passes_ = true;

    /* ---- sequence assembly -------------------------------------------- */

    /* Render one LINE: run the passes, then assemble its children.  The passes
     * mutate the child list (they move sibling content into empty slots and
     * drop the siblings), which is why the list is taken as mutable. */
    std::string line_xml(const LineNode& line) {
        if (line.isNull) return std::string();
        NodeList& children = const_cast<NodeList&>(line.children);
        if (run_passes_) pipeline_.process(children, depth_, prodVer_);
        depth_++;
        std::string s = seq_xml(children);
        depth_--;
        return s;
    }

    /* Assemble a sibling sequence, giving base-less scripts the part before
     * them and big operators everything after them. */
    std::string seq_xml(const NodeList& list, size_t start = 0) {
        std::vector<std::string> parts;
        for (size_t i = start; i < list.size(); ++i) {
            const Node* n = list[i].get();
            if (!n) continue;

            if (n->tag() == Node::kScript) {
                const auto& s = static_cast<const ScriptNode&>(*n);
                if (s.base.empty() && !parts.empty()) {
                    std::string base = parts.back();
                    parts.pop_back();
                    parts.push_back(script_xml(s, base));
                    continue;
                }
            }

            /* A big operator's operand is normally in its own body slot.  When
             * the passes could not fill it, the operand is still the sibling
             * run that follows, and it has to be absorbed -- an empty m:e draws
             * a placeholder box next to the operator. */
            if (n->tag() == Node::kIntegral || n->tag() == Node::kBigOp) {
                /* EQNEDT32 writes limits the template itself does not carry as
                 * a base-less script right after the operator.  Left alone that
                 * becomes a script with an empty base -- a placeholder box next
                 * to the operator instead of its limits. */
                const ScriptNode* limits = nullptr;
                if (nary_limits_empty(*n) && i + 1 < list.size() &&
                    list[i + 1] && list[i + 1]->tag() == Node::kScript) {
                    const auto& s = static_cast<const ScriptNode&>(*list[i + 1]);
                    if (s.base.empty()) { limits = &s; ++i; }
                }
                std::string body = nary_body(*n);
                if (body.empty()) {
                    parts.push_back(nary_xml(*n, seq_xml(list, i + 1), limits));
                    break;
                }
                parts.push_back(nary_xml(*n, body, limits));
                continue;
            }

            std::string part = node_xml(*n);
            if (!part.empty()) parts.push_back(part);
        }
        std::string out;
        for (const auto& p : parts) out += p;
        return out;
    }

    std::string list_xml(const NodeList& list) { return seq_xml(list); }

    /* An OMML slot always expects element content; an empty slot must still be
     * present or Word treats the equation as malformed. */
    std::string slot(const char* tag, const NodeList& list) {
        return std::string("<") + tag + '>' + list_xml(list) + "</" + tag + '>';
    }
    std::string slot_raw(const char* tag, const std::string& xml) {
        return std::string("<") + tag + '>' + xml + "</" + tag + '>';
    }

    /* ---- individual nodes ---------------------------------------------- */

    std::string script_xml(const ScriptNode& s, const std::string& base_xml) {
        std::ostringstream o;
        if (s.hasSub && s.hasSup) {
            o << "<m:sSubSup><m:sSubSupPr/>" << slot_raw("m:e", base_xml)
              << slot("m:sub", s.sub) << slot("m:sup", s.sup) << "</m:sSubSup>";
        } else if (s.hasSub) {
            o << "<m:sSub><m:sSubPr/>" << slot_raw("m:e", base_xml)
              << slot("m:sub", s.sub) << "</m:sSub>";
        } else if (s.hasSup) {
            o << "<m:sSup><m:sSupPr/>" << slot_raw("m:e", base_xml)
              << slot("m:sup", s.sup) << "</m:sSup>";
        } else {
            o << base_xml;
        }
        return o.str();
    }

    std::string nary_body(const Node& n) {
        if (n.tag() == Node::kIntegral)
            return list_xml(static_cast<const IntegralNode&>(n).body);
        return list_xml(static_cast<const BigOpNode&>(n).body);
    }

    static bool nary_limits_empty(const Node& n) {
        if (n.tag() == Node::kIntegral) {
            const auto& i = static_cast<const IntegralNode&>(n);
            return !i.displayLower && !i.displayUpper && !i.hasLower && !i.hasUpper;
        }
        const auto& b = static_cast<const BigOpNode&>(n);
        return !b.displayLower && !b.displayUpper && !b.hasLower && !b.hasUpper;
    }

    /* Big operators and integrals share <m:nary>.  The limits the passes
     * recovered (displayLower/displayUpper) win over the template's own slots,
     * which EQNEDT32 usually leaves empty -- and those slots count only when
     * the MTEF variation bits say they are present. */
    std::string nary_xml(const Node& n, const std::string& body_xml,
                         const ScriptNode* limits = nullptr) {
        int selector; bool stacked;
        std::string lower, upper;
        if (n.tag() == Node::kIntegral) {
            const auto& i = static_cast<const IntegralNode&>(n);
            selector = i.selector; stacked = i.hasLimits;
            lower = i.displayLower ? line_xml(*i.displayLower)
                  : i.hasLower    ? list_xml(i.lower) : std::string();
            upper = i.displayUpper ? line_xml(*i.displayUpper)
                  : i.hasUpper    ? list_xml(i.upper) : std::string();
        } else {
            const auto& b = static_cast<const BigOpNode&>(n);
            selector = b.selector; stacked = b.hasLimits;
            lower = b.displayLower ? line_xml(*b.displayLower)
                  : b.hasLower    ? list_xml(b.lower) : std::string();
            upper = b.displayUpper ? line_xml(*b.displayUpper)
                  : b.hasUpper    ? list_xml(b.upper) : std::string();
        }
        if (limits) {
            if (limits->hasSub) lower = list_xml(limits->sub);
            if (limits->hasSup) upper = list_xml(limits->sup);
        }

        std::ostringstream o;
        o << "<m:nary><m:naryPr><m:chr m:val=\""
          << esc(utf8_of(nary_char(selector))) << "\"/>"
          << "<m:limLoc m:val=\"" << (stacked ? "undOvr" : "subSup") << "\"/>";
        if (lower.empty()) o << "<m:subHide m:val=\"1\"/>";
        if (upper.empty()) o << "<m:supHide m:val=\"1\"/>";
        o << "</m:naryPr>"
          << slot_raw("m:sub", lower) << slot_raw("m:sup", upper)
          << slot_raw("m:e", body_xml) << "</m:nary>";
        return o.str();
    }

    std::string text_run(const std::string& utf8, int typeface) {
        const char* sty = nullptr;
        if (typeface == TF_VECTOR)            sty = "bi";   /* bold italic */
        else if (upright_typeface(typeface))  sty = "p";    /* upright     */
        std::string s = "<m:r>";
        if (sty && opt_.italic_variables)
            s += std::string("<m:rPr><m:sty m:val=\"") + sty + "\"/></m:rPr>";
        s += "<m:t>" + esc(utf8) + "</m:t></m:r>";
        return s;
    }

    std::string node_xml(const Node& n) {
        std::ostringstream o;
        switch (n.tag()) {
            case Node::kLine:
                o << line_xml(static_cast<const LineNode&>(n));
                break;

            case Node::kSize:
            case Node::kFont:
                break;   /* Office sizes the equation; MTEF size records drop out */

            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(n);
                uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));
                if (!cp) break;
                o << text_run(utf8_of(cp), c.typeface);
                break;
            }
            case Node::kScript: {
                const auto& s = static_cast<const ScriptNode&>(n);
                o << script_xml(s, list_xml(s.base));
                break;
            }
            case Node::kFrac: {
                const auto& f = static_cast<const FracNode&>(n);
                o << "<m:f><m:fPr>";
                if (f.slashed) o << "<m:type m:val=\"skw\"/>";
                o << "</m:fPr>" << slot("m:num", f.numer)
                  << slot("m:den", f.denom) << "</m:f>";
                break;
            }
            case Node::kSqrt: {
                const auto& s = static_cast<const SqrtNode&>(n);
                o << "<m:rad><m:radPr>";
                if (!s.hasIndex) o << "<m:degHide m:val=\"1\"/>";
                o << "</m:radPr>";
                if (s.hasIndex) o << slot("m:deg", s.index);
                else            o << "<m:deg/>";
                o << slot("m:e", s.content) << "</m:rad>";
                break;
            }
            case Node::kFence: {
                const auto& f = static_cast<const FenceNode&>(n);
                Fence ch = fence_chars(f.selector);
                o << "<m:d><m:dPr>";
                /* variation 1 = left delimiter only, 2 = right only */
                if (f.variation == 2)      o << "<m:begChr m:val=\"\"/>";
                else if (ch.beg)           o << "<m:begChr m:val=\"" << esc(ch.beg) << "\"/>";
                if (f.variation == 1)      o << "<m:endChr m:val=\"\"/>";
                else if (ch.end)           o << "<m:endChr m:val=\"" << esc(ch.end) << "\"/>";
                o << "</m:dPr>" << slot("m:e", f.content) << "</m:d>";
                break;
            }
            case Node::kIntegral:
            case Node::kBigOp:
                o << nary_xml(n, std::string());
                break;

            case Node::kMatrix: {
                const auto& m = static_cast<const MatrixNode&>(n);
                o << "<m:m><m:mPr><m:mcs><m:mc><m:mcPr><m:count m:val=\""
                  << m.cols << "\"/><m:mcJc m:val=\"center\"/></m:mcPr></m:mc>"
                  << "</m:mcs></m:mPr>";
                for (int r = 0; r < m.rows; ++r) {
                    o << "<m:mr>";
                    for (int c = 0; c < m.cols; ++c) {
                        size_t idx = size_t(r) * size_t(m.cols) + size_t(c);
                        o << "<m:e>";
                        if (idx < m.elements.size() && m.elements[idx])
                            o << node_xml(*m.elements[idx]);
                        o << "</m:e>";
                    }
                    o << "</m:mr>";
                }
                o << "</m:m>";
                break;
            }
            case Node::kPile: {
                const auto& p = static_cast<const PileNode&>(n);
                o << "<m:eqArr><m:eqArrPr/>";
                for (const auto& ln : p.lines) {
                    o << "<m:e>";
                    if (ln) o << node_xml(*ln);
                    o << "</m:e>";
                }
                o << "</m:eqArr>";
                break;
            }
            case Node::kEmbell: {
                const auto& e = static_cast<const EmbellNode&>(n);
                o << "<m:acc><m:accPr>";
                if (uint32_t g = accent_char(e.embellType))
                    o << "<m:chr m:val=\"" << esc(utf8_of(g)) << "\"/>";
                o << "</m:accPr>" << slot("m:e", e.content) << "</m:acc>";
                break;
            }
            case Node::kDecoration: {
                const auto& d = static_cast<const DecorationNode&>(n);
                const bool over = (d.selector == tmOBAR);
                o << "<m:bar><m:barPr><m:pos m:val=\""
                  << (over ? "top" : "bot") << "\"/></m:barPr>"
                  << slot("m:e", d.content) << "</m:bar>";
                break;
            }
            case Node::kBraceDeco: {
                const auto& b = static_cast<const BraceDecoNode&>(n);
                o << list_xml(b.content);
                break;
            }
            default:
                break;
        }
        return o.str();
    }

    PassPipeline pipeline_;
};

}  // namespace

std::string render_omml(const LineNode& root, const OmmlOptions& opt,
                        bool run_passes) {
    OmmlWriter w(opt);
    return w.run(root, 3, run_passes);
}

std::string mtef_to_omml(const uint8_t* data, size_t len, const OmmlOptions& opt) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return std::string();
    OmmlWriter w(opt);
    return w.run(*res.root, res.prodVer, /*run_passes=*/true);
}

std::string tex_to_omml(const std::string& latex, const OmmlOptions& opt) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_omml(*root, opt, /*run_passes=*/false);
}

}  // namespace mtef
