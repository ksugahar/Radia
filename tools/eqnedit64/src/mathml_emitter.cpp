#include "mathml_emitter.h"

#include "tex_parser.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <memory>
#include <sstream>

namespace eqnedit {
namespace {

std::string utf8_of(uint32_t cp) {
    std::string out;
    if (cp <= 0x7F) {
        out.push_back(char(cp));
    } else if (cp <= 0x7FF) {
        out.push_back(char(0xC0 | (cp >> 6)));
        out.push_back(char(0x80 | (cp & 0x3F)));
    } else if (cp <= 0xFFFF) {
        out.push_back(char(0xE0 | (cp >> 12)));
        out.push_back(char(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(char(0x80 | (cp & 0x3F)));
    } else if (cp <= 0x10FFFF) {
        out.push_back(char(0xF0 | (cp >> 18)));
        out.push_back(char(0x80 | ((cp >> 12) & 0x3F)));
        out.push_back(char(0x80 | ((cp >> 6) & 0x3F)));
        out.push_back(char(0x80 | (cp & 0x3F)));
    }
    return out;
}

std::string xml_text(const std::string& value) {
    std::string out;
    out.reserve(value.size());
    for (char c : value) {
        switch (c) {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        default: out += c; break;
        }
    }
    return out;
}

std::string element(const char* name, const std::string& content,
                    const char* attributes = nullptr) {
    std::string out = "<";
    out += name;
    if (attributes && *attributes) {
        out += " ";
        out += attributes;
    }
    out += ">";
    out += content;
    out += "</";
    out += name;
    out += ">";
    return out;
}

std::string row(const std::string& content) {
    return content.empty() ? "<mrow/>" : element("mrow", content);
}

class MathMlEmitter {
public:
    std::string emit(const LineNode& root, double pointSize) {
        if (!std::isfinite(pointSize) || pointSize <= 0.0) pointSize = 24.0;
        std::ostringstream size;
        size << std::fixed << std::setprecision(3) << pointSize;
        std::string digits = size.str();
        while (!digits.empty() && digits.back() == '0') digits.pop_back();
        if (!digits.empty() && digits.back() == '.') digits.pop_back();
        /* The browser editor publishes inline MathML to Office.  Keep the
         * native root identical so PowerPoint does not turn one product into
         * a centred equation paragraph and the other into inline Office Math. */
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" "
               "display=\"inline\" mathsize=\"" + digits + "pt\">" +
               sequence_content(root.children) + "</math>";
    }

private:
    std::string sequence_content(const NodeList& nodes) {
        std::string content;
        for (const auto& node : nodes)
            if (node) content += emit_node(*node);
        return content;
    }

    std::string sequence(const NodeList& nodes) {
        return row(sequence_content(nodes));
    }

    std::string token(const CharNode& ch) {
        /* Explicit spacing has a width and no character, which is exactly
         * what <mspace> is for; an empty <mi> would collapse in Word. */
        if (ch.typeface == TF_SPACE) {
            std::ostringstream width;
            width << std::fixed << std::setprecision(4)
                  << space_width_em(ch.latex.c_str());
            return "<mspace width=\"" + width.str() + "em\"/>";
        }
        const std::string value = xml_text(utf8_of(ch.charCode));
        std::string result;
        if (ch.typeface == TF_TEXT) {
            result = element("mtext", value);
        } else if (ch.typeface == TF_NUMBER) {
            result = element("mn", value);
        } else if (ch.typeface == TF_FUNCTION || ch.typeface == TF_ROMAN) {
            result = element("mi", value, "mathvariant=\"normal\"");
        } else if (ch.typeface == TF_MATH_ITALIC) {
            result = element("mi", value, "mathvariant=\"italic\"");
        } else if (ch.typeface == TF_VECTOR ||
                   ch.typeface == TF_BOLD_SYMBOL) {
            result = element("mi", value, "mathvariant=\"bold\"");
        } else if (ch.typeface == TF_MATH_SANS) {
            result = element("mi", value, "mathvariant=\"sans-serif\"");
        } else if (ch.typeface == TF_MATH_MONO) {
            result = element("mi", value, "mathvariant=\"monospace\"");
        } else if (ch.typeface == TF_MATH_SCRIPT) {
            result = element("mi", value, "mathvariant=\"script\"");
        } else if (ch.typeface == TF_MATH_DOUBLE) {
            result = element("mi", value, "mathvariant=\"double-struck\"");
        } else if (ch.typeface == TF_MATH_FRAKTUR) {
            result = element("mi", value, "mathvariant=\"fraktur\"");
        } else if (ch.typeface == TF_SYMBOL ||
                   (ch.charCode < 128 && !std::isalnum(int(ch.charCode)))) {
            result = element("mo", value);
        } else {
            result = element("mi", value);
        }
        for (Embellishment embell : ch.embells) result = apply_embell(result, embell);
        return result;
    }

    std::string apply_embell(const std::string& base, Embellishment embell) {
        const char* mark = nullptr;
        bool under = false;
        switch (embell) {
        case EM_NONE: return base;   /* the field's default */
        case EM_DOT: mark = "&#x02D9;"; break;
        case EM_DDOT: mark = "&#x00A8;"; break;
        case EM_TDOT: mark = "&#x20DB;"; break;
        case EM_TILDE: mark = "~"; break;
        case EM_HAT: mark = "^"; break;
        case EM_NOT: mark = "&#x0338;"; break;
        case EM_RARROW: mark = "&#x2192;"; break;
        case EM_LARROW: mark = "&#x2190;"; break;
        case EM_BARROW: mark = "&#x2194;"; break;
        case EM_R1ARROW: mark = "&#x21C0;"; break;
        case EM_L1ARROW: mark = "&#x21BC;"; break;
        case EM_MBAR: mark = "_"; under = true; break;
        case EM_OBAR: mark = "&#x00AF;"; break;
        case EM_FROWN: mark = "&#x2322;"; break;
        case EM_SMILE: mark = "&#x2323;"; under = true; break;
        /* The prime family is a suffix, not something to put over the base. */
        case EM_PRIME: return base + element("mo", "&#x2032;");
        case EM_DPRIME: return base + element("mo", "&#x2033;");
        case EM_BPRIME:
        case EM_TPRIME: return base + element("mo", "&#x2034;");
        /* No default, deliberately: /W4 /WX turns a new embellishment that
         * nobody mapped into C4062 at compile time. `default: return base`
         * dropped the decoration silently instead, so an equation pasted into
         * Office quietly lost it. */
        }
        return element(under ? "munder" : "mover",
                       base + element("mo", mark, "stretchy=\"true\""));
    }

    std::string script(const ScriptNode& value) {
        const std::string base = sequence(value.base);
        if (value.hasSub && value.hasSup)
            return element("msubsup", base + sequence(value.sub) +
                                      sequence(value.sup));
        if (value.hasSub)
            return element("msub", base + sequence(value.sub));
        if (value.hasSup)
            return element("msup", base + sequence(value.sup));
        return base;
    }

    std::string limited_operator(const std::string& symbol,
                                 const NodeList& lower, bool hasLower,
                                 const NodeList& upper, bool hasUpper) {
        const std::string op = element("mo", symbol,
            "largeop=\"true\" movablelimits=\"true\"");
        if (hasLower && hasUpper)
            return element("munderover", op + sequence(lower) + sequence(upper));
        if (hasLower) return element("munder", op + sequence(lower));
        if (hasUpper) return element("mover", op + sequence(upper));
        return op;
    }

    std::string integral_operator(const std::string& symbol,
                                  const NodeList& lower, bool hasLower,
                                  const NodeList& upper, bool hasUpper) {
        const std::string op = element("mo", symbol,
            "largeop=\"true\" movablelimits=\"true\"");
        /* TeX (and MathJax.tex2mml) puts ordinary integral limits beside the
         * glyph.  munderover is correct for sums/products, but using it for
         * integrals was the visible native/Web PowerPoint mismatch. */
        if (hasLower && hasUpper)
            return element("msubsup", op + sequence(lower) + sequence(upper));
        if (hasLower) return element("msub", op + sequence(lower));
        if (hasUpper) return element("msup", op + sequence(upper));
        return op;
    }

    std::string table(const MatrixNode& matrix) {
        const int rowCount = std::max(0, matrix.rows);
        const int columnCount = std::max(1, matrix.cols);
        std::string rows;
        for (int r = 0; r < rowCount; ++r) {
            std::string cells;
            for (int c = 0; c < columnCount; ++c) {
                const size_t at = size_t(r * columnCount + c);
                const std::string cell = at < matrix.elements.size() &&
                        matrix.elements[at]
                    ? emit_node(*matrix.elements[at]) : "<mrow/>";
                cells += element("mtd", cell);
            }
            rows += element("mtr", cells);
        }
        const char* alignment =
            matrix.layoutKind == MatrixNode::kAlignedLayout
                ? (matrix.cols <= 1 ? "columnalign=\"left\""
                                    : "columnalign=\"right left\"")
            : matrix.layoutKind == MatrixNode::kCasesLayout
                ? "columnalign=\"left\""
                : nullptr;
        return element("mtable", rows, alignment);
    }

    std::string pile(const PileNode& value) {
        std::string rows;
        for (const auto& line : value.lines)
            rows += element("mtr", element("mtd",
                line ? emit_node(*line) : "<mrow/>"));
        return element("mtable", rows,
            value.halign == 1 || value.kind == 1
                ? "columnalign=\"left\"" : "columnalign=\"center\"");
    }

    std::string emit_node(const Node& node) {
        switch (node.tag()) {
        case Node::kLine:
            return sequence(static_cast<const LineNode&>(node).children);
        case Node::kChar:
            return token(static_cast<const CharNode&>(node));
        case Node::kFrac: {
            const auto& value = static_cast<const FracNode&>(node);
            return element("mfrac", sequence(value.numer) + sequence(value.denom),
                           value.slashed ? "bevelled=\"true\"" : nullptr);
        }
        case Node::kSqrt: {
            const auto& value = static_cast<const SqrtNode&>(node);
            return value.hasIndex
                ? element("mroot", sequence(value.content) + sequence(value.index))
                : element("msqrt", sequence(value.content));
        }
        case Node::kScript:
            return script(static_cast<const ScriptNode&>(node));
        case Node::kFence: {
            const auto& value = static_cast<const FenceNode&>(node);
            static const char* left[] = {
                "&#x27E8;", "(", "{", "[", "|", "&#x2016;", "&#x230A;",
                "&#x2308;", "[", "]", "]", "[", "("
            };
            static const char* right[] = {
                "&#x27E9;", ")", "}", "]", "|", "&#x2016;", "&#x230B;",
                "&#x2309;", "[", "]", "[", ")", "]"
            };
            auto clamp = [](int s) {
                return s >= 0 && s <= 12 ? s : int(tmPAREN);
            };
            const int selector = clamp(value.selector);
            const int closing = clamp(value.right_selector());
            std::string content;
            if (value.variation != 2)
                content += element("mo", left[selector],
                    "fence=\"true\" stretchy=\"true\"");
            content += sequence(value.content);
            if (value.variation != 1)
                content += element("mo", right[closing],
                    "fence=\"true\" stretchy=\"true\"");
            return row(content);
        }
        case Node::kIntegral: {
            const auto& value = static_cast<const IntegralNode&>(node);
            const char* symbol = "&#x222B;";
            if (value.selector == tmDINT) symbol = "&#x222C;";
            else if (value.selector == tmTINT) symbol = "&#x222D;";
            else if (value.selector == tmSSINT) symbol = "&#x222E;";
            else if (value.selector == tmDSINT) symbol = "&#x222F;";
            else if (value.selector == tmTSINT) symbol = "&#x2230;";
            return integral_operator(symbol, value.lower, value.hasLower,
                value.upper, value.hasUpper) + sequence_content(value.body);
        }
        case Node::kBigOp: {
            const auto& value = static_cast<const BigOpNode&>(node);
            const char* symbol = "&#x2211;";
            if (value.selector == tmPROD || value.selector == tmIPROD)
                symbol = "&#x220F;";
            else if (value.selector == tmCOPROD || value.selector == tmICOPROD)
                symbol = "&#x2210;";
            else if (value.selector == tmUNION || value.selector == tmIUNION)
                symbol = "&#x22C3;";
            else if (value.selector == tmINTER || value.selector == tmIINTER)
                symbol = "&#x22C2;";
            return limited_operator(symbol, value.lower, value.hasLower,
                value.upper, value.hasUpper) + sequence_content(value.body);
        }
        case Node::kDecoration: {
            const auto& value = static_cast<const DecorationNode&>(node);
            bool under = value.selector == tmUBAR;
            const char* mark = under ? "_" : "&#x00AF;";
            if (value.selector == tmRARROW) mark = "&#x2192;";
            else if (value.selector == tmLARROW) mark = "&#x2190;";
            else if (value.selector == tmBARROW) mark = "&#x2194;";
            return element(under ? "munder" : "mover", sequence(value.content) +
                element("mo", mark, "stretchy=\"true\""));
        }
        case Node::kEmbell: {
            const auto& value = static_cast<const EmbellNode&>(node);
            return apply_embell(sequence(value.content), value.embellType);
        }
        case Node::kBraceDeco: {
            const auto& value = static_cast<const BraceDecoNode&>(node);
            const bool over = value.selector == tmUHBRACE;
            const std::string braced = element(over ? "mover" : "munder",
                sequence(value.content) + element("mo", over ? "&#x23DE;" :
                    "&#x23DF;", "stretchy=\"true\""));
            return element(over ? "mover" : "munder",
                           braced + sequence(value.label));
        }
        case Node::kDirac: {
            const auto& value = static_cast<const DiracNode&>(node);
            const std::string left = element("mo", "&#x27E8;",
                "fence=\"true\" stretchy=\"true\"");
            const std::string middle = element("mo", "|", "stretchy=\"true\"");
            if (value.variation != 0)
                return row(left + sequence(value.bra) + middle);
            return row(left + sequence(value.bra) + middle + sequence(value.ket) +
                element("mo", "&#x27E9;", "fence=\"true\" stretchy=\"true\""));
        }
        case Node::kLim: {
            const auto& value = static_cast<const LimNode&>(node);
            const std::string op = element("mi", "lim", "mathvariant=\"normal\"");
            return value.content.empty() ? op
                : element("munder", op + sequence(value.content));
        }
        case Node::kPile:
            return pile(static_cast<const PileNode&>(node));
        case Node::kMatrix:
            return table(static_cast<const MatrixNode&>(node));
        case Node::kFunction: {
            const auto& value = static_cast<const FunctionNode&>(node);
            return element("mi", xml_text(value.name), "mathvariant=\"normal\"");
        }
        case Node::kText:
            return element("mtext", xml_text(static_cast<const TextNode&>(node).text));
        case Node::kMathbf:
            return element("mstyle", sequence(static_cast<const MathbfNode&>(node).content),
                           "mathvariant=\"bold\"");
        case Node::kGroup:
            return sequence(static_cast<const GroupNode&>(node).children);
        case Node::kPrime: {
            const int count = std::max(1, static_cast<const PrimeNode&>(node).count);
            std::string primes;
            for (int i = 0; i < count; ++i) primes += "&#x2032;";
            return element("mo", primes);
        }
        case Node::kDegree:
            return element("mo", "&#x00B0;");
        case Node::kOverset: {
            const auto& value = static_cast<const OversetNode&>(node);
            return element(value.under ? "munder" : "mover",
                           sequence(value.base) + sequence(value.over));
        }
        case Node::kRM: {
            const char ch = static_cast<const RMNode&>(node).ch;
            return element("mtext", xml_text(std::string(1, ch)));
        }
        case Node::kSize:
        case Node::kFont:
            return {};
        }
        return {};
    }
};

}  // namespace

std::string tree_to_mathml(const LineNode& root, double pointSize) {
    return MathMlEmitter().emit(root, pointSize);
}

std::string latex_to_mathml(const std::string& latex, double pointSize) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    return root ? tree_to_mathml(*root, pointSize) : std::string();
}

namespace {

bool collect_unanchored_aligned_rows(
        const LineNode& line, std::vector<const LineNode*>& rows) {
    if (line.children.size() == 1 && line.children[0] &&
        line.children[0]->tag() == Node::kMatrix) {
        const auto& matrix =
            static_cast<const MatrixNode&>(*line.children[0]);
        if (matrix.layoutKind == MatrixNode::kAlignedLayout &&
            matrix.cols <= 1 && matrix.rows > 0 &&
            matrix.elements.size() >= size_t(matrix.rows)) {
            for (int rowIndex = 0; rowIndex < matrix.rows; ++rowIndex) {
                const Node* rowNode = matrix.elements[size_t(rowIndex)].get();
                if (!rowNode || rowNode->tag() != Node::kLine ||
                    !collect_unanchored_aligned_rows(
                        static_cast<const LineNode&>(*rowNode), rows)) {
                    return false;
                }
            }
            return true;
        }
    }
    rows.push_back(&line);
    return true;
}

}  // namespace

std::string latex_to_office_mathml_fragment(
        const std::string& latex, double pointSize) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return {};
    std::vector<const LineNode*> rows;
    if (!collect_unanchored_aligned_rows(*root, rows)) return {};
    if (rows.size() > 1) {
        std::string fragment;
        for (const LineNode* row : rows) {
            if (!fragment.empty()) fragment += "<br>";
            fragment += tree_to_mathml(*row, pointSize);
        }
        return fragment;
    }
    return tree_to_mathml(*root, pointSize);
}

}  // namespace eqnedit
