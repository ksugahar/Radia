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
        return "<math xmlns=\"http://www.w3.org/1998/Math/MathML\" "
               "display=\"block\" mathsize=\"" + digits + "pt\">" +
               sequence(root.children) + "</math>";
    }

private:
    std::string sequence(const NodeList& nodes) {
        std::string content;
        for (const auto& node : nodes)
            if (node) content += emit_node(*node);
        return row(content);
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
        } else if (ch.typeface == TF_FUNCTION) {
            result = element("mi", value, "mathvariant=\"normal\"");
        } else if (ch.typeface == TF_VECTOR) {
            result = element("mi", value, "mathvariant=\"bold\"");
        } else if (ch.typeface == TF_SYMBOL ||
                   (ch.charCode < 128 && !std::isalnum(int(ch.charCode)))) {
            result = element("mo", value);
        } else {
            result = element("mi", value);
        }
        for (int embell : ch.embells) result = apply_embell(result, embell);
        return result;
    }

    std::string apply_embell(const std::string& base, int embell) {
        const char* mark = nullptr;
        bool under = false;
        switch (embell) {
        case EM_DOT: mark = "&#x02D9;"; break;
        case EM_DDOT: mark = "&#x00A8;"; break;
        case EM_TDOT: mark = "&#x20DB;"; break;
        case EM_PRIME: mark = "&#x2032;"; break;
        case EM_DPRIME: mark = "&#x2033;"; break;
        case EM_TPRIME: mark = "&#x2034;"; break;
        case EM_TILDE: mark = "~"; break;
        case EM_HAT: mark = "^"; break;
        case EM_NOT: mark = "&#x0338;"; break;
        case EM_RARROW: mark = "&#x2192;"; break;
        case EM_LARROW: mark = "&#x2190;"; break;
        case EM_BARROW: mark = "&#x2194;"; break;
        case EM_MBAR: mark = "_"; under = true; break;
        case EM_OBAR: mark = "&#x00AF;"; break;
        case EM_FROWN: mark = "&#x2322;"; break;
        case EM_SMILE: mark = "&#x2323;"; under = true; break;
        default: return base;
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

    std::string table(const MatrixNode& matrix) {
        std::string rows;
        const int rowCount = std::max(0, matrix.rows);
        const int columnCount = std::max(1, matrix.cols);
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
                ? "columnalign=\"right left\""
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
            return row(limited_operator(symbol, value.lower, value.hasLower,
                value.upper, value.hasUpper) + sequence(value.body));
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
            return row(limited_operator(symbol, value.lower, value.hasLower,
                value.upper, value.hasUpper) + sequence(value.body));
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

}  // namespace eqnedit
