/*
 * mtef_mathml.cpp -- MathML spelling of an equation
 *
 * Where OMML and RTF share element names, MathML does not, which is why the
 * shared writer names meanings rather than tags: a root takes its arguments the
 * other way round, delimiters are ordinary operators inside a row, and a run is
 * an identifier, a number or an operator depending on what it contains.
 */
#include "mtef_mathml.h"
#include "math_writer.h"
#include "mtef_parser.h"
#include "tex_parser.h"

#include <string>

namespace mtef {
namespace {

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

std::string elem(const char* name, const std::string& inner) {
    return std::string("<") + name + '>' + inner + "</" + name + '>';
}

std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }

bool all_digits(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s)
        if (c < '0' || c > '9') return false;
    return true;
}

/* An empty MathML slot still has to be an element, or the row arity is wrong
 * and the renderer drops the construct. */
std::string row_of(const std::string& inner) {
    return elem("mrow", inner);
}

class MathMLSyntax : public MathSyntax {
public:
    explicit MathMLSyntax(const MathMLOptions& opt) : opt_(opt) {}

    /* MathML splits a run three ways, and the split carries the spacing: an
     * <mo> gets operator spacing, an <mi> is italicised as a variable, an <mn>
     * is a number. */
    std::string run(const std::string& utf8, int style) const override {
        if (all_digits(utf8)) return elem("mn", esc(utf8));
        if (style == 0) return elem("mo", esc(utf8));
        const char* tag = "mi";
        if (style == 2)
            return std::string("<mi mathvariant=\"bold-italic\">") + esc(utf8) +
                   "</mi>";
        return elem(tag, esc(utf8));
    }

    std::string row(const std::string& inner) const override {
        return row_of(inner);
    }

    std::string fraction(const std::string& num, const std::string& den,
                         bool slashed) const override {
        if (slashed)
            return std::string("<mfrac bevelled=\"true\">") + num + den +
                   "</mfrac>";
        return elem("mfrac", num + den);
    }

    /* <mroot> takes the body first and the index second -- the reverse of the
     * order the tree stores them in. */
    std::string radical(const std::string& body, const std::string& index,
                        bool has_index) const override {
        if (has_index) return elem("mroot", body + index);
        return elem("msqrt", body);
    }

    std::string script(const std::string& base, const std::string& sub,
                       const std::string& sup,
                       bool has_sub, bool has_sup,
                       bool limits) const override {
        if (limits) {
            if (has_sub && has_sup) return elem("munderover", base + sub + sup);
            if (has_sub)            return elem("munder", base + sub);
            return elem("mover", base + sup);
        }
        if (has_sub && has_sup) return elem("msubsup", base + sub + sup);
        if (has_sub)            return elem("msub", base + sub);
        return elem("msup", base + sup);
    }

    /* MathML has no fence element in current use: the delimiters are operators
     * inside a row, which is what Word writes. */
    std::string fence(const std::string& body, const char* beg,
                      const char* end) const override {
        std::string s;
        if (beg && *beg) s += elem("mo", esc(beg));
        s += body;
        if (end && *end) s += elem("mo", esc(end));
        return row_of(s);
    }

    std::string nary(uint32_t chr, const std::string& lower,
                     const std::string& upper, const std::string& body,
                     bool stacked, bool has_lower, bool has_upper) const override {
        std::string op = std::string("<mo stretchy=\"false\">") +
                         esc(utf8_of(chr)) + "</mo>";
        std::string head;
        if (has_lower && has_upper)
            head = elem(stacked ? "munderover" : "msubsup", op + lower + upper);
        else if (has_lower)
            head = elem(stacked ? "munder" : "msub", op + lower);
        else if (has_upper)
            head = elem(stacked ? "mover" : "msup", op + upper);
        else
            head = op;
        return row_of(head + body);
    }

    std::string matrix(int rows, int cols,
                       const std::vector<std::string>& cells) const override {
        std::string body;
        for (int r = 0; r < rows; ++r) {
            std::string row_xml;
            for (int c = 0; c < cols; ++c) {
                size_t i = size_t(r) * size_t(cols) + size_t(c);
                row_xml += elem("mtd", i < cells.size() ? cells[i] : std::string());
            }
            body += elem("mtr", row_xml);
        }
        return elem("mtable", body);
    }

    std::string stack(const std::vector<std::string>& lines) const override {
        std::string body;
        for (const auto& l : lines) body += elem("mtr", elem("mtd", l));
        return elem("mtable", body);
    }

    std::string accent(uint32_t chr, const std::string& body) const override {
        uint32_t g = chr ? chr : 0x0302;      /* default circumflex */
        return std::string("<mover accent=\"true\">") + body +
               "<mo>" + esc(utf8_of(g)) + "</mo></mover>";
    }

    std::string bar(const std::string& body, bool over) const override {
        const char* tag = over ? "mover" : "munder";
        uint32_t g = over ? 0x00AF : 0x005F;  /* macron / low line */
        return std::string("<") + tag + " accent=\"true\">" + body +
               "<mo stretchy=\"true\">" + esc(utf8_of(g)) + "</mo></" + tag + ">";
    }

    std::string document(const std::string& inner, bool display) const override {
        std::string s = "<math";
        if (opt_.declare_namespace)
            s += " xmlns=\"http://www.w3.org/1998/Math/MathML\"";
        s += std::string(" display=\"") + (display ? "block" : "inline") + "\">";
        s += inner;
        s += "</math>";
        return s;
    }

private:
    const MathMLOptions& opt_;
};

}  // namespace

std::string render_mathml(const LineNode& root, const MathMLOptions& opt,
                          bool run_passes) {
    MathMLSyntax syn(opt);
    return write_math(root, syn, opt.display, run_passes);
}

std::string tex_to_mathml(const std::string& latex, const MathMLOptions& opt) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_mathml(*root, opt, /*run_passes=*/false);
}

std::string mtef_to_mathml(const uint8_t* data, size_t len,
                           const MathMLOptions& opt) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return std::string();
    return render_mathml(*res.root, opt, /*run_passes=*/true);
}

}  // namespace mtef
