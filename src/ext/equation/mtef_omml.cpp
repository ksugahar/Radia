/*
 * mtef_omml.cpp -- OMML, the spelling Office uses inside a .docx or .pptx
 *
 * Only the spelling lives here; the tree walk that decides what to emit is
 * shared with the RTF and MathML outputs (see math_writer.h).
 */
#include "mtef_omml.h"
#include "math_writer.h"
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
    return std::string("<m:") + name + '>' + inner + "</m:" + name + '>';
}
std::string attr(const char* name, const std::string& value) {
    return std::string("<m:") + name + " m:val=\"" + esc(value) + "\"/>";
}
std::string flag(const char* name) {
    return std::string("<m:") + name + " m:val=\"1\"/>";
}

class OmmlSyntax : public MathSyntax {
public:
    explicit OmmlSyntax(const OmmlOptions& opt) : opt_(opt) {}

    std::string run(const std::string& utf8, int style) const override {
        const char* sty = style == 2 ? "bi" : style == 1 ? nullptr : "p";
        std::string s = "<m:r>";
        if (sty && opt_.italic_variables)
            s += std::string("<m:rPr><m:sty m:val=\"") + sty + "\"/></m:rPr>";
        s += "<m:t>" + esc(utf8) + "</m:t></m:r>";
        return s;
    }

    /* OMML slots are named by their parent, so a row is just its contents; the
     * parent wraps it. */
    std::string row(const std::string& inner) const override { return inner; }

    std::string fraction(const std::string& num, const std::string& den,
                         bool slashed) const override {
        std::string pr = slashed ? attr("type", "skw") : std::string();
        return elem("f", elem("fPr", pr) + elem("num", num) + elem("den", den));
    }

    std::string radical(const std::string& body, const std::string& index,
                        bool has_index) const override {
        std::string pr = has_index ? std::string() : flag("degHide");
        return elem("rad", elem("radPr", pr) +
                           elem("deg", has_index ? index : std::string()) +
                           elem("e", body));
    }

    std::string script(const std::string& base, const std::string& sub,
                       const std::string& sup,
                       bool has_sub, bool has_sup,
                       bool limits) const override {
        if (limits) {
            /* Word draws \lim_{x \to 0} as a lower limit, not a subscript. */
            std::string out = base;
            if (has_sub)
                out = elem("limLow", elem("limLowPr", "") + elem("e", out) +
                                     elem("lim", sub));
            if (has_sup)
                out = elem("limUpp", elem("limUppPr", "") + elem("e", out) +
                                     elem("lim", sup));
            return out;
        }
        if (has_sub && has_sup)
            return elem("sSubSup", elem("sSubSupPr", "") + elem("e", base) +
                                   elem("sub", sub) + elem("sup", sup));
        if (has_sub)
            return elem("sSub", elem("sSubPr", "") + elem("e", base) +
                                elem("sub", sub));
        return elem("sSup", elem("sSupPr", "") + elem("e", base) +
                            elem("sup", sup));
    }

    std::string fence(const std::string& body, const char* beg,
                      const char* end) const override {
        std::string pr;
        /* OMML defaults to parentheses, so only other shapes are spelled out. */
        if (beg && std::string(beg) != "(") pr += attr("begChr", beg);
        if (end && std::string(end) != ")") pr += attr("endChr", end);
        return elem("d", elem("dPr", pr) + elem("e", body));
    }

    std::string nary(uint32_t chr, const std::string& lower,
                     const std::string& upper, const std::string& body,
                     bool stacked, bool has_lower, bool has_upper) const override {
        std::string pr = attr("chr", utf8_of(chr)) +
                         attr("limLoc", stacked ? "undOvr" : "subSup");
        if (!has_lower) pr += flag("subHide");
        if (!has_upper) pr += flag("supHide");
        return elem("nary", elem("naryPr", pr) + elem("sub", lower) +
                            elem("sup", upper) + elem("e", body));
    }

    std::string matrix(int rows, int cols,
                       const std::vector<std::string>& cells) const override {
        std::string pr = elem("mcs", elem("mc", elem("mcPr",
                             attr("count", std::to_string(cols)) +
                             attr("mcJc", "center"))));
        std::string body = elem("mPr", pr);
        for (int r = 0; r < rows; ++r) {
            std::string row_xml;
            for (int c = 0; c < cols; ++c) {
                size_t i = size_t(r) * size_t(cols) + size_t(c);
                row_xml += elem("e", i < cells.size() ? cells[i] : std::string());
            }
            body += elem("mr", row_xml);
        }
        return elem("m", body);
    }

    std::string stack(const std::vector<std::string>& lines) const override {
        std::string body = elem("eqArrPr", "");
        for (const auto& l : lines) body += elem("e", l);
        return elem("eqArr", body);
    }

    std::string accent(uint32_t chr, const std::string& body) const override {
        std::string pr = chr ? attr("chr", utf8_of(chr)) : std::string();
        return elem("acc", elem("accPr", pr) + elem("e", body));
    }

    std::string bar(const std::string& body, bool over) const override {
        return elem("bar", elem("barPr", attr("pos", over ? "top" : "bot")) +
                           elem("e", body));
    }

    std::string document(const std::string& inner, bool display) const override {
        const char* ns = opt_.declare_namespace
            ? " xmlns:m=\"http://schemas.openxmlformats.org/officeDocument/2006/math\""
            : "";
        if (display)
            return std::string("<m:oMathPara") + ns + "><m:oMath>" + inner +
                   "</m:oMath></m:oMathPara>";
        return std::string("<m:oMath") + ns + '>' + inner + "</m:oMath>";
    }

private:
    const OmmlOptions& opt_;

    static std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }
};

}  // namespace

std::string render_omml(const LineNode& root, const OmmlOptions& opt,
                        bool run_passes) {
    OmmlSyntax syn(opt);
    return write_math(root, syn, opt.display, run_passes);
}

std::string tex_to_omml(const std::string& latex, const OmmlOptions& opt) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_omml(*root, opt, /*run_passes=*/false);
}

}  // namespace mtef
