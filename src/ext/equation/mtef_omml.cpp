/*
 * mtef_omml.cpp -- OMML, the spelling Office uses inside a .docx or .pptx
 *
 * Only the spelling lives here; the tree walk that decides what to emit is
 * shared with the RTF clipboard output (see math_writer.h).
 */
#include "mtef_omml.h"
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

class OmmlSyntax : public MathSyntax {
public:
    explicit OmmlSyntax(const OmmlOptions& opt) : opt_(opt) {}

    std::string group(const char* name, const std::string& inner) const override {
        return std::string("<m:") + name + '>' + inner + "</m:" + name + '>';
    }
    std::string prop(const char* name, const std::string& value) const override {
        return std::string("<m:") + name + " m:val=\"" + esc(value) + "\"/>";
    }
    std::string flag(const char* name) const override {
        return std::string("<m:") + name + " m:val=\"1\"/>";
    }
    /* OMML property groups need no placeholder; Word writes one in RTF only. */
    std::string ctrl() const override { return std::string(); }

    std::string run(const std::string& utf8, int style) const override {
        const char* sty = style == 2 ? "bi" : style == 1 ? nullptr : "p";
        std::string s = "<m:r>";
        if (sty && opt_.italic_variables)
            s += std::string("<m:rPr><m:sty m:val=\"") + sty + "\"/></m:rPr>";
        s += "<m:t>" + esc(utf8) + "</m:t></m:r>";
        return s;
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
};

}  // namespace

std::string render_omml(const LineNode& root, const OmmlOptions& opt,
                        bool run_passes) {
    OmmlSyntax syn(opt);
    return write_math(root, syn, opt.display, run_passes);
}

std::string mtef_to_omml(const uint8_t* data, size_t len, const OmmlOptions& opt) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return std::string();
    return render_omml(*res.root, opt, /*run_passes=*/true);
}

std::string tex_to_omml(const std::string& latex, const OmmlOptions& opt) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return render_omml(*root, opt, /*run_passes=*/false);
}

}  // namespace mtef
