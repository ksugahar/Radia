/*
 * mtef_rtf.cpp -- RTF spelling of an equation
 *
 * Transcribed from what Word puts on the clipboard, construct by construct:
 *
 *   <m:f><m:fPr/>             {\mf{\mfPr{\mctrlPr}}
 *   <m:num>                   {\mnum ...}
 *   <m:rad><m:degHide/>       {\mrad{\mradPr{\mdegHide on}{\mctrlPr}}{\mdeg}
 *   <m:d><m:begChr m:val="["> {\md{\mdPr{\mbegChr [}{\mendChr ]}{\mctrlPr}}
 *   <m:nary><m:chr m:val="S"> {\mnary{\mnaryPr{\mchr S}{\mlimLoc undOvr}...}}
 *   <m:r><m:sty m:val="p"/>   {\mr\mscr0\msty0 text}
 *
 * The run style values are Word's: 0 plain, 1 bold, 2 italic, 3 bold italic.
 *
 * This is Word's clipboard route only: PowerPoint reads the same RTF as plain
 * text (measured -- it arrives as a text box, not an equation), which is what
 * the MathML output is for.
 */
#include "mtef_rtf.h"
#include "math_writer.h"
#include "mtef_parser.h"
#include "tex_parser.h"

#include <cstdio>
#include <string>

namespace mtef {
namespace {

constexpr int kTextFont = 0;
constexpr int kMathFont = 34;      /* Cambria Math, the number Word itself uses */

/* Decode UTF-8 so each code point can become an RTF \u escape; RTF is a byte
 * format with no notion of UTF-8. */
std::u32string decode_utf8(const std::string& s) {
    std::u32string out;
    for (size_t i = 0; i < s.size();) {
        unsigned char c = s[i];
        uint32_t cp;
        int n;
        if (c < 0x80)                { cp = c;         n = 1; }
        else if ((c & 0xE0) == 0xC0) { cp = c & 0x1F;  n = 2; }
        else if ((c & 0xF0) == 0xE0) { cp = c & 0x0F;  n = 3; }
        else if ((c & 0xF8) == 0xF0) { cp = c & 0x07;  n = 4; }
        else { ++i; continue; }
        if (i + n > s.size()) break;
        for (int k = 1; k < n; ++k) cp = (cp << 6) | (uint32_t(s[i + k]) & 0x3F);
        out.push_back(char32_t(cp));
        i += n;
    }
    return out;
}

std::string rtf_text(const std::string& utf8) {
    std::string o;
    for (char32_t cp : decode_utf8(utf8)) {
        if (cp == U'\\' || cp == U'{' || cp == U'}') {
            o += '\\';
            o += char(cp);
        } else if (cp < 0x80) {
            o += char(cp);
        } else if (cp < 0x10000) {
            char buf[32];
            /* \uN is signed 16-bit; Word writes values above 32767 negative. */
            int v = int(cp);
            if (v > 32767) v -= 65536;
            std::snprintf(buf, sizeof(buf), "\\u%d?", v);
            o += buf;
        } else {
            uint32_t v = uint32_t(cp) - 0x10000;
            char buf[64];
            std::snprintf(buf, sizeof(buf), "\\u%d?\\u%d?",
                          int(0xD800 + (v >> 10)) - 65536,
                          int(0xDC00 + (v & 0x3FF)) - 65536);
            o += buf;
        }
    }
    return o;
}

/* A control word is separated from following text by one space, which the
 * reader consumes; an immediately following brace needs none. */
std::string cw(const char* name, const std::string& inner) {
    std::string s = std::string("{\\m") + name;
    if (!inner.empty() && inner[0] != '{' && inner[0] != '\\') s += ' ';
    s += inner;
    s += '}';
    return s;
}
std::string val(const char* name, const std::string& value) {
    return std::string("{\\m") + name + ' ' + rtf_text(value) + '}';
}
std::string on(const char* name) {
    return std::string("{\\m") + name + " on}";
}
const char* kCtrl = "{\\mctrlPr}";

std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }

class RtfSyntax : public MathSyntax {
public:
    std::string run(const std::string& utf8, int style) const override {
        const int sty = style == 2 ? 3 : style == 1 ? 2 : 0;
        char head[64];
        std::snprintf(head, sizeof(head), "{\\mr\\mscr0\\msty%d ", sty);
        return std::string(head) + rtf_text(utf8) + '}';
    }

    /* Named by their parent, as in OMML. */
    std::string row(const std::string& inner) const override { return inner; }

    std::string fraction(const std::string& num, const std::string& den,
                         bool slashed) const override {
        std::string pr = slashed ? val("type", "skw") : std::string();
        return cw("f", cw("fPr", pr + kCtrl) + cw("num", num) + cw("den", den));
    }

    std::string radical(const std::string& body, const std::string& index,
                        bool has_index) const override {
        std::string pr = has_index ? std::string() : on("degHide");
        return cw("rad", cw("radPr", pr + kCtrl) +
                         cw("deg", has_index ? index : std::string()) +
                         cw("e", body));
    }

    std::string script(const std::string& base, const std::string& sub,
                       const std::string& sup,
                       bool has_sub, bool has_sup,
                       bool limits) const override {
        if (limits) {
            std::string out = base;
            if (has_sub)
                out = cw("limLow", cw("limLowPr", kCtrl) + cw("e", out) +
                                   cw("lim", sub));
            if (has_sup)
                out = cw("limUpp", cw("limUppPr", kCtrl) + cw("e", out) +
                                   cw("lim", sup));
            return out;
        }
        if (has_sub && has_sup)
            return cw("sSubSup", cw("sSubSupPr", kCtrl) + cw("e", base) +
                                 cw("sub", sub) + cw("sup", sup));
        if (has_sub)
            return cw("sSub", cw("sSubPr", kCtrl) + cw("e", base) +
                              cw("sub", sub));
        return cw("sSup", cw("sSupPr", kCtrl) + cw("e", base) + cw("sup", sup));
    }

    std::string fence(const std::string& body, const char* beg,
                      const char* end) const override {
        std::string pr;
        if (beg && std::string(beg) != "(") pr += val("begChr", beg);
        if (end && std::string(end) != ")") pr += val("endChr", end);
        return cw("d", cw("dPr", pr + kCtrl) + cw("e", body));
    }

    std::string nary(uint32_t chr, const std::string& lower,
                     const std::string& upper, const std::string& body,
                     bool stacked, bool has_lower, bool has_upper) const override {
        std::string pr = val("chr", utf8_of(chr)) +
                         val("limLoc", stacked ? "undOvr" : "subSup");
        if (!has_lower) pr += on("subHide");
        if (!has_upper) pr += on("supHide");
        return cw("nary", cw("naryPr", pr + kCtrl) + cw("sub", lower) +
                          cw("sup", upper) + cw("e", body));
    }

    std::string matrix(int rows, int cols,
                       const std::vector<std::string>& cells) const override {
        std::string mpr = cw("mcs", cw("mc", cw("mcPr",
                              val("count", std::to_string(cols)) +
                              val("mcJc", "center"))));
        std::string body = cw("mPr", mpr + kCtrl);
        for (int r = 0; r < rows; ++r) {
            std::string row_rtf;
            for (int c = 0; c < cols; ++c) {
                size_t i = size_t(r) * size_t(cols) + size_t(c);
                row_rtf += cw("e", i < cells.size() ? cells[i] : std::string());
            }
            body += cw("mr", row_rtf);
        }
        return cw("m", body);
    }

    std::string stack(const std::vector<std::string>& lines) const override {
        std::string body = cw("eqArrPr", kCtrl);
        for (const auto& l : lines) body += cw("e", l);
        return cw("eqArr", body);
    }

    std::string accent(uint32_t chr, const std::string& body) const override {
        std::string pr = chr ? val("chr", utf8_of(chr)) : std::string();
        return cw("acc", cw("accPr", pr + kCtrl) + cw("e", body));
    }

    std::string bar(const std::string& body, bool over) const override {
        return cw("bar", cw("barPr", val("pos", over ? "top" : "bot") + kCtrl) +
                         cw("e", body));
    }

    std::string document(const std::string& inner, bool display) const override {
        std::string body = display
            ? std::string("{\\*\\moMathPara {\\*\\moMath ") + inner + "}}"
            : std::string("{\\*\\moMath ") + inner + "}";
        char font[32];
        std::snprintf(font, sizeof(font), "\\f%d ", kMathFont);
        return std::string("{\\mmath") + font + body + '}';
    }
};

/* The smallest RTF file that carries one equation.  A bare fragment is not
 * something Word will paste. */
std::string rtf_document(const std::string& math, double size_pt) {
    char size[32];
    std::snprintf(size, sizeof(size), "\\fs%d", int(size_pt * 2.0 + 0.5));
    return std::string("{\\rtf1\\ansi\\ansicpg1252\\uc1\\deff") +
           std::to_string(kTextFont) +
           "{\\fonttbl{\\f" + std::to_string(kTextFont) +
           "\\froman\\fcharset0 Times New Roman;}{\\f" +
           std::to_string(kMathFont) + "\\froman\\fcharset0 Cambria Math;}}" +
           "{\\mmathPr\\mmathFont" + std::to_string(kMathFont) +
           "\\mdispDef1\\mwrapIndent1440}" +
           "\\pard\\plain" + size + " " + math + "}";
}

}  // namespace

std::string render_rtf_math(const LineNode& root, const RtfOptions& opt,
                            bool run_passes) {
    RtfSyntax syn;
    return write_math(root, syn, opt.display, run_passes);
}

std::string tex_to_rtf(const std::string& latex, const RtfOptions& opt) {
    std::unique_ptr<LineNode> root = parse_latex(latex);
    if (!root) return std::string();
    return rtf_document(render_rtf_math(*root, opt, /*run_passes=*/false),
                        opt.font_size_pt);
}

std::string mtef_to_rtf(const uint8_t* data, size_t len, const RtfOptions& opt) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return std::string();
    return rtf_document(render_rtf_math(*res.root, opt, /*run_passes=*/true),
                        opt.font_size_pt);
}

}  // namespace mtef
