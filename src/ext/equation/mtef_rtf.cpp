/*
 * mtef_rtf.cpp -- RTF spelling of an equation
 *
 * Transcribed from what Word puts on the clipboard, construct by construct:
 *
 *   <m:f><m:fPr/>            {\mf{\mfPr{\mctrlPr}}
 *   <m:num>                  {\mnum ...}
 *   <m:rad><m:degHide/>      {\mrad{\mradPr{\mdegHide on}{\mctrlPr}}{\mdeg}
 *   <m:d><m:begChr m:val="["> {\md{\mdPr{\mbegChr [}{\mendChr ]}{\mctrlPr}}
 *   <m:nary><m:chr m:val="Σ"> {\mnary{\mnaryPr{\mchr Σ}{\mlimLoc undOvr}...}}
 *   <m:r><m:sty m:val="p"/>   {\mr\mscr0\msty0 text}
 *
 * The run style values are Word's: 0 plain, 1 bold, 2 italic, 3 bold italic.
 */
#include "mtef_rtf.h"
#include "math_writer.h"
#include "mtef_parser.h"
#include "tex_parser.h"

#include <cstdio>
#include <string>

namespace mtef {
namespace {

/* Font numbers in the header written by rtf_document below. */
constexpr int kTextFont = 0;
constexpr int kMathFont = 34;      /* Cambria Math, the number Word itself uses */

/* Decode UTF-8 into code points so each can be written as an RTF \u escape.
 * RTF is a byte format with no notion of UTF-8. */
std::u32string decode_utf8(const std::string& s) {
    std::u32string out;
    for (size_t i = 0; i < s.size();) {
        unsigned char c = s[i];
        uint32_t cp;
        int n;
        if (c < 0x80)             { cp = c;          n = 1; }
        else if ((c & 0xE0) == 0xC0) { cp = c & 0x1F; n = 2; }
        else if ((c & 0xF0) == 0xE0) { cp = c & 0x0F; n = 3; }
        else if ((c & 0xF8) == 0xF0) { cp = c & 0x07; n = 4; }
        else { ++i; continue; }
        if (i + n > s.size()) break;
        for (int k = 1; k < n; ++k) cp = (cp << 6) | (uint32_t(s[i + k]) & 0x3F);
        out.push_back(char32_t(cp));
        i += n;
    }
    return out;
}

/* RTF text: braces and backslashes are syntax, and anything past ASCII needs
 * the \uN escape with an ASCII stand-in for readers that ignore it. */
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

class RtfSyntax : public MathSyntax {
public:
    std::string group(const char* name, const std::string& inner) const override {
        /* A control word is separated from following text by one space, which
         * the reader consumes; an immediately following brace needs none. */
        std::string s = std::string("{\\m") + name;
        if (!inner.empty() && inner[0] != '{' && inner[0] != '\\') s += ' ';
        s += inner;
        s += '}';
        return s;
    }
    std::string prop(const char* name, const std::string& value) const override {
        return std::string("{\\m") + name + ' ' + rtf_text(value) + '}';
    }
    std::string flag(const char* name) const override {
        return std::string("{\\m") + name + " on}";
    }
    std::string ctrl() const override { return "{\\mctrlPr}"; }

    std::string run(const std::string& utf8, int style) const override {
        /* \mscr0 = roman script; \msty 0 plain / 2 italic / 3 bold italic. */
        const int sty = style == 2 ? 3 : style == 1 ? 2 : 0;
        char head[64];
        std::snprintf(head, sizeof(head), "{\\mr\\mscr0\\msty%d ", sty);
        return std::string(head) + rtf_text(utf8) + '}';
    }

    std::string document(const std::string& inner, bool display) const override {
        /* Just the destination; rtf_document adds the file around it. */
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
