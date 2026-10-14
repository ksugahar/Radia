#include "tex_document.h"

#include <algorithm>
#include <cctype>

namespace eqnedit {
namespace {

bool ascii_space(unsigned char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

std::string trim(std::string text) {
    if (text.size() >= 3 && text.compare(0, 3, "\xEF\xBB\xBF") == 0)
        text.erase(0, 3);
    while (!text.empty() && ascii_space(static_cast<unsigned char>(text.front())))
        text.erase(text.begin());
    while (!text.empty() && ascii_space(static_cast<unsigned char>(text.back())))
        text.pop_back();
    return text;
}

bool unwrap(std::string& text, const char* open, const char* close) {
    const std::string a(open), b(close);
    if (text.size() < a.size() + b.size() || text.compare(0, a.size(), a) != 0 ||
        text.compare(text.size() - b.size(), b.size(), b) != 0)
        return false;
    text = trim(text.substr(a.size(), text.size() - a.size() - b.size()));
    return true;
}

std::string unwrap_display_math(std::string text) {
    text = trim(std::move(text));
    bool changed = true;
    while (changed) {
        changed = unwrap(text, "$$", "$$") ||
                  unwrap(text, "\\[", "\\]") ||
                  unwrap(text, "\\(", "\\)") ||
                  unwrap(text, "$", "$");
    }
    return text;
}

}  // namespace

TexDocument parse_tex_document(const std::string& input) {
    TexDocument doc;
    const std::string plainBegin = "\\begin{equation}";
    const std::string plainEnd = "\\end{equation}";
    const std::string starBegin = "\\begin{equation*}";
    const std::string starEnd = "\\end{equation*}";

    size_t plain = input.find(plainBegin);
    size_t star = input.find(starBegin);
    size_t at = std::string::npos;
    if (plain != std::string::npos && star != std::string::npos)
        at = std::min(plain, star);
    else
        at = plain != std::string::npos ? plain : star;

    if (at == std::string::npos) {
        doc.body = unwrap_display_math(input);
        return doc;
    }

    const bool starred = at == star;
    const std::string& begin = starred ? starBegin : plainBegin;
    const std::string& end = starred ? starEnd : plainEnd;
    const size_t close = input.find(end, at + begin.size());
    if (close == std::string::npos) {
        doc.body = unwrap_display_math(input);
        return doc;
    }

    doc.prefix = input.substr(0, at);
    doc.suffix = input.substr(close + end.size());
    doc.body = unwrap_display_math(
        input.substr(at + begin.size(), close - at - begin.size()));
    doc.numbered = !starred;
    doc.hadEquationEnvironment = true;
    return doc;
}

std::string normalize_tex_paste(const std::string& text) {
    return parse_tex_document(text).body;
}

std::string compose_tex_document(const std::string& body, bool numbered,
                                 const std::string& prefix,
                                 const std::string& suffix) {
    const char* env = numbered ? "equation" : "equation*";
    std::string result = prefix + "\\begin{" + env + "}\n" + trim(body) +
                         "\n\\end{" + env + "}";
    if (suffix.empty()) result.push_back('\n');
    else result += suffix;
    return result;
}

}  // namespace eqnedit
