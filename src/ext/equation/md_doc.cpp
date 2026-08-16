/*
 * md_doc.cpp -- scan Markdown into text and math spans
 *
 * One left-to-right pass.  At each position the scanner asks, in order:
 * am I at a fenced code block, an inline code span, an escape, or a math
 * delimiter?  Anything else accumulates into the current text run.
 */
#include "md_doc.h"

#include <cctype>

namespace mtef {
namespace {

bool at_line_start(const std::string& s, size_t i) {
    return i == 0 || s[i - 1] == '\n';
}

/* A fence is three or more backticks or tildes at the start of a line. */
size_t fence_run(const std::string& s, size_t i, char c) {
    size_t n = 0;
    while (i + n < s.size() && s[i + n] == c) ++n;
    return n >= 3 ? n : 0;
}

size_t line_end(const std::string& s, size_t i) {
    size_t nl = s.find('\n', i);
    return nl == std::string::npos ? s.size() : nl + 1;
}

bool is_space(char c) { return c == ' ' || c == '\t' || c == '\n' || c == '\r'; }

/* Pandoc's inline rule: no space just inside the delimiters, and the closing
 * one is not followed by a digit -- so "$5 and $6" is prose, not math. */
bool inline_math_at(const std::string& s, size_t i, size_t& body_begin,
                    size_t& body_end, size_t& next) {
    if (s[i] != '$' || i + 1 >= s.size()) return false;
    if (is_space(s[i + 1])) return false;
    size_t j = i + 1;
    while (j < s.size()) {
        if (s[j] == '\\' && j + 1 < s.size()) { j += 2; continue; }
        if (s[j] == '\n' && j + 1 < s.size() && s[j + 1] == '\n') return false;
        /* A code span cannot live inside inline math.  Without this, prose
         * like "$5 and $6, where `$HOME` is set" has its second dollar open a
         * span that closes on the one inside the code span. */
        if (s[j] == '`') return false;
        if (s[j] == '$') break;
        ++j;
    }
    if (j >= s.size() || s[j] != '$') return false;
    if (j == i + 1) return false;                       /* "$$" is not inline */
    if (is_space(s[j - 1])) return false;
    if (j + 1 < s.size() && std::isdigit((unsigned char)s[j + 1])) return false;
    body_begin = i + 1;
    body_end = j;
    next = j + 1;
    return true;
}

}  // namespace

void MarkdownDoc::load(const std::string& s) {
    segs_.clear();
    std::string run;
    size_t i = 0;

    auto flush = [&]() {
        if (run.empty()) return;
        MdSegment seg;
        seg.kind = MdSegment::kText;
        seg.body = run;
        segs_.push_back(seg);
        run.clear();
    };

    while (i < s.size()) {
        /* ---- fenced code block ---------------------------------------- */
        if (at_line_start(s, i) && (s[i] == '`' || s[i] == '~')) {
            size_t n = fence_run(s, i, s[i]);
            if (n) {
                char c = s[i];
                size_t start = i;
                size_t j = line_end(s, i);
                while (j < s.size()) {
                    if (at_line_start(s, j)) {
                        size_t k = j;
                        while (k < s.size() && (s[k] == ' ')) ++k;
                        if (fence_run(s, k, c) >= n) { j = line_end(s, k); break; }
                    }
                    j = line_end(s, j);
                }
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kCodeBlock;
                seg.body = s.substr(start, j - start);
                segs_.push_back(seg);
                i = j;
                continue;
            }
        }

        /* ---- inline code span ------------------------------------------ */
        if (s[i] == '`') {
            size_t n = 0;
            while (i + n < s.size() && s[i + n] == '`') ++n;
            size_t j = i + n;
            size_t close = std::string::npos;
            while (j < s.size()) {
                if (s[j] == '`') {
                    size_t m = 0;
                    while (j + m < s.size() && s[j + m] == '`') ++m;
                    if (m == n) { close = j; break; }
                    j += m;
                    continue;
                }
                ++j;
            }
            if (close != std::string::npos) {
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kCodeSpan;
                seg.body = s.substr(i, close + n - i);
                segs_.push_back(seg);
                i = close + n;
                continue;
            }
        }

        /* ---- escapes: \$ is a dollar sign, not a delimiter -------------- */
        if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == '$') {
            run += s[i]; run += s[i + 1];
            i += 2;
            continue;
        }

        /* ---- display math ---------------------------------------------- */
        if (s[i] == '$' && i + 1 < s.size() && s[i + 1] == '$') {
            size_t j = s.find("$$", i + 2);
            if (j != std::string::npos) {
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kDisplayMath;
                seg.open = "$$";
                seg.body = s.substr(i + 2, j - (i + 2));
                seg.close = "$$";
                segs_.push_back(seg);
                i = j + 2;
                continue;
            }
        }
        if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == '[') {
            size_t j = s.find("\\]", i + 2);
            if (j != std::string::npos) {
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kDisplayMath;
                seg.open = "\\[";
                seg.body = s.substr(i + 2, j - (i + 2));
                seg.close = "\\]";
                segs_.push_back(seg);
                i = j + 2;
                continue;
            }
        }

        /* ---- inline math ------------------------------------------------ */
        if (s[i] == '\\' && i + 1 < s.size() && s[i + 1] == '(') {
            size_t j = s.find("\\)", i + 2);
            if (j != std::string::npos) {
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kInlineMath;
                seg.open = "\\(";
                seg.body = s.substr(i + 2, j - (i + 2));
                seg.close = "\\)";
                segs_.push_back(seg);
                i = j + 2;
                continue;
            }
        }
        {
            size_t b, e, next;
            if (inline_math_at(s, i, b, e, next)) {
                flush();
                MdSegment seg;
                seg.kind = MdSegment::kInlineMath;
                seg.open = "$";
                seg.body = s.substr(b, e - b);
                seg.close = "$";
                segs_.push_back(seg);
                i = next;
                continue;
            }
        }

        run += s[i];
        ++i;
    }
    flush();
}

std::string MarkdownDoc::text() const {
    std::string out;
    for (const MdSegment& s : segs_) out += s.source();
    return out;
}

int MarkdownDoc::math_count() const {
    int n = 0;
    for (const MdSegment& s : segs_) if (s.is_math()) ++n;
    return n;
}

int MarkdownDoc::math_segment_index(int i) const {
    int n = 0;
    for (size_t k = 0; k < segs_.size(); ++k) {
        if (!segs_[k].is_math()) continue;
        if (n == i) return int(k);
        ++n;
    }
    return -1;
}

std::string MarkdownDoc::math_latex(int i) const {
    int k = math_segment_index(i);
    return k < 0 ? std::string() : segs_[size_t(k)].body;
}

bool MarkdownDoc::math_is_display(int i) const {
    int k = math_segment_index(i);
    return k >= 0 && segs_[size_t(k)].kind == MdSegment::kDisplayMath;
}

bool MarkdownDoc::set_math_latex(int i, const std::string& latex) {
    int k = math_segment_index(i);
    if (k < 0) return false;
    segs_[size_t(k)].body = latex;
    return true;
}

}  // namespace mtef
