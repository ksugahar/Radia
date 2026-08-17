/*
 * md_blocks.cpp -- split Markdown into blocks
 *
 * Line-oriented and single-pass.  Every line ends up in exactly one block and
 * every block keeps its original text, so concatenating the sources rebuilds
 * the file byte for byte -- the same guarantee MarkdownDoc gives for spans, and
 * for the same reason: an editor that cannot reproduce a file it did not write
 * is not one to trust with it.
 */
#include "md_blocks.h"

#include <cctype>

namespace mtef {
namespace {

/* A line, and where it ends in the source, so the block can keep it verbatim. */
struct Line {
    size_t begin = 0, end = 0;      /* end includes the newline, if any */
    std::string text;               /* without the newline */
};

std::vector<Line> split_lines(const std::string& s) {
    std::vector<Line> out;
    size_t i = 0;
    while (i < s.size()) {
        size_t nl = s.find('\n', i);
        Line l;
        l.begin = i;
        if (nl == std::string::npos) {
            l.end = s.size();
            l.text = s.substr(i);
            out.push_back(l);
            break;
        }
        l.end = nl + 1;
        l.text = s.substr(i, nl - i);
        if (!l.text.empty() && l.text.back() == '\r') l.text.pop_back();
        out.push_back(l);
        i = nl + 1;
    }
    return out;
}

bool blank(const std::string& s) {
    for (char c : s)
        if (!std::isspace((unsigned char)c)) return false;
    return true;
}

size_t indent_of(const std::string& s) {
    size_t n = 0;
    while (n < s.size() && (s[n] == ' ' || s[n] == '\t')) ++n;
    return n;
}

/* A fence is three or more backticks or tildes, optionally indented. */
size_t fence_run(const std::string& s, char& marker) {
    size_t i = indent_of(s);
    if (i >= s.size()) return 0;
    char c = s[i];
    if (c != '`' && c != '~') return 0;
    size_t n = 0;
    while (i + n < s.size() && s[i + n] == c) ++n;
    if (n < 3) return 0;
    marker = c;
    return n;
}

bool heading_of(const std::string& s, int& level, std::string& text) {
    size_t i = indent_of(s);
    if (i > 3 || i >= s.size() || s[i] != '#') return false;
    size_t n = 0;
    while (i + n < s.size() && s[i + n] == '#') ++n;
    if (n > 6) return false;
    size_t j = i + n;
    if (j >= s.size() || (s[j] != ' ' && s[j] != '\t')) return false;
    while (j < s.size() && (s[j] == ' ' || s[j] == '\t')) ++j;
    level = int(n);
    text = s.substr(j);
    return true;
}

bool bullet_of(const std::string& s, int& level, std::string& text) {
    size_t i = indent_of(s);
    if (i >= s.size()) return false;
    char c = s[i];
    if (c != '-' && c != '*' && c != '+') return false;
    if (i + 1 >= s.size() || (s[i + 1] != ' ' && s[i + 1] != '\t')) return false;
    size_t j = i + 1;
    while (j < s.size() && (s[j] == ' ' || s[j] == '\t')) ++j;
    level = int(i);
    text = s.substr(j);
    return true;
}

bool numbered_of(const std::string& s, int& level, std::string& text) {
    size_t i = indent_of(s);
    size_t j = i;
    while (j < s.size() && std::isdigit((unsigned char)s[j])) ++j;
    if (j == i || j >= s.size()) return false;
    if (s[j] != '.' && s[j] != ')') return false;
    ++j;
    if (j >= s.size() || (s[j] != ' ' && s[j] != '\t')) return false;
    while (j < s.size() && (s[j] == ' ' || s[j] == '\t')) ++j;
    level = int(i);
    text = s.substr(j);
    return true;
}

bool starts_something_else(const std::string& s) {
    int lvl;
    std::string txt;
    char marker;
    return blank(s) || heading_of(s, lvl, txt) || bullet_of(s, lvl, txt) ||
           numbered_of(s, lvl, txt) || fence_run(s, marker) != 0;
}

}  // namespace

std::vector<MdBlock> md_blocks(const std::string& markdown) {
    std::vector<MdBlock> out;
    std::vector<Line> lines = split_lines(markdown);
    const std::string& s = markdown;

    size_t i = 0;
    while (i < lines.size()) {
        MdBlock b;
        size_t first = i;

        /* ---- fenced code ------------------------------------------------ */
        char marker = 0;
        size_t open = fence_run(lines[i].text, marker);
        if (open) {
            b.kind = MdBlock::kCode;
            size_t after = indent_of(lines[i].text) + open;
            b.info = lines[i].text.substr(after);
            ++i;
            std::string body;
            while (i < lines.size()) {
                char m2 = 0;
                size_t close = fence_run(lines[i].text, m2);
                if (close >= open && m2 == marker) { ++i; break; }
                if (!body.empty()) body += "\n";
                body += lines[i].text;
                ++i;
            }
            b.text = body;
            b.source = s.substr(lines[first].begin, lines[i - 1].end - lines[first].begin);
            out.push_back(b);
            continue;
        }

        /* ---- blank run --------------------------------------------------- */
        if (blank(lines[i].text)) {
            b.kind = MdBlock::kBlank;
            while (i < lines.size() && blank(lines[i].text)) ++i;
            b.source = s.substr(lines[first].begin, lines[i - 1].end - lines[first].begin);
            out.push_back(b);
            continue;
        }

        /* ---- one-line kinds ---------------------------------------------- */
        int level = 0;
        std::string text;
        if (heading_of(lines[i].text, level, text)) {
            b.kind = MdBlock::kHeading;
        } else if (bullet_of(lines[i].text, level, text)) {
            b.kind = MdBlock::kBullet;
        } else if (numbered_of(lines[i].text, level, text)) {
            b.kind = MdBlock::kNumbered;
        } else {
            /* ---- paragraph: runs until something else starts ------------- */
            b.kind = MdBlock::kParagraph;
            std::string body;
            while (i < lines.size() && !starts_something_else(lines[i].text)) {
                if (!body.empty()) body += "\n";
                body += lines[i].text;
                ++i;
            }
            b.text = body;
            b.source = s.substr(lines[first].begin, lines[i - 1].end - lines[first].begin);
            out.push_back(b);
            continue;
        }

        b.level = level;
        b.text = text;
        b.source = s.substr(lines[i].begin, lines[i].end - lines[i].begin);
        ++i;
        out.push_back(b);
    }
    return out;
}

}  // namespace mtef
