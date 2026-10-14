#include "tex_document.h"

#include <iostream>
#include <string>

namespace {

int checks = 0;

bool expect(const std::string& actual, const std::string& wanted,
            const char* label) {
    ++checks;
    if (actual == wanted) return true;
    std::cerr << "FAIL " << label << "\nexpected: [" << wanted
              << "]\nactual:   [" << actual << "]\n";
    return false;
}

}  // namespace

int main() {
    bool ok = true;
    ok &= expect(eqnedit::normalize_tex_paste("  plain text  "), "plain text", "plain");
    ok &= expect(eqnedit::normalize_tex_paste("$x+y$"), "x+y", "inline dollars");
    ok &= expect(eqnedit::normalize_tex_paste("$$x+y$$"), "x+y", "display dollars");
    ok &= expect(eqnedit::normalize_tex_paste("\\(x+y\\)"), "x+y", "inline delimiters");
    ok &= expect(eqnedit::normalize_tex_paste("\\[x+y\\]"), "x+y", "display delimiters");

    const std::string aligned = "\\begin{aligned}a&=b\\\\c&=d\\end{aligned}";
    ok &= expect(eqnedit::normalize_tex_paste(aligned), aligned, "aligned retained");

    auto numbered = eqnedit::parse_tex_document(
        "before\n\\begin{equation}\n a+b \n\\end{equation}\nafter\n");
    ++checks;
    if (!numbered.hadEquationEnvironment || !numbered.numbered ||
        numbered.body != "a+b" || numbered.prefix != "before\n" ||
        numbered.suffix != "\nafter\n") {
        std::cerr << "FAIL numbered equation envelope\n";
        ok = false;
    }

    auto unnumbered = eqnedit::parse_tex_document(
        "\\begin{equation*}\nq\\end{equation*}");
    ++checks;
    if (!unnumbered.hadEquationEnvironment || unnumbered.numbered ||
        unnumbered.body != "q") {
        std::cerr << "FAIL starred equation envelope\n";
        ok = false;
    }

    ok &= expect(eqnedit::compose_tex_document("E=mc^{2}", true),
        "\\begin{equation}\nE=mc^{2}\n\\end{equation}\n", "numbered save");
    ok &= expect(eqnedit::compose_tex_document(aligned, true),
        "\\begin{equation}\n" + aligned +
        "\n\\end{equation}\n", "new multiline equation save");
    ok &= expect(eqnedit::compose_tex_document(aligned, false),
        "\\begin{equation*}\n" + aligned +
        "\n\\end{equation*}\n", "existing starred multiline save");

    if (!ok) return 1;
    std::cout << "PASS: " << checks << " TeX document/paste checks\n";
    return 0;
}
