#include "latex_emitter.h"
#include "tex_parser.h"
#include <iostream>

int main() {
    for (const char* input : {"f'(x)", "x''", "x'''", "x'^{2}"}) {
        if (eqnedit::tree_to_latex(*eqnedit::parse_latex(input)) != input) {
            std::cerr << "Legacy prime spelling changed: " << input << '\n';
            return 3;
        }
    }
    for (const auto& sample : {
            std::pair<const char*, const char*>{"a^{2}{}'", "a^{2}{}^{\\prime }"},
            {"{a+b}'", "{a+b}^{\\prime }"},
            {"{a+b}''", "{a+b}^{\\prime \\prime }"},
            {"{a^{2}}'''", "{a^{2}}^{\\prime \\prime \\prime }"},
            {"{a+b}'_{i}", "{a+b}_{i}^{\\prime }"},
            {"{a+b}'^{2}", "{a+b}^{\\prime 2}"}}) {
        const auto first = eqnedit::tree_to_latex(*eqnedit::parse_latex(sample.first));
        const auto second = eqnedit::tree_to_latex(*eqnedit::parse_latex(first));
        if (first != sample.second || first != second) {
            std::cerr << sample.first << " -> " << first << " -> " << second << '\n';
            return 2;
        }
    }
    for (auto kind : {eqnedit::EM_PRIME, eqnedit::EM_DPRIME,
                      eqnedit::EM_TPRIME, eqnedit::EM_BPRIME}) {
        for (const char* input : {"a^{2}", "x_{i}^{n}", "\\frac{a}{b}",
                                  "a+b", "\\mathbf{x}", "x", "\\alpha"}) {
            auto base = eqnedit::parse_latex(input);
            auto mark = std::make_unique<eqnedit::EmbellNode>();
            mark->embellType = kind;
            mark->content = std::move(base->children);
            eqnedit::LineNode line;
            line.children.push_back(std::move(mark));
            const auto tex = eqnedit::tree_to_latex(line);
            const auto reopened = eqnedit::tree_to_latex(*eqnedit::parse_latex(tex));
            if (tex.find("^{\\prime") == std::string::npos ||
                tex.find('\'') != std::string::npos || tex != reopened) {
                std::cerr << input << " -> " << tex << " -> " << reopened << '\n';
                return 1;
            }
            std::cout << tex << '\n';
        }
    }
    std::cout << "PASS: 28 prime bases, 6 grouped imports, 4 legacy spellings\n";
}
