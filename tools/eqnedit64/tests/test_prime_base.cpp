#include "latex_emitter.h"
#include "tex_parser.h"
#include <iostream>

int main() {
    for (auto kind : {eqnedit::EM_PRIME, eqnedit::EM_DPRIME,
                      eqnedit::EM_TPRIME, eqnedit::EM_BPRIME}) {
        for (const char* input : {"a^{2}", "x_{i}^{n}", "\\frac{a}{b}",
                                  "a+b", "\\mathbf{x}"}) {
            auto base = eqnedit::parse_latex(input);
            auto mark = std::make_unique<eqnedit::EmbellNode>();
            mark->embellType = kind;
            mark->content = std::move(base->children);
            eqnedit::LineNode line;
            line.children.push_back(std::move(mark));
            const auto tex = eqnedit::tree_to_latex(line);
            const auto reopened = eqnedit::tree_to_latex(*eqnedit::parse_latex(tex));
            if (tex.find("}^{\\prime") == std::string::npos || tex != reopened) {
                std::cerr << input << " -> " << tex << " -> " << reopened << '\n';
                return 1;
            }
            std::cout << tex << '\n';
        }
    }
    std::cout << "PASS: 20 compound prime bases survive first reparse\n";
}
