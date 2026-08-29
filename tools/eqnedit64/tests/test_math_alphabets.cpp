#include "latex_emitter.h"
#include "mathml_emitter.h"
#include "tex_parser.h"

#include <iostream>
#include <string>
#include <vector>

namespace {

std::string normalize(const std::string& tex) {
    return eqnedit::tree_to_latex(*eqnedit::parse_latex(tex));
}

}  // namespace

int main() {
    const std::vector<std::string> cases = {
        "\\mathrm{x}", "\\mathit{x}", "\\mathbf{x}",
        "\\mathsf{S}", "\\mathtt{T}", "\\mathcal{C}",
        "\\mathbb{R}", "\\mathfrak{F}", "\\bm{\\alpha}",
        "\\boldsymbol{\\beta}", "\\bm{\\wedge beta}",
        "\\bm{\\leq x}", "\\mathnormal{x}",
    };
    for (const std::string& input : cases) {
        const std::string once = normalize(input);
        const std::string twice = normalize(once);
        if (once != twice) {
            std::cerr << "not fixed: " << input << " -> " << once
                      << " -> " << twice << "\n";
            return 1;
        }
    }
    const std::string alias = normalize("\\boldsymbol{\\alpha}");
    if (alias.rfind("\\bm{", 0) != 0) {
        std::cerr << "boldsymbol did not canonicalize to bm: " << alias << "\n";
        return 2;
    }
    const std::string mathml = eqnedit::latex_to_mathml(
        "\\bm{\\alpha}\\mathsf{S}\\mathtt{T}\\mathcal{C}"
        "\\mathbb{R}\\mathfrak{F}");
    for (const char* variant : {"bold", "sans-serif", "monospace", "script",
                                "double-struck", "fraktur"}) {
        if (mathml.find(std::string("mathvariant=\"") + variant + "\"") ==
            std::string::npos) {
            std::cerr << "MathML lost variant " << variant << ": " << mathml
                      << "\n";
            return 3;
        }
    }
    std::cout << "PASS: " << cases.size()
              << " math alphabet normalization cases\n";
    return 0;
}
