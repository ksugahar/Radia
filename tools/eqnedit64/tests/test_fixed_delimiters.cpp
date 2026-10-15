#include "tex_parser.h"
#include "latex_emitter.h"
#include "mathml_emitter.h"
#include "office_html.h"
#include <iostream>
#include <stdexcept>
#include <string>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
int main() {
    using namespace eqnedit;
    for (const std::string base : {"\\big", "\\Big", "\\bigg", "\\Bigg"}) {
        for (const std::string suffix : {"", "l", "r", "m"}) {
            const std::string command = base + suffix;
            for (const std::string delimiter : {"(", ")", "[", "]", "|", "\\{", "\\}", "\\langle", "\\rangle"}) {
                const auto input = command + " " + delimiter;
                std::string error;
                auto tree = parse_latex(input, nullptr, &error, true);
                require(bool(tree), "strict sized delimiter parse failed");
                const auto once = tree_to_latex(*tree);
                require(once.find(command + " " + delimiter) != std::string::npos,
                        "exact delimiter command lost");
                auto reopened = parse_latex(once, nullptr, &error, true);
                require(bool(reopened) && tree_to_latex(*reopened) == once,
                        "sized delimiter not a fixed point");
                require(tree_to_mathml(*tree).find("minsize=") != std::string::npos,
                        "fixed size absent from MathML");
                require(latex_to_powerpoint_math(input).find("<m:d>") == std::string::npos,
                        "fixed token became grow-to-content delimiter");
            }
            for (const std::string argument : {"x", "\\alpha", ".", "", "{x}"}) {
                std::string error;
                require(!parse_latex(command + " " + argument, nullptr, &error, true) && !error.empty(),
                        "invalid sized delimiter argument accepted");
            }
            require(bool(parse_latex(command)), "half-typed delimiter broke interactive parsing");
            auto partial = parse_latex(command + " x");
            require(bool(partial) && tree_to_latex(*partial).find('x') != std::string::npos &&
                    tree_to_mathml(*partial).find("minsize=") == std::string::npos,
                    "interactive invalid argument lost its ordinary token");
        }
    }
    for (const std::string tex : {"\\mathbf{(x)}", "\\color{red}{(x)}"}) {
        const auto math = latex_to_powerpoint_math(tex);
        require(math.find("<m:d>") != std::string::npos,
                "styled parentheses lost their delimiter structure");
        require(math.find("<m:t") != std::string::npos &&
                math.find("<a:t") == std::string::npos,
                "PowerPoint math text is not canonical m:t");
    }
    std::cout << "PASS: fixed delimiter commands, TeX fixed points, styled fences and canonical math text\n";
}
