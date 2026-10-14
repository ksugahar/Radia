#include "equation_edit.h"
#include "latex_emitter.h"
#include "mathml_emitter.h"
#include "tex_parser.h"
#include "named_colors.h"
#include <iostream>
#include <stdexcept>
#include <string>

void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

int main() {
    using namespace eqnedit;
    for (const char* name : {"red", "blue", "green", "black", "cyan", "magenta",
            "yellow", "white", "gray", "orange", "purple", "brown", "lime",
            "olive", "pink", "teal", "violet", "darkgray", "lightgray",
            "Red", "Green", "Blue", "Apricot", "YellowGreen"}) {
        for (const std::string& tex : {
                "\\color{" + std::string(name) + "}{\\frac{x}{y}}+z",
                "{\\color{" + std::string(name) + "}x+y}z",
                "\\textcolor{" + std::string(name) + "}{\\int_a^b x\\,dx}"}) {
            auto root = parse_latex(tex);
            require(bool(root), "colour parse failed");
            const auto normalized = tree_to_latex(*root);
            auto reopened = parse_latex(normalized);
            require(bool(reopened), "colour reopen failed");
            require(tree_to_latex(*reopened) == normalized, "colour not fixed point");
            const std::string attribute = "mathcolor=\"" + named_color_hex(name) + "\"";
            require(tree_to_mathml(*root).find(attribute) != std::string::npos,
                    "registered MathML lost colour");
            require(latex_to_office_mathml_fragment(normalized).find(attribute) != std::string::npos,
                    "CF_HTML lost colour");
            Equation equation;
            require(equation.load_latex(tex), "editor colour load failed");
            equation.select_all();
            require(equation.selection_latex() == normalized, "selection lost colour");
            require(equation.insert_latex("\\textcolor{blue}{q}"), "edit colour failed");
            equation.undo();
            require(equation.latex() == normalized, "undo lost colour");
        }
    }
    for (const char* input : {"\\color{not-a-colour}{x}", "\\color[rgb]{1,0,0}{x}",
            "\\textcolor[HTML]{FF0000}{x}", "\\color{", "\\textcolor{red}",
            "\\sqrt\\color{unsupported}{x}"}) {
        std::string error;
        require(!parse_latex(input, nullptr, &error), "unsupported colour accepted");
        require(!error.empty(), "unsupported colour failed silently");
        Equation equation;
        equation.load_latex("x");
        require(!equation.replace_latex(input, false), "bad source accepted");
        require(equation.latex() == "x", "rejected source changed model");
        require(!equation.last_error().empty(), "editor error missing");
    }
    const auto nested = latex_to_mathml("\\textcolor{red}{a+\\textcolor{blue}{b}}+c");
    require(nested.find("mathcolor=\"#FF0000\"") != std::string::npos &&
            nested.find("mathcolor=\"#0000FF\"") != std::string::npos, "nested colour missing");
    const auto switchTree = parse_latex("\\color{red}{x}+y");
    require(switchTree->children.size() == 1 &&
            switchTree->children[0]->tag() == Node::kGroup, "switch colour not scoped");
    require(static_cast<GroupNode&>(*switchTree->children[0]).children.size() == 3,
            "switch colour failed to include following atoms");
    const auto bounded = parse_latex("\\textcolor{red}{x}+y");
    require(bounded->children.size() == 3, "textcolor leaked over following atoms");
    const auto rows = latex_to_office_mathml_fragment(
        "\\textcolor{red}{\\begin{aligned}x\\\\y\\end{aligned}}");
    require(rows.find("<br>") != std::string::npos, "coloured independent rows merged");
    require(rows.find("<math mathcolor=\"#FF0000\"") != std::string::npos,
            "independent rows lost scoped colour");
    Equation partial;
    partial.load_latex("\\textcolor{red}{xy}");
    partial.move_left(); // Enter the colour body at its end.
    partial.select_step_left();
    require(partial.selection_latex() == "\\textcolor{red}{y}",
            "copy of an inner selection lost inherited colour");
    require(latex_to_mathml("\\color{ red }{x}").find("#FF0000") != std::string::npos,
            "spaced CSS colour name rejected");
    require(named_color_hex(" Red ") == "#FF0000", "spaced CSS Red confused with dvips Red");
    std::cout << "PASS: named colours, scope, editor undo, fixed point, MathML and rejection\n";
}
