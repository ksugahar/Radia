#include "equation_render.h"
#include <iostream>
#include <stdexcept>
#include <string>

template<class T>
bool rejects(const char* message) {
    eqnedit::LineNode root;
    root.children.push_back(std::make_unique<T>());
    try {
        eqnedit::render_svg(root);
    } catch (const std::invalid_argument& error) {
        return std::string(error.what()) == message;
    }
    return false;
}

int main() {
    if (!rejects<eqnedit::FontNode>("unsupported legacy font node") ||
        !rejects<eqnedit::RMNode>("unsupported legacy RM node")) {
        std::cerr << "legacy node silently rendered or raised the wrong error\n";
        return 1;
    }
    // Size is a state transition, not unsupported content.
    eqnedit::LineNode root;
    root.children.push_back(std::make_unique<eqnedit::SizeNode>());
    try {
        eqnedit::render_svg(root);
    } catch (const std::exception& error) {
        std::cerr << "size state failed: " << error.what() << '\n';
        return 2;
    }
    std::cout << "PASS: legacy nodes reject; size state is accepted\n";
}
