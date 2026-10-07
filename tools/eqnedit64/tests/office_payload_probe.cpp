#include "office_html.h"
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
int main(int argc, char **argv) {
  try {
    if (argc != 2)
      throw std::runtime_error("usage: office_payload_probe input.tex");
    std::ifstream input(argv[1], std::ios::binary);
    if (!input)
      throw std::runtime_error("cannot open input");
    const std::string tex((std::istreambuf_iterator<char>(input)),
                          std::istreambuf_iterator<char>());
    std::cout << eqnedit::latex_to_office_html(tex);
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
