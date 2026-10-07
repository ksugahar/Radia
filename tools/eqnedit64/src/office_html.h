#ifndef EQNEDIT_OFFICE_HTML_H
#define EQNEDIT_OFFICE_HTML_H
#include <string>
namespace eqnedit {
// Normal-copy Office payload. Throws on unsupported structures before the
// clipboard is opened. Registered MathML must not compete with this format.
std::string latex_to_office_html(const std::string &latex);
} // namespace eqnedit
#endif
