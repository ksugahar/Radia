#ifndef EQNEDIT_OFFICE_HTML_H
#define EQNEDIT_OFFICE_HTML_H
#include <string>
namespace eqnedit {
// Normal-copy Office payload. Throws on unsupported structures before the
// clipboard is opened. Registered MathML must not compete with this format.
std::string latex_to_office_html(const std::string &latex);
// Strict PowerPoint OOXML math export: m:* structure/text with a:rPr formatting.
// Shares the clipboard emitter's structure; this is not RichEdit HTML.
std::string latex_to_powerpoint_math(const std::string &latex);
} // namespace eqnedit
#endif
