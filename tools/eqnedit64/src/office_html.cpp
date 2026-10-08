#include "office_html.h"
#include "mathml_emitter.h"
#include "named_colors.h"
#include "office_alphabets.h"
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <map>
#include <memory>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>
namespace eqnedit {
namespace {
std::string utf8(uint32_t c) {
  std::string s;
  if (c < 0x80)
    s += char(c);
  else if (c < 0x800) {
    s += char(0xC0 | (c >> 6));
    s += char(0x80 | (c & 63));
  } else if (c < 0x10000) {
    s += char(0xE0 | (c >> 12));
    s += char(0x80 | ((c >> 6) & 63));
    s += char(0x80 | (c & 63));
  } else {
    s += char(0xF0 | (c >> 18));
    s += char(0x80 | ((c >> 12) & 63));
    s += char(0x80 | ((c >> 6) & 63));
    s += char(0x80 | (c & 63));
  }
  return s;
}
std::string escape(const std::string &s) {
  std::string out;
  for (char c : s) {
    if (c == '&')
      out += "&amp;";
    else if (c == '<')
      out += "&lt;";
    else if (c == '>')
      out += "&gt;";
    else if (c == '\"')
      out += "&quot;";
    else
      out += c;
  }
  return out;
}
std::string decode(const std::string &s) {
  std::string out;
  for (size_t i = 0; i < s.size(); i++) {
    if (s[i] != '&') {
      out += s[i];
      continue;
    }
    size_t end = s.find(';', i);
    if (end == std::string::npos)
      throw std::runtime_error("bad XML entity");
    const auto e = s.substr(i + 1, end - i - 1);
    if (e == "amp")
      out += '&';
    else if (e == "lt")
      out += '<';
    else if (e == "gt")
      out += '>';
    else if (e == "quot")
      out += '\"';
    else if (e == "apos")
      out += '\'';
    else if (e.size() > 1 && e[0] == '#')
      out += utf8(uint32_t(std::stoul(e.substr(e[1] == 'x' ? 2 : 1), nullptr,
                                      e[1] == 'x' ? 16 : 10)));
    else
      throw std::runtime_error("unknown XML entity");
    i = end;
  }
  return out;
}
struct Xml {
  std::string name, text;
  std::map<std::string, std::string> attrs;
  std::vector<std::unique_ptr<Xml>> children;
  Xml *parent = nullptr;
  std::string attr(const std::string &k) const {
    auto it = attrs.find(k);
    return it == attrs.end() ? "" : it->second;
  }
};
// Deliberately restricted to XML emitted by our MathML emitter: no
// declarations, DTD, external entities or mixed-content rewriting. Not an
// external XML reader.
class Parser {
  const std::string &s;
  size_t at = 0;
  void ws() {
    while (at < s.size() && std::isspace(static_cast<unsigned char>(s[at])))
      at++;
  }
  void expect(char c) {
    if (at >= s.size() || s[at++] != c)
      throw std::runtime_error("invalid internal MathML");
  }
  std::string name() {
    size_t start = at;
    while (at < s.size() && (std::isalnum(static_cast<unsigned char>(s[at])) ||
                             s[at] == '-' || s[at] == ':'))
      at++;
    if (start == at)
      throw std::runtime_error("missing XML name");
    return s.substr(start, at - start);
  }
  std::unique_ptr<Xml> node(Xml *parent, int depth) {
    if (depth > 256)
      throw std::runtime_error("Office tree too deep");
    expect('<');
    auto n = std::make_unique<Xml>();
    n->parent = parent;
    n->name = name();
    ws();
    while (at < s.size() && s[at] != '>' && s[at] != '/') {
      auto k = name();
      ws();
      expect('=');
      ws();
      expect('\"');
      size_t end = s.find('\"', at);
      if (end == std::string::npos)
        throw std::runtime_error("invalid XML attribute");
      n->attrs[k] = decode(s.substr(at, end - at));
      at = end + 1;
      ws();
    }
    if (at < s.size() && s[at] == '/') {
      at++;
      expect('>');
      return n;
    }
    expect('>');
    while (at < s.size()) {
      if (s.compare(at, 2, "</") == 0) {
        at += 2;
        if (name() != n->name)
          throw std::runtime_error("XML tag mismatch");
        expect('>');
        return n;
      }
      if (s[at] == '<')
        n->children.push_back(node(n.get(), depth + 1));
      else {
        size_t end = s.find('<', at);
        if (end == std::string::npos)
          throw std::runtime_error("unterminated XML");
        n->text += decode(s.substr(at, end - at));
        at = end;
      }
    }
    throw std::runtime_error("unterminated XML");
  }

public:
  explicit Parser(const std::string &value) : s(value) {}
  std::unique_ptr<Xml> parse() {
    auto root = node(nullptr, 0);
    ws();
    if (at != s.size())
      throw std::runtime_error("XML trailing data");
    return root;
  }
};
std::string inherited(const Xml &n, const std::string &key) {
  for (auto p = &n; p; p = p->parent) {
    auto it = p->attrs.find(key);
    if (it != p->attrs.end())
      return it->second;
  }
  return {};
}
std::string content(const Xml &n) {
  std::string s = n.text;
  for (auto &c : n.children)
    s += content(*c);
  return s;
}
const Xml &unwrapped(const Xml &n) {
  if ((n.name == "mrow" || n.name == "mstyle") && n.children.size() == 1)
    return unwrapped(*n.children[0]);
  return n;
}
std::string alphabet(const std::string &text, const std::string &variant) {
  auto found = office_alphabets().find(variant);
  if (found == office_alphabets().end())
    return text;
  const auto &a = found->second;
  std::string out;
  for (size_t i = 0; i < text.size();) {
    uint32_t cp = static_cast<unsigned char>(text[i++]);
    int extra = 0;
    if (cp >= 0xF0) {
      cp &= 7;
      extra = 3;
    } else if (cp >= 0xE0) {
      cp &= 15;
      extra = 2;
    } else if (cp >= 0xC0) {
      cp &= 31;
      extra = 1;
    }
    while (extra--) {
      if (i >= text.size())
        throw std::runtime_error("invalid UTF8");
      cp = (cp << 6) | (static_cast<unsigned char>(text[i++]) & 63);
    }
    auto hole = a.holes.find(cp);
    if (hole != a.holes.end())
      cp = hole->second;
    else if (cp >= 'A' && cp <= 'Z')
      cp = a.upper + cp - 'A';
    else if (cp >= 'a' && cp <= 'z')
      cp = a.lower + cp - 'a';
    else if (cp >= '0' && cp <= '9' && a.digit)
      cp = a.digit + cp - '0';
    else if (a.greekUpper) {
      if (cp >= 0x391 && cp <= 0x3A9 && cp != 0x3A2)
        cp = a.greekUpper + cp - 0x391;
      else if (cp >= 0x3B1 && cp <= 0x3C9)
        cp = a.greekLower + cp - 0x3B1;
      else if (cp == 0x3F4)
        cp = a.greekUpper + 17;
      else if (cp == 0x2207)
        cp = a.greekUpper + 25;
      else {
        const uint32_t symbols[] = {0x2202, 0x3F5, 0x3D1, 0x3F0,
                                    0x3D5,  0x3F1, 0x3D6};
        for (uint32_t j = 0; j < 7; j++)
          if (cp == symbols[j]) {
            cp = a.greekLower + 25 + j;
            break;
          }
      }
    }
    out += utf8(cp);
  }
  return out;
}
const std::string style = "font-size:18.0pt;font-family:'Cambria "
                          "Math';mso-ascii-font-family:'Cambria "
                          "Math';mso-font-kerning:12.0pt;color:black";
std::string control() {
  return "<span style=\"" + style +
         ";font-style:normal\"><m:ctrlPr></m:ctrlPr></span>";
}
std::string emit(const Xml &);
std::string run(const Xml &n) {
  if (content(n).empty())
    return {};
  // MathJax marks named-function application with nonprinting U+2061;
  // native FunctionNode has no such token. Neither route emits it in OMML.
  if (n.name == "mo" && n.text == "\xE2\x81\xA1")
    return {};
  const auto variant = inherited(n, "mathvariant");
  const bool mapped = office_alphabets().count(variant) != 0;
  if (!variant.empty() && variant != "normal" && variant != "italic" &&
      variant != "bold" && variant != "bold-italic" && !mapped)
    throw std::runtime_error("unsupported Office mathvariant");
  bool normal = mapped || variant == "normal" || n.name != "mi";
  std::string colour = inherited(n, "mathcolor");
  if (!colour.empty()) {
    auto canonical = named_color_hex(colour);
    if (canonical.empty() && colour.size() == 7 && colour[0] == '#')
      canonical = colour;
    colour = canonical;
    if (colour.empty())
      throw std::runtime_error("unsupported Office colour");
  }
  return "<m:r><span style=\"" + style +
         (colour.empty() ? "" : ";color:" + escape(colour)) +
         ";font-style:" + (normal ? "normal" : "italic") +
         (!mapped && variant.find("bold") != std::string::npos
              ? ";font-weight:bold"
              : "") +
         "\">" + escape(alphabet(content(n), variant)) + "</span></m:r>";
}
bool scripted(const Xml &n) {
  return n.name == "msub" || n.name == "msup" || n.name == "msubsup" ||
         n.name == "munder" || n.name == "mover" || n.name == "munderover";
}
bool large(const Xml &n) {
  const auto &x = unwrapped(n);
  return x.name == "mo" &&
         std::string("∑∏∐⋃⋂∫∬∭∮∯∰").find(x.text) != std::string::npos &&
         !x.text.empty();
}
std::string nary(const Xml &n, const std::string &body) {
  const auto &base = unwrapped(*n.children.at(0));
  const bool lo = n.name == "msub" || n.name == "msubsup" ||
                  n.name == "munder" || n.name == "munderover";
  const bool hi = n.name == "msup" || n.name == "msubsup" ||
                  n.name == "mover" || n.name == "munderover";
  const bool under =
      n.name == "munder" || n.name == "mover" || n.name == "munderover";
  return "<m:nary><m:naryPr><m:chr m:val=\"" + escape(base.text) +
         "\"/><m:limLoc m:val=\"" + (under ? "undOvr" : "subSup") +
         "\"/><m:grow m:val=\"on\"/>" +
         (!lo ? "<m:subHide m:val=\"1\"/>" : "") +
         (!hi ? "<m:supHide m:val=\"1\"/>" : "") + control() +
         "</m:naryPr><m:sub>" + (lo ? emit(*n.children.at(1)) : "") +
         "</m:sub><m:sup>" + (hi ? emit(*n.children.at(lo ? 2 : 1)) : "") +
         "</m:sup><m:e>" + body + "</m:e></m:nary>";
}
// A fence is recognized by its visible delimiter, on both producers, rather
// than their different fence/stretchy/largeop attributes. Atomic mrows keep
// their own boundary; a singleton wrapper may still contain a fence token.
std::string fenceClose(const std::string &token) {
  static const std::map<std::string, std::string> pairs = {
      {"(", ")"}, {"[", "]"}, {"{", "}"}, {"⟨", "⟩"}, {"⌈", "⌉"},
      {"⌊", "⌋"}, {"⟦", "⟧"}, {"|", "|"}, {"‖", "‖"}};
  auto it = pairs.find(token);
  return it == pairs.end() ? "" : it->second;
}
bool closingFence(const std::string &token) {
  return token == ")" || token == "]" || token == "}" || token == "⟩" ||
         token == "⌉" || token == "⌋" || token == "⟧";
}
void advanceFence(const Xml &node, std::vector<std::string> &stack) {
  const auto &token = unwrapped(node);
  if (token.name != "mo")
    return;
  if (!stack.empty() && stack.back() == token.text) {
    stack.pop_back();
  } else {
    const auto close = fenceClose(token.text);
    if (!close.empty())
      stack.push_back(close);
  }
}
bool operandBoundary(const Xml &node, const std::vector<std::string> &inside,
                     const std::vector<std::string> &outside) {
  const auto &token = unwrapped(node);
  if (token.name != "mo")
    return false;
  const auto &text = token.text;
  if (closingFence(text))
    return inside.empty() || inside.back() != text;
  if ((text == "|" || text == "‖") && inside.empty() && !outside.empty() &&
      outside.back() == text)
    return true;
  return inside.empty() &&
         (text == "+" || text == "=" || text == "," || text == ";" ||
          text == "<" || text == ">" || text == "≤" || text == "≥" ||
          text == "≠" || text == "−" || text == "-");
}
std::string sequence(const Xml &n) {
  std::string out;
  for (size_t i = 0; i < n.children.size(); i++) {
    const auto &c = *n.children[i];
    if (scripted(c) && !c.children.empty() && large(*c.children[0])) {
      std::string operand;
      std::vector<std::string> outside, inside;
      for (size_t prefix = 0; prefix < i; ++prefix)
        advanceFence(*n.children[prefix], outside);
      while (i + 1 < n.children.size()) {
        const auto &next = *n.children[i + 1];
        if (operandBoundary(next, inside, outside))
          break;
        advanceFence(next, inside);
        operand += emit(next);
        i++;
      }
      out += nary(c, operand);
    } else
      out += emit(c);
  }
  return out;
}
std::string wrapped(const std::string &tag, const std::string &body) {
  return "<m:" + tag + ">" + body + "</m:" + tag + ">";
}
std::string emit(const Xml &n) {
  const auto &name = n.name;
  auto child = [&](size_t i) { return emit(*n.children.at(i)); };
  if (name == "math" || name == "mrow" || name == "mstyle" ||
      name == "semantics" || name == "mtd")
    return sequence(n);
  if (name == "annotation" || name == "annotation-xml")
    return {};
  if (name == "mi" || name == "mn" || name == "mo" || name == "mtext")
    return run(n);
  if (name == "mspace") {
    if (std::stod(n.attr("width")) == 0.0)
      return {};
    return "<m:r><span style=\"" + style +
           ";font-style:normal\">&#160;</span></m:r>";
  }
  if (scripted(n) && !n.children.empty() && large(*n.children[0]))
    return nary(n, "");
  if (name == "mfrac")
    return wrapped("f", wrapped("fPr", (n.attr("linethickness") == "0"
                                            ? "<m:type m:val=\"noBar\"/>"
                                            : "") +
                                           control()) +
                            wrapped("num", child(0)) +
                            wrapped("den", child(1)));
  if (name == "msqrt" || name == "mroot")
    return wrapped(
        "rad",
        wrapped("radPr", (name == "msqrt" ? "<m:degHide m:val=\"on\"/>" : "") +
                             control()) +
            wrapped("deg", name == "mroot" ? child(1) : "") +
            wrapped("e", name == "msqrt" ? sequence(n) : child(0)));
  if (name == "msup" || name == "msub" || name == "msubsup") {
    const std::string tag = name == "msup"   ? "sSup"
                            : name == "msub" ? "sSub"
                                             : "sSubSup";
    std::string body = wrapped(tag + "Pr", control()) + wrapped("e", child(0));
    if (name != "msup")
      body += wrapped("sub", child(1));
    if (name != "msub")
      body += wrapped("sup", child(name == "msup" ? 1 : 2));
    return wrapped(tag, body);
  }
  if (name == "munder" || name == "mover") {
    const auto &mark = *n.children.at(1);
    const bool bar = mark.name == "mo" &&
                     (mark.text == "―" || mark.text == "¯" || mark.text == "_");
    if (bar)
      return wrapped(
          "bar",
          wrapped("barPr", "<m:pos m:val=\"" +
                               std::string(name == "munder" ? "bot" : "top") +
                               "\"/>" + control()) +
              wrapped("e", child(0)));
    if (name == "mover" && mark.name == "mo")
      return wrapped("acc",
                     wrapped("accPr", "<m:chr m:val=\"" + escape(mark.text) +
                                          "\"/>" + control()) +
                         wrapped("e", child(0)));
    const std::string tag = name == "munder" ? "limLow" : "limUpp";
    return wrapped(tag, wrapped(tag + "Pr", control()) +
                            wrapped("e", child(0)) + wrapped("lim", child(1)));
  }
  if (name == "mtable") {
    if (n.attr("columnalign") == "right left") {
      std::string rows;
      for (auto &r : n.children)
        rows += wrapped("e", sequence(*r));
      return wrapped(
          "eqArr", "<m:eqArrPr><m:baseJc m:val=\"left\"/></m:eqArrPr>" + rows);
    }
    return wrapped("m", wrapped("mPr", control()) + sequence(n));
  }
  if (name == "mtr") {
    std::string cells;
    for (auto &c : n.children)
      cells += wrapped("e", emit(*c));
    return wrapped("mr", cells);
  }
  if (name == "mphantom")
    return wrapped("phant",
                   wrapped("phantPr", "<m:show m:val=\"0\"/>" + control()) +
                       wrapped("e", sequence(n)));
  if (name == "menclose") {
    const auto notation = n.attr("notation");
    if (notation == "box")
      return wrapped("borderBox", wrapped("borderBoxPr", control()) +
                                      wrapped("e", sequence(n)));
    const std::map<std::string, std::string> strikes = {
        {"updiagonalstrike", "strikeBLTR"},
        {"downdiagonalstrike", "strikeTLBR"},
        {"horizontalstrike", "strikeH"},
        {"verticalstrike", "strikeV"}};
    auto found = strikes.find(notation);
    if (found == strikes.end())
      throw std::runtime_error("unsupported Office enclosure");
    return wrapped(
        "borderBox",
        wrapped("borderBoxPr",
                "<m:hideTop m:val=\"1\"/><m:hideBot m:val=\"1\"/><m:hideLeft "
                "m:val=\"1\"/><m:hideRight m:val=\"1\"/><m:" +
                    found->second + " m:val=\"1\"/>" + control()) +
            wrapped("e", sequence(n)));
  }
  throw std::runtime_error("unsupported Office element: " + name);
}
const Xml *outerTable(const Xml &n) {
  const auto &x = unwrapped(n);
  if (x.name == "math" && x.children.size() == 1)
    return outerTable(*x.children[0]);
  return x.name == "mtable" ? &x : nullptr;
}
bool leafRows(const Xml &table, std::vector<const Xml *> &rows) {
  if (table.attr("columnalign") != "left")
    return false;
  for (auto &row : table.children) {
    if (row->children.size() != 1)
      return false;
    const auto &cell = *row->children[0];
    const Xml *nested = outerTable(cell);
    if (!nested && cell.children.size() == 1)
      nested = outerTable(*cell.children[0]);
    if (nested) {
      if (!leafRows(*nested, rows))
        return false;
    } else
      rows.push_back(&cell);
  }
  return true;
}
} // namespace
std::string latex_to_office_html(const std::string &latex) {
  const auto mml = latex_to_mathml(latex, 18);
  if (mml.empty())
    throw std::runtime_error("invalid TeX");
  auto root = Parser(mml).parse();
  std::vector<const Xml *> rows;
  const auto table = outerTable(*root);
  std::string body;
  if (table && leafRows(*table, rows) && rows.size() >= 2) {
    for (auto row : rows)
      body += wrapped("e", emit(*row));
    body = wrapped("eqArr",
                   "<m:eqArrPr><m:baseJc m:val=\"left\"/></m:eqArrPr>" + body);
  } else
    body = emit(*root);
  const std::regex adjacent(
      R"rx(<m:r><span style="([^"]*)">([^<]*)</span></m:r><m:r><span style="\1">([^<]*)</span></m:r>)rx");
  for (;;) {
    auto merged = std::regex_replace(
        body, adjacent, "<m:r><span style=\"$1\">$2$3</span></m:r>");
    if (merged == body)
      break;
    body = std::move(merged);
  }
  const std::string math =
      "<m:oMathPara "
      "xmlns:m=\"http://schemas.microsoft.com/office/2004/12/"
      "omml\"><m:oMathParaPr><m:jc m:val=\"left\"/></m:oMathParaPr><m:oMath>" +
      body + "</m:oMath></m:oMathPara>";
  return "<p style=\"margin:0;font-size:18pt;font-family:Cambria "
         "Math\"><!--[if gte msEquation 12]>" +
         math + "<![endif]--><![if !msEquation]>" + mml + "<![endif]></p>";
}
} // namespace eqnedit
