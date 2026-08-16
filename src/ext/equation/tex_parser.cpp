/*
 * tex_parser.cpp -- LaTeX math -> shared node tree
 *
 * Recursive descent, one pass, no backtracking.  The grammar is small because
 * LaTeX math is: a sequence of atoms, each optionally carrying scripts, with a
 * handful of commands that take braced arguments.
 *
 *   seq   := atom*
 *   atom  := primary scripts?
 *   scripts := ('_' group | '^' group)*
 *   primary := char | '{' seq '}' | command
 *
 * A big operator is the one construct that reaches forward: its limits come
 * from the scripts and its body is the rest of the enclosing sequence, which
 * is what OMML's <m:nary> expects in <m:e>.
 */
#include "tex_parser.h"
#include "tex2mtef.h"
#include "mtef_common.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

namespace mtef {
namespace {

/* MTEF embellishment codes, shared with the LaTeX emitter's EMBELL_MAP. */
enum {
    EM_DOT = 2, EM_DDOT = 3, EM_TDOT = 4, EM_TILDE = 8, EM_HAT = 9,
    EM_RARROW = 11, EM_LARROW = 12, EM_BARROW = 13, EM_OBAR = 17,
};

bool is_letter(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}
bool is_digit(char c) { return c >= '0' && c <= '9'; }

/* Typeface for a symbol that came from the command table. */
int typeface_for_code(uint16_t code) {
    if (code >= 0x0391 && code <= 0x03A9) return TF_UCGREEK;
    if (code >= 0x03B1 && code <= 0x03C9) return TF_LCGREEK;
    return TF_SYMBOL;
}

std::unique_ptr<CharNode> make_char(int typeface, uint16_t code, char ascii = 0) {
    auto c = std::make_unique<CharNode>();
    c->typeface = typeface;
    c->charCode = code;
    c->ch = ascii;
    return c;
}

/* Functions LaTeX sets upright and spaces as operators. */
bool is_function_name(const std::string& name) {
    static const char* kNames[] = {
        "arccos","arcsin","arctan","arg","cos","cosh","cot","coth","csc",
        "deg","det","dim","exp","gcd","hom","inf","ker","lg","lim","liminf",
        "limsup","ln","log","max","min","Pr","sec","sin","sinh","sup","tan",
        "tanh","curl","div","grad","rot","sgn","tr","diag","mod","Re","Im",
        nullptr
    };
    for (int i = 0; kNames[i]; ++i)
        if (name == kNames[i]) return true;
    return false;
}

/* Delimiter token -> fence selector.  -1 when the token is not a delimiter. */
int fence_selector(const std::string& tok) {
    if (tok == "(" || tok == ")") return tmPAREN;
    if (tok == "[" || tok == "]") return tmBRACK;
    if (tok == "\\{" || tok == "\\}") return tmBRACE;
    if (tok == "|") return tmBAR;
    if (tok == "\\|" || tok == "\\Vert") return tmDBAR;
    if (tok == "\\langle" || tok == "\\rangle") return tmANGLE;
    if (tok == "\\lfloor" || tok == "\\rfloor") return tmFLOOR;
    if (tok == "\\lceil" || tok == "\\rceil") return tmCEIL;
    return -1;
}

class TexParser {
public:
    explicit TexParser(std::string s) : s_(std::move(s)), p_(0) {}

    std::unique_ptr<LineNode> parse() {
        auto root = std::make_unique<LineNode>();
        root->children = parse_seq(kTop);
        return root;
    }

private:
    /* Stop conditions for a sequence.  A matrix cell also ends at & and \\. */
    enum Ctx { kTop, kBrace, kCell };

    const std::string s_;    /* owned: the caller passes a stripped temporary */
    size_t p_;

    /* ---- lexer ---------------------------------------------------------- */

    bool eof() const { return p_ >= s_.size(); }
    char peek(size_t off = 0) const {
        return (p_ + off < s_.size()) ? s_[p_ + off] : '\0';
    }
    void skip_space() {
        while (!eof() && (s_[p_] == ' ' || s_[p_] == '\t' || s_[p_] == '\n' ||
                          s_[p_] == '\r'))
            ++p_;
    }

    /* Read a command starting at the backslash: \alpha, \{, \\ . */
    std::string read_command() {
        std::string cmd = "\\";
        ++p_;                                  /* the backslash */
        if (eof()) return cmd;
        if (!is_letter(peek())) { cmd += s_[p_++]; return cmd; }
        while (!eof() && is_letter(peek())) cmd += s_[p_++];
        return cmd;
    }

    /* Peek the next command without consuming it. */
    std::string peek_command() {
        size_t save = p_;
        skip_space();
        std::string cmd;
        if (peek() == '\\') cmd = read_command();
        p_ = save;
        return cmd;
    }

    bool eat_command(const char* name) {
        size_t save = p_;
        skip_space();
        if (peek() != '\\') { p_ = save; return false; }
        std::string cmd = read_command();
        if (cmd == name) return true;
        p_ = save;
        return false;
    }

    bool at_stop(Ctx ctx) {
        skip_space();
        if (eof()) return true;
        if (ctx == kBrace && peek() == '}') return true;
        if (ctx == kCell) {
            if (peek() == '&' || peek() == '}') return true;
            if (peek() == '\\' && peek(1) == '\\') return true;
            std::string cmd = peek_command();
            if (cmd == "\\end") return true;
        }
        /* \right closes the fence that opened this sequence. */
        if (peek() == '\\') {
            std::string cmd = peek_command();
            if (cmd == "\\right" || cmd == "\\end") return true;
        }
        return false;
    }

    /* ---- sequences ------------------------------------------------------ */

    NodeList parse_seq(Ctx ctx) {
        NodeList out;
        while (!at_stop(ctx)) {
            size_t before = p_;
            NodePtr a = parse_atom_with_scripts(ctx);
            if (!a) {
                if (p_ == before) ++p_;        /* never spin on bad input */
                continue;
            }
            /* A big operator takes the rest of this sequence as its body. */
            if (a->tag() == Node::kIntegral) {
                static_cast<IntegralNode&>(*a).body = parse_seq(ctx);
                out.push_back(std::move(a));
                break;
            }
            if (a->tag() == Node::kBigOp) {
                static_cast<BigOpNode&>(*a).body = parse_seq(ctx);
                out.push_back(std::move(a));
                break;
            }
            out.push_back(std::move(a));
        }
        return out;
    }

    /* A braced argument, or a single atom when the argument is unbraced. */
    NodeList parse_arg() {
        skip_space();
        if (peek() == '{') {
            ++p_;
            NodeList l = parse_seq(kBrace);
            skip_space();
            if (peek() == '}') ++p_;
            return l;
        }
        NodeList l;
        if (NodePtr a = parse_atom(kTop)) l.push_back(std::move(a));
        return l;
    }

    /* An optional [..] argument, e.g. the index of \sqrt. */
    bool parse_optional(NodeList& out) {
        skip_space();
        if (peek() != '[') return false;
        ++p_;
        while (!eof() && peek() != ']') {
            size_t before = p_;
            if (NodePtr a = parse_atom_with_scripts(kTop)) out.push_back(std::move(a));
            else if (p_ == before) ++p_;
        }
        if (peek() == ']') ++p_;
        return true;
    }

    /* ---- atoms ---------------------------------------------------------- */

    NodePtr parse_atom_with_scripts(Ctx ctx) {
        NodePtr base = parse_atom(ctx);
        if (!base) return nullptr;

        if (base->tag() == Node::kIntegral || base->tag() == Node::kBigOp) {
            attach_limits(*base);
            return base;
        }

        std::unique_ptr<ScriptNode> sc;
        for (;;) {
            skip_space();
            char c = peek();
            if (c != '_' && c != '^') break;
            ++p_;
            if (!sc) {
                sc = std::make_unique<ScriptNode>();
                sc->base.push_back(std::move(base));
            }
            if (c == '_') { sc->sub = parse_arg(); sc->hasSub = true; }
            else          { sc->sup = parse_arg(); sc->hasSup = true; }
        }
        if (sc) return sc;
        return base;
    }

    /* Limits on a big operator go in its own slots.  Sums stack their limits
     * above and below by default, integrals set them beside; \limits and
     * \nolimits override that. */
    void attach_limits(Node& n) {
        const bool is_int = (n.tag() == Node::kIntegral);
        bool stacked = !is_int;
        for (;;) {
            skip_space();
            if (eat_command("\\limits"))   { stacked = true;  continue; }
            if (eat_command("\\nolimits")) { stacked = false; continue; }
            char c = peek();
            if (c != '_' && c != '^') break;
            ++p_;
            NodeList arg = parse_arg();
            if (n.tag() == Node::kIntegral) {
                auto& i = static_cast<IntegralNode&>(n);
                if (c == '_') { i.lower = std::move(arg); i.hasLower = true; }
                else          { i.upper = std::move(arg); i.hasUpper = true; }
            } else {
                auto& b = static_cast<BigOpNode&>(n);
                if (c == '_') { b.lower = std::move(arg); b.hasLower = true; }
                else          { b.upper = std::move(arg); b.hasUpper = true; }
            }
        }
        if (is_int) static_cast<IntegralNode&>(n).hasLimits = stacked;
        else        static_cast<BigOpNode&>(n).hasLimits = stacked;
    }

    NodePtr parse_atom(Ctx ctx) {
        skip_space();
        if (eof()) return nullptr;
        char c = peek();

        if (c == '{') {
            ++p_;
            auto line = std::make_unique<LineNode>();
            line->children = parse_seq(kBrace);
            skip_space();
            if (peek() == '}') ++p_;
            return line;
        }
        if (c == '}' || c == '&') return nullptr;
        if (c == '\\') return parse_command(ctx);
        if (c == '_' || c == '^') return nullptr;   /* handled by the caller */

        ++p_;
        if (is_digit(c))  return make_char(TF_NUMBER,   uint16_t(c), c);
        if (is_letter(c)) return make_char(TF_VARIABLE, uint16_t(c), c);
        /* Office sets a binary minus as U+2212, not the hyphen on the key. */
        if (c == '-')     return make_char(TF_SYMBOL, 0x2212, c);
        return make_char(TF_SYMBOL, uint16_t((unsigned char)c), c);
    }

    /* ---- commands ------------------------------------------------------- */

    NodePtr parse_command(Ctx ctx) {
        size_t save = p_;
        std::string cmd = read_command();

        /* Spacing and layout commands carry no glyph. */
        if (cmd == "\\," || cmd == "\\;" || cmd == "\\:" || cmd == "\\!" ||
            cmd == "\\ " || cmd == "\\quad" || cmd == "\\qquad" ||
            cmd == "\\displaystyle" || cmd == "\\textstyle" ||
            cmd == "\\limits" || cmd == "\\nolimits" || cmd == "\\left" ||
            cmd == "\\right") {
            if (cmd == "\\left")  { p_ = save; return parse_fence(); }
            if (cmd == "\\right") { p_ = save; return nullptr; }
            return nullptr;
        }

        if (cmd == "\\frac" || cmd == "\\dfrac" || cmd == "\\tfrac" ||
            cmd == "\\cfrac") {
            auto f = std::make_unique<FracNode>();
            f->numer = parse_arg();
            f->denom = parse_arg();
            f->display = (cmd == "\\dfrac");
            return f;
        }
        if (cmd == "\\sqrt") {
            auto s = std::make_unique<SqrtNode>();
            s->hasIndex = parse_optional(s->index);
            s->content = parse_arg();
            return s;
        }
        if (cmd == "\\begin") return parse_environment();

        if (int sel = big_op_selector(cmd), isint = integral_selector(cmd);
            sel >= 0 || isint >= 0) {
            if (isint >= 0) {
                auto i = std::make_unique<IntegralNode>();
                i->selector = isint;
                return i;
            }
            auto b = std::make_unique<BigOpNode>();
            b->selector = sel;
            return b;
        }

        if (int em = accent_code(cmd); em >= 0) {
            auto e = std::make_unique<EmbellNode>();
            e->embellType = em;
            e->content = parse_arg();
            return e;
        }
        if (cmd == "\\overline" || cmd == "\\underline") {
            auto d = std::make_unique<DecorationNode>();
            d->selector = (cmd == "\\overline") ? tmOBAR : tmUBAR;
            d->content = parse_arg();
            return d;
        }
        if (cmd == "\\overrightarrow" || cmd == "\\overleftarrow" ||
            cmd == "\\overleftrightarrow") {
            auto d = std::make_unique<DecorationNode>();
            d->selector = (cmd == "\\overrightarrow") ? tmRARROW
                        : (cmd == "\\overleftarrow")  ? tmLARROW : tmBARROW;
            d->content = parse_arg();
            return d;
        }

        if (cmd == "\\text" || cmd == "\\mathrm" || cmd == "\\textrm" ||
            cmd == "\\operatorname")
            return styled_group(TF_TEXT);
        if (cmd == "\\mathbf" || cmd == "\\bm" || cmd == "\\boldsymbol")
            return styled_group(TF_VECTOR);
        if (cmd == "\\mathit")
            return styled_group(TF_VARIABLE);

        /* \sin, \log, ... -- upright, no braces. */
        if (cmd.size() > 1 && is_function_name(cmd.substr(1))) {
            auto line = std::make_unique<LineNode>();
            for (size_t i = 1; i < cmd.size(); ++i)
                line->children.push_back(
                    make_char(TF_FUNCTION, uint16_t(cmd[i]), cmd[i]));
            return line;
        }

        /* Symbols share the table the MTEF writer uses. */
        int code = tex_command_to_unicode(cmd.c_str());
        if (code >= 0) {
            auto ch = make_char(typeface_for_code(uint16_t(code)), uint16_t(code));
            ch->latex = cmd;
            return ch;
        }

        /* Escaped literals: \{ \} \% \& \$ \# \_ */
        if (cmd.size() == 2 && !is_letter(cmd[1]))
            return make_char(TF_SYMBOL, uint16_t((unsigned char)cmd[1]), cmd[1]);

        /* Unknown: show the name rather than losing the equation. */
        auto line = std::make_unique<LineNode>();
        for (size_t i = 1; i < cmd.size(); ++i)
            line->children.push_back(make_char(TF_TEXT, uint16_t(cmd[i]), cmd[i]));
        return line;
    }

    /* A braced group whose characters all take one typeface. */
    NodePtr styled_group(int typeface) {
        NodeList inner = parse_arg();
        retype(inner, typeface);
        auto line = std::make_unique<LineNode>();
        line->children = std::move(inner);
        return line;
    }

    static void retype(NodeList& list, int typeface) {
        for (auto& n : list) {
            if (!n) continue;
            if (n->tag() == Node::kChar)
                static_cast<CharNode&>(*n).typeface = typeface;
            else if (n->tag() == Node::kLine)
                retype(static_cast<LineNode&>(*n).children, typeface);
        }
    }

    static int accent_code(const std::string& cmd) {
        if (cmd == "\\hat" || cmd == "\\widehat")   return EM_HAT;
        if (cmd == "\\tilde" || cmd == "\\widetilde") return EM_TILDE;
        if (cmd == "\\dot")   return EM_DOT;
        if (cmd == "\\ddot")  return EM_DDOT;
        if (cmd == "\\dddot") return EM_TDOT;
        if (cmd == "\\vec")   return EM_RARROW;
        if (cmd == "\\bar")   return EM_OBAR;
        return -1;
    }

    static int integral_selector(const std::string& cmd) {
        if (cmd == "\\int")    return tmSINT;
        if (cmd == "\\iint")   return tmDINT;
        if (cmd == "\\iiint")  return tmTINT;
        if (cmd == "\\oint")   return tmSSINT;
        if (cmd == "\\oiint")  return tmDSINT;
        if (cmd == "\\oiiint") return tmTSINT;
        return -1;
    }

    static int big_op_selector(const std::string& cmd) {
        if (cmd == "\\sum")     return tmSUM;
        if (cmd == "\\prod")    return tmPROD;
        if (cmd == "\\coprod")  return tmCOPROD;
        if (cmd == "\\bigcup")  return tmUNION;
        if (cmd == "\\bigcap")  return tmINTER;
        return -1;
    }

    /* ---- fences --------------------------------------------------------- */

    /* Read the delimiter token after \left or \right. */
    std::string read_delimiter() {
        skip_space();
        if (eof()) return ".";
        if (peek() == '\\') return read_command();
        char c = s_[p_++];
        return std::string(1, c);
    }

    NodePtr parse_fence() {
        eat_command("\\left");
        std::string open = read_delimiter();
        auto f = std::make_unique<FenceNode>();
        int sel = fence_selector(open);
        f->selector = (sel >= 0) ? sel : tmPAREN;
        f->content = parse_seq(kBrace);

        std::string close = ".";
        if (eat_command("\\right")) close = read_delimiter();

        /* variation 1 = left delimiter only, 2 = right only */
        if (open == "." && close != ".") {
            f->variation = 2;
            int cs = fence_selector(close);
            if (cs >= 0) f->selector = cs;
        } else if (close == "." && open != ".") {
            f->variation = 1;
        }
        return f;
    }

    /* ---- environments --------------------------------------------------- */

    std::string read_env_name() {
        skip_space();
        std::string name;
        if (peek() != '{') return name;
        ++p_;
        while (!eof() && peek() != '}') name += s_[p_++];
        if (peek() == '}') ++p_;
        return name;
    }

    NodePtr parse_environment() {
        std::string env = read_env_name();

        int fence = -1;
        if (env == "pmatrix") fence = tmPAREN;
        else if (env == "bmatrix") fence = tmBRACK;
        else if (env == "Bmatrix") fence = tmBRACE;
        else if (env == "vmatrix") fence = tmBAR;
        else if (env == "Vmatrix") fence = tmDBAR;

        const bool is_matrix = (env == "matrix" || fence >= 0);
        std::vector<std::vector<NodeList>> rows;
        rows.emplace_back();

        for (;;) {
            skip_space();
            if (eof()) break;
            if (peek_command() == "\\end") { eat_command("\\end"); read_env_name(); break; }

            rows.back().push_back(parse_seq(kCell));

            skip_space();
            if (peek() == '&') { ++p_; continue; }
            if (peek() == '\\' && peek(1) == '\\') { p_ += 2; rows.emplace_back(); continue; }
            if (peek_command() == "\\end") { eat_command("\\end"); read_env_name(); break; }
            if (eof()) break;
            ++p_;                                   /* never spin */
        }
        while (!rows.empty() && rows.back().empty()) rows.pop_back();

        size_t cols = 0;
        for (const auto& r : rows) cols = std::max(cols, r.size());

        if (is_matrix) {
            auto m = std::make_unique<MatrixNode>();
            m->rows = int(rows.size());
            m->cols = int(cols);
            for (auto& r : rows) {
                for (size_t c = 0; c < cols; ++c) {
                    auto cell = std::make_unique<LineNode>();
                    if (c < r.size()) cell->children = std::move(r[c]);
                    m->elements.push_back(std::move(cell));
                }
            }
            if (fence < 0) return m;
            auto f = std::make_unique<FenceNode>();
            f->selector = fence;
            f->content.push_back(std::move(m));
            return f;
        }

        /* cases / aligned / gathered -- one stacked line per row. */
        auto pile = std::make_unique<PileNode>();
        pile->ncols = int(cols);
        for (auto& r : rows) {
            auto line = std::make_unique<LineNode>();
            for (auto& cell : r)
                for (auto& n : cell) line->children.push_back(std::move(n));
            pile->lines.push_back(std::move(line));
        }
        if (env != "cases") return pile;
        auto f = std::make_unique<FenceNode>();
        f->selector = tmBRACE;
        f->variation = 1;
        f->content.push_back(std::move(pile));
        return f;
    }
};

/* Strip $...$, $$...$$, \(..\), \[..\] wrappers. */
std::string strip_delimiters(const std::string& in) {
    std::string s = in;
    size_t a = s.find_first_not_of(" \t\r\n");
    size_t b = s.find_last_not_of(" \t\r\n");
    if (a == std::string::npos) return std::string();
    s = s.substr(a, b - a + 1);

    auto strip = [&s](const char* open, const char* close) {
        size_t no = strlen(open), nc = strlen(close);
        if (s.size() >= no + nc && s.compare(0, no, open) == 0 &&
            s.compare(s.size() - nc, nc, close) == 0)
            s = s.substr(no, s.size() - no - nc);
    };
    strip("$$", "$$");
    strip("\\[", "\\]");
    strip("\\(", "\\)");
    strip("$", "$");
    return s;
}

}  // namespace

std::unique_ptr<LineNode> parse_latex(const std::string& latex) {
    TexParser p(strip_delimiters(latex));
    return p.parse();
}

}  // namespace mtef
