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
 * from scripts and its body is the rest of the enclosing sequence.
 */
#include "tex_parser.h"
#include "math_symbols.h"

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>

namespace eqnedit {
namespace {

bool is_letter(char c) {
    return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}
bool is_digit(char c) { return c >= '0' && c <= '9'; }

std::unique_ptr<CharNode> make_char(int typeface, uint32_t code, char ascii = 0) {
    auto c = std::make_unique<CharNode>();
    c->typeface = typeface;
    c->charCode = code;
    c->ch = ascii;
    return c;
}

/* Functions LaTeX sets upright and spaces as operators.  Names that also
 * exist as glyph commands (\div, \Re, \Im) belong in the symbol table, not
 * here: parse_command() consults that table first. */
bool is_function_name(const std::string& name) {
    static const char* kNames[] = {
        "arccos","arcsin","arctan","arg","cos","cosh","cot","coth","csc",
        "deg","det","dim","exp","gcd","hom","inf","ker","lg","lim","liminf",
        "limsup","ln","log","max","min","Pr","sec","sin","sinh","sup","tan",
        "tanh","curl","grad","rot","sgn","tr","diag","mod",
        /* \lim is handled as its own node so its subscript sets under it in
         * display style, so it is deliberately not in this upright-name list. */
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

/* The glyph a \left / \right / \middle delimiter token stands for.  Null for
 * the invisible delimiter "." and for anything unrecognisable. */
NodePtr delimiter_char(const std::string& tok) {
    if (tok.empty() || tok == ".") return nullptr;
    const int code = latex_symbol_codepoint(tok == "\\|" ? "\\Vert" : tok);
    if (code >= 0) {
        auto ch = make_char(typeface_for_code(uint32_t(code)), uint32_t(code));
        ch->latex = tok;
        return ch;
    }
    if (tok[0] == '\\' && tok.size() == 2)          /* \|  \{  \}  */
        return make_char(TF_SYMBOL, uint32_t(uint8_t(tok[1])), tok[1]);
    if (tok.size() == 1)
        return make_char(TF_SYMBOL, uint32_t(uint8_t(tok[0])), tok[0]);
    return nullptr;
}

class TexParser {
public:
    explicit TexParser(std::string s) : s_(std::move(s)), p_(0) {}

    std::unique_ptr<LineNode> parse() {
        auto root = std::make_unique<LineNode>();
        root->children = parse_seq(kTop);
        return root;
    }

    bool depth_exceeded() const { return depthExceeded_; }

private:
    /* Stop conditions for a sequence.  A matrix cell also ends at & and \\;
     * a fence body also ends at \middle, so the separator stays visible to
     * parse_fence() instead of being swallowed as one more atom. */
    enum Ctx { kTop, kBrace, kCell, kFenceBody };

    const std::string s_;    /* owned: the caller passes a stripped temporary */
    size_t p_;
    int depth_ = 0;          /* recursion depth, to bound pathological nesting */
    bool depthExceeded_ = false;
    int literalSpaceDepth_ = 0; /* inside text-like braced arguments */

    /* Parse, layout, and emit are all recursive over the tree, so a deeply
     * nested paste -- \sqrt{\sqrt{...}} thousands deep -- overflowed the
     * stack and took the process down (measured: fine to ~1200, crash by
     * ~1400).  Equation Editor 3.0 caps nesting and refuses with a message
     * rather than crashing (its string 16044); this is the same finite cap,
     * set well above any real equation and well below the crash.  Content
     * past the cap is dropped, not built. */
    /* Consume the body of the current group without building it, so the
     * parser still terminates and brace nesting stays balanced. */
    void skip_group_body(Ctx ctx) {
        int brace = 0;
        while (!eof()) {
            char c = peek();
            if (brace == 0 && at_stop(ctx)) return;
            if (c == '{') { ++brace; ++p_; }
            else if (c == '}') { if (brace == 0) return; --brace; ++p_; }
            else ++p_;
        }
    }

    /* ---- lexer ---------------------------------------------------------- */

    bool eof() const { return p_ >= s_.size(); }
    char peek(size_t off = 0) const {
        return (p_ + off < s_.size()) ? s_[p_ + off] : '\0';
    }

    /* Consume one UTF-8 scalar.  Invalid sequences become U+FFFD, consuming
     * one byte, so malformed pasted text can never stall the parser. */
    uint32_t read_codepoint() {
        if (eof()) return 0;
        const unsigned char lead = static_cast<unsigned char>(s_[p_]);
        uint32_t cp = 0;
        size_t count = 1;
        uint32_t minimum = 0;
        if (lead < 0x80) {
            ++p_;
            return lead;
        }
        if ((lead & 0xE0) == 0xC0) {
            cp = lead & 0x1F; count = 2; minimum = 0x80;
        } else if ((lead & 0xF0) == 0xE0) {
            cp = lead & 0x0F; count = 3; minimum = 0x800;
        } else if ((lead & 0xF8) == 0xF0) {
            cp = lead & 0x07; count = 4; minimum = 0x10000;
        } else {
            ++p_;
            return 0xFFFD;
        }
        if (p_ + count > s_.size()) {
            ++p_;
            return 0xFFFD;
        }
        for (size_t i = 1; i < count; ++i) {
            const unsigned char next = static_cast<unsigned char>(s_[p_ + i]);
            if ((next & 0xC0) != 0x80) {
                ++p_;
                return 0xFFFD;
            }
            cp = (cp << 6) | (next & 0x3F);
        }
        p_ += count;
        if (cp < minimum || cp > 0x10FFFF ||
            (cp >= 0xD800 && cp <= 0xDFFF))
            return 0xFFFD;
        return cp;
    }
    void skip_space() {
        for (;;) {
            if (literalSpaceDepth_ == 0) {
                while (!eof() &&
                       (s_[p_] == ' ' || s_[p_] == '\t' ||
                        s_[p_] == '\n' || s_[p_] == '\r'))
                    ++p_;
            }
            /* An unescaped percent starts a TeX comment.  Escaped `\\%`
             * reaches parse_command() and remains a visible percent sign. */
            if (eof() || s_[p_] != '%') return;
            while (!eof() && s_[p_] != '\n' && s_[p_] != '\r') ++p_;
            /* A TeX comment consumes its line ending even in a text-like
             * argument, where ordinary spaces are otherwise content. */
            if (literalSpaceDepth_ > 0) {
                if (!eof() && s_[p_] == '\r') ++p_;
                if (!eof() && s_[p_] == '\n') ++p_;
            }
        }
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
        if ((ctx == kBrace || ctx == kFenceBody) && peek() == '}') return true;
        if (ctx == kFenceBody && peek() == '\\' &&
            peek_command() == "\\middle")
            return true;
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
        if (depth_ >= kMaxNestingDepth && !at_stop(ctx)) {
            /* The empty slot inside the deepest accepted template still
             * enters parse_seq once.  It consumes no structural level, so it
             * is valid; only actual content beyond the cap is rejected. */
            depthExceeded_ = true;
            skip_group_body(ctx);
            return out;
        }
        ++depth_;
        struct Pop { int& d; ~Pop() { --d; } } pop{depth_};
        while (!at_stop(ctx)) {
            size_t before = p_;
            NodePtr a = parse_atom_with_scripts();
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
        if (NodePtr a = parse_atom()) l.push_back(std::move(a));
        return l;
    }

    /* An optional [..] argument, e.g. the index of \sqrt. */
    bool parse_optional(NodeList& out) {
        skip_space();
        if (peek() != '[') return false;
        ++p_;
        for (;;) {
            /* emitLim() terminates its control word with a space.  Test the
             * optional delimiter after ignoring that TeX whitespace, or the
             * next atom consumes `]` and the radicand into the index. */
            skip_space();
            if (eof() || peek() == ']') break;
            size_t before = p_;
            if (NodePtr a = parse_atom_with_scripts()) out.push_back(std::move(a));
            else if (p_ == before) ++p_;
        }
        if (peek() == ']') ++p_;
        return true;
    }

    /* ---- atoms ---------------------------------------------------------- */

    NodePtr parse_atom_with_scripts() {
        NodePtr base = parse_atom();
        if (!base) return nullptr;

        if (base->tag() == Node::kIntegral || base->tag() == Node::kBigOp) {
            attach_limits(*base);
            return base;
        }

        /* \lim's following _{...} is its condition, set underneath, not a
         * subscript beside it. */
        if (base->tag() == Node::kLim) {
            skip_space();
            if (peek() == '_') {
                ++p_;
                static_cast<LimNode&>(*base).content = parse_arg();
            }
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

    NodePtr parse_atom() {
        skip_space();
        if (eof()) return nullptr;
        char c = peek();

        /* In text-like arguments TeX turns any run of ordinary source
         * whitespace into one visible word space.  Normal math parsing has
         * already skipped it before reaching this point. */
        if (literalSpaceDepth_ > 0 &&
            (c == ' ' || c == '\t' || c == '\n' || c == '\r')) {
            do { ++p_; }
            while (!eof() &&
                   (peek() == ' ' || peek() == '\t' ||
                    peek() == '\n' || peek() == '\r'));
            return make_char(TF_TEXT, uint32_t(' '), ' ');
        }

        if (c == '{') {
            ++p_;
            auto line = std::make_unique<LineNode>();
            line->children = parse_seq(kBrace);
            skip_space();
            if (peek() == '}') ++p_;
            return line;
        }
        if (c == '}' || c == '&') return nullptr;
        if (c == '\\') return parse_command();
        if (c == '_' || c == '^') return nullptr;   /* handled by the caller */

        /* In TeX source a bare tilde is a non-breaking space, not the
         * mathematical relation U+223C.  Preserve the source spelling while
         * giving the canvas and Office Math a real, inkless advance. */
        if (c == '~') {
            ++p_;
            auto space = make_char(TF_SPACE, 0);
            space->latex = "~";
            return space;
        }

        uint32_t cp = read_codepoint();
        if (is_digit(c))  return make_char(TF_NUMBER,   cp, c);
        if (is_letter(c)) return make_char(TF_VARIABLE, cp, c);
        /* Render a binary minus as U+2212, not the keyboard hyphen. */
        if (c == '-')     return make_char(TF_SYMBOL, 0x2212, c);
        if (cp >= 0x80) {
            const int tf = typeface_for_code(cp);
            /* Known mathematical Unicode is normalized to its TeX command.
             * Other Unicode is text, preserving Japanese and supplementary
             * characters in a stable \text{...} group. */
            return make_char(tf == TF_SYMBOL && !is_latex_symbol_codepoint(cp)
                                 ? TF_TEXT : tf,
                             cp);
        }
        return make_char(TF_SYMBOL, cp, c);
    }

    /* ---- commands ------------------------------------------------------- */

    NodePtr parse_command() {
        size_t save = p_;
        /* A pasted backslash can land immediately before Unicode text.  TeX
         * control words are ASCII, so consume the complete scalar instead of
         * treating its first UTF-8 byte as an escaped character. */
        if (peek() == '\\' &&
            static_cast<unsigned char>(peek(1)) >= 0x80) {
            ++p_;
            uint32_t cp = read_codepoint();
            int tf = typeface_for_code(cp);
            if (tf == TF_SYMBOL && !is_latex_symbol_codepoint(cp)) tf = TF_TEXT;
            return make_char(tf, cp);
        }
        std::string cmd = read_command();

        /* Explicit spacing.  These used to be discarded, so `a \quad b` came
         * back as `ab`: spacing an author had deliberately written was lost
         * the first time the file was opened and saved. */
        if (cmd == "\\," || cmd == "\\;" || cmd == "\\:" || cmd == "\\!" ||
            cmd == "\\ " || cmd == "\\quad" || cmd == "\\qquad") {
            auto space = make_char(TF_SPACE, 0);
            space->latex = cmd;
            return space;
        }

        /* Layout commands carry neither glyph nor width. */
        if (cmd == "\\displaystyle" || cmd == "\\textstyle" ||
            cmd == "\\limits" || cmd == "\\nolimits" || cmd == "\\left" ||
            cmd == "\\right") {
            if (cmd == "\\left")  { p_ = save; return parse_fence(); }
            if (cmd == "\\right") { p_ = save; return nullptr; }
            return nullptr;
        }

        /* \middle<delim> is a stretching separator inside \left...\right.
         * The tree has no such node, so keep the delimiter itself: an
         * ordinary bar is right for a bra-ket, and far better than the word
         * "middle" that the unknown-command fallback used to print. */
        if (cmd == "\\middle") return delimiter_char(read_delimiter());

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

        /* Document metadata controls do not draw mathematical content.  A
         * paper equation pasted with these commands must not display the
         * words "label" or "nonumber" on the canvas. */
        if (cmd == "\\label" || cmd == "\\tag") {
            skip_space();
            if (cmd == "\\tag" && peek() == '*') ++p_;
            (void)parse_arg();
            return nullptr;
        }
        if (cmd == "\\nonumber" || cmd == "\\notag") return nullptr;

        /* \lim takes its condition underneath in display style, so it is its
         * own node; parse_atom_with_scripts pulls the following _{...} into
         * that node rather than making an ordinary subscript beside it. */
        if (cmd == "\\lim") return std::make_unique<LimNode>();

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
        /* Horizontal braces.  The emitter has always written these, but
         * nothing ever read them back, so a horizontal brace saved and
         * reopened returned as the literal text "overbrace".  TeX labels the
         * brace on its own side: \overbrace{x}^{n}, \underbrace{x}_{n}. */
        if (cmd == "\\overbrace" || cmd == "\\underbrace") {
            const bool over = (cmd == "\\overbrace");
            auto b = std::make_unique<BraceDecoNode>();
            b->selector = over ? tmUHBRACE : tmLHBRACE;
            b->content = parse_arg();
            skip_space();
            if (peek() == (over ? '^' : '_')) { ++p_; b->label = parse_arg(); }
            return b;
        }
        /* The strike-through embellishment emits \cancel; with no reader it
         * came back as the text "cancel" followed by its own content. */
        if (cmd == "\\cancel") {
            auto e = std::make_unique<EmbellNode>();
            e->embellType = EM_MBAR;
            e->content = parse_arg();
            return e;
        }
        if (cmd == "\\overset" || cmd == "\\underset") {
            auto o = std::make_unique<OversetNode>();
            o->under = (cmd == "\\underset");
            o->over = parse_arg();
            o->base = parse_arg();
            return o;
        }

        if (cmd == "\\text" || cmd == "\\textrm")
            return styled_group(TF_TEXT);
        if (cmd == "\\mathrm")
            return styled_group(TF_ROMAN);
        if (cmd == "\\textbackslash") return make_char(TF_TEXT, '\\', '\\');
        if (cmd == "\\textasciicircum") return make_char(TF_TEXT, '^', '^');
        if (cmd == "\\textasciitilde") return make_char(TF_TEXT, '~', '~');
        if (cmd == "\\operatorname") return styled_group(TF_FUNCTION);
        if (cmd == "\\mathbf")
            return styled_group(TF_VECTOR);
        if (cmd == "\\bm" || cmd == "\\boldsymbol")
            return styled_group(TF_BOLD_SYMBOL);
        if (cmd == "\\mathit")
            return styled_group(TF_MATH_ITALIC);
        if (cmd == "\\mathsf")
            return styled_group(TF_MATH_SANS);
        if (cmd == "\\mathtt")
            return styled_group(TF_MATH_MONO);
        if (cmd == "\\mathcal")
            return styled_group(TF_MATH_SCRIPT);
        if (cmd == "\\mathbb")
            return styled_group(TF_MATH_DOUBLE);
        if (cmd == "\\mathfrak")
            return styled_group(TF_MATH_FRAKTUR);
        if (cmd == "\\mathnormal") {
            auto line = std::make_unique<LineNode>();
            line->children = parse_arg();
            return line;
        }

        /* Musical accidentals are also the standard differential-geometry
         * flat/sharp operators.  They are used inside the shared native/Web
         * geometry palette as ^{\\flat} and ^{\\sharp}.  Treating every
         * otherwise-unknown control word as invisible must not make these
         * two supported palette entries disappear. */
        if (cmd == "\\flat" || cmd == "\\sharp") {
            const uint32_t code = cmd == "\\flat" ? 0x266D : 0x266F;
            auto ch = make_char(TF_SYMBOL, code);
            ch->latex = cmd;
            return ch;
        }

        /* Named glyphs win over operator names.  \div is the division sign,
         * not the word "div"; the same holds for \Re and \Im.  Keeping the
         * command in `latex` is what lets the emitter write back exactly
         * what was read instead of guessing from the code point. */
        /* \not applied to a relation is one negated character, and it is how
         * LaTeX spells the ones that have no command of their own.  Fold the
         * pair here so \not\subset comes back as U+2284 rather than as a
         * combining mark followed by a subset sign. */
        if (cmd == "\\not") {
            const size_t afterNot = p_;
            skip_space();
            if (peek() == '\\') {
                const std::string relation = read_command();
                if (latex_symbol_codepoint(cmd + relation) >= 0) {
                    const int code = latex_symbol_codepoint(cmd + relation);
                    auto ch = make_char(typeface_for_code(uint32_t(code)),
                                        uint32_t(code));
                    ch->latex = cmd + relation;
                    return ch;
                }
            }
            p_ = afterNot;
        }

        /* \| is TeX's short spelling of \Vert; canonicalise so the double bar
         * survives instead of decaying to a single one. */
        const std::string named = (cmd == "\\|") ? std::string("\\Vert") : cmd;
        int code = latex_symbol_codepoint(named);
        if (code >= 0) {
            auto ch = make_char(typeface_for_code(uint32_t(code)), uint32_t(code));
            ch->latex = named;
            return ch;
        }

        /* \sin, \log, ... -- upright, no braces. */
        if (cmd.size() > 1 && is_function_name(cmd.substr(1))) {
            auto line = std::make_unique<LineNode>();
            for (size_t i = 1; i < cmd.size(); ++i)
                line->children.push_back(
                    make_char(TF_FUNCTION, uint32_t(uint8_t(cmd[i])), cmd[i]));
            return line;
        }

        /* Escaped literals: \{ \} \% \& \$ \# \_ */
        if (cmd.size() == 2 && !is_letter(cmd[1])) {
            if (cmd[1] == '-') return make_char(TF_SYMBOL, 0x2212, '-');
            return make_char(TF_SYMBOL, uint32_t(uint8_t(cmd[1])), cmd[1]);
        }

        /* Unknown control words are not prose.  Rendering their names as
         * letters silently invents mathematical content (`\\foo{x}` became
         * "foox").  Ignore the unsupported operator itself; any following
         * braced argument remains in the input stream and is still editable. */
        return nullptr;
    }

    /* A braced group whose characters all take one typeface. */
    NodePtr styled_group(int typeface) {
        NodeList inner;
        const bool preserveSpaces =
            typeface == TF_TEXT || typeface == TF_FUNCTION;
        if (preserveSpaces) {
            /* The delimiter after a control word is syntax, but whitespace
             * inside its braces is content.  Disable ordinary math-space
             * skipping only after the opening brace has been consumed. */
            skip_space();
            if (peek() == '{') {
                ++p_;
                ++literalSpaceDepth_;
                inner = parse_seq(kBrace);
                --literalSpaceDepth_;
                if (peek() == '}') ++p_;
            } else {
                inner = parse_arg();
            }
        } else {
            inner = parse_arg();
        }
        retype(inner, typeface);
        auto line = std::make_unique<LineNode>();
        line->children = std::move(inner);
        return line;
    }

    static void retype(NodeList& list, int typeface) {
        for (auto& n : list) {
            if (!n) continue;
            if (n->tag() == Node::kChar) {
                auto& ch = static_cast<CharNode&>(*n);
                /* Explicit TeX spaces retain their width/non-breaking
                 * semantics and split the surrounding text run. */
                if (ch.typeface != TF_SPACE) ch.typeface = typeface;
                ch.automaticFunction = false;
            } else {
                for (NodeList* slot : node_slots(*n))
                    if (slot) retype(*slot, typeface);
            }
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

    /* A separator that reads as a bra-ket bar rather than as a glyph in its
     * own right. */
    static bool is_braket_bar(const std::string& tok) {
        return tok == "|" || tok == "\\vert" || tok == "\\mid";
    }
    static bool is_angle(const std::string& tok) {
        return tok == "\\langle" || tok == "\\rangle";
    }

    NodePtr parse_fence() {
        eat_command("\\left");
        const std::string open = read_delimiter();
        auto f = std::make_unique<FenceNode>();

        /* Collect the \middle-separated segments of the body.  Almost every
         * fence has exactly one. */
        std::vector<NodeList> parts;
        std::vector<std::string> separators;
        parts.push_back(parse_seq(kFenceBody));
        while (peek_command() == "\\middle") {
            eat_command("\\middle");
            separators.push_back(read_delimiter());
            parts.push_back(parse_seq(kFenceBody));
        }

        std::string close = ".";
        if (eat_command("\\right")) close = read_delimiter();

        /* \left\langle a \middle| b \right\rangle is a bra-ket, and the tree
         * has a node for it whose bar stretches to the taller side. */
        if (parts.size() == 2 && is_angle(open) && is_angle(close) &&
            is_braket_bar(separators.front())) {
            auto dirac = std::make_unique<DiracNode>();
            dirac->bra = std::move(parts[0]);
            dirac->ket = std::move(parts[1]);
            return dirac;
        }

        /* Any other separator keeps its delimiter inline.  That loses the
         * stretching, but it is the glyph the author asked for. */
        f->content = std::move(parts[0]);
        for (size_t i = 0; i + 1 < parts.size(); ++i) {
            if (NodePtr sep = delimiter_char(separators[i]))
                f->content.push_back(std::move(sep));
            for (auto& n : parts[i + 1])
                if (n) f->content.push_back(std::move(n));
        }

        const bool hasOpen = open != ".";
        const bool hasClose = close != ".";
        const int openSelector = fence_selector(open);
        const int closeSelector = fence_selector(close);

        /* \left. x \right. is a sizing group with no visible delimiter.  It
         * has no fence to draw, so hand back the content rather than the
         * parentheses the old tmPAREN fallback invented. */
        if (!hasOpen && !hasClose) {
            auto line = std::make_unique<LineNode>();
            line->children = std::move(f->content);
            return line;
        }

        if (!hasOpen) {                     /* variation 2: right only */
            f->variation = 2;
            f->selector = closeSelector >= 0 ? closeSelector : tmPAREN;
            return f;
        }
        f->selector = openSelector >= 0 ? openSelector : tmPAREN;
        if (!hasClose) {                    /* variation 1: left only */
            f->variation = 1;
            return f;
        }
        /* Both sides visible.  A mismatched closing delimiter is meaningful --
         * \left( x \right] is a half-open interval -- so record it instead of
         * mirroring the opening one. */
        const int right = closeSelector >= 0 ? closeSelector : tmPAREN;
        if (right != f->selector) f->rightSelector = right;
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

        /* cases is a left brace around a table, and its second column is the
         * condition.  Flattening the row into one line silently merged
         * "x" and "x > 0", so it goes through the matrix path. */
        const bool is_cases = (env == "cases");
        if (is_cases) fence = tmBRACE;

        const bool is_equation =
            (env == "equation" || env == "equation*" ||
             env == "displaymath");
        const bool is_aligned =
            (env == "aligned" || env == "align" || env == "align*" ||
             env == "split" || env == "eqnarray" || env == "eqnarray*");
        const bool is_matrix =
            (env == "matrix" || fence >= 0 || is_aligned || is_cases);
        std::vector<std::vector<NodeList>> rows;
        rows.emplace_back();
        bool cellExpectedAfterAmpersand = false;

        for (;;) {
            skip_space();
            if (eof()) break;
            if (peek_command() == "\\end") {
                if (cellExpectedAfterAmpersand)
                    rows.back().emplace_back();
                eat_command("\\end"); read_env_name(); break;
            }

            rows.back().push_back(parse_seq(kCell));
            cellExpectedAfterAmpersand = false;

            skip_space();
            if (peek() == '&') {
                ++p_;
                cellExpectedAfterAmpersand = true;
                continue;
            }
            if (peek() == '\\' && peek(1) == '\\') { p_ += 2; rows.emplace_back(); continue; }
            if (peek_command() == "\\end") { eat_command("\\end"); read_env_name(); break; }
            if (eof()) break;
            ++p_;                                   /* never spin */
        }
        while (!rows.empty() && rows.back().empty()) rows.pop_back();

        size_t cols = 0;
        for (const auto& r : rows) cols = std::max(cols, r.size());

        if (is_equation) {
            auto line = std::make_unique<LineNode>();
            for (auto& row : rows)
                for (auto& cell : row)
                    for (auto& node : cell)
                        if (node) line->children.push_back(std::move(node));
            return line;
        }

        if (is_matrix) {
            auto m = std::make_unique<MatrixNode>();
            m->rows = int(rows.size());
            m->cols = int(cols);
            m->layoutKind = is_aligned ? MatrixNode::kAlignedLayout
                          : is_cases   ? MatrixNode::kCasesLayout
                                       : MatrixNode::kMatrixLayout;
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
            /* cases opens a brace and never closes it. */
            if (is_cases) f->variation = 1;
            f->content.push_back(std::move(m));
            return f;
        }

        /* gather/gather* and unknown environments -- one stacked line per
         * row.  Unknown wrappers are discarded rather than printed. */
        auto pile = std::make_unique<PileNode>();
        pile->ncols = int(cols);
        for (auto& r : rows) {
            auto line = std::make_unique<LineNode>();
            for (auto& cell : r)
                for (auto& n : cell) line->children.push_back(std::move(n));
            pile->lines.push_back(std::move(line));
        }
        return pile;
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

/* A bare TeX fragment such as `a\\\\b` is the natural source-pane spelling
 * of two equation rows.  At top level the parser previously treated the token
 * as a visible backslash command, even though the same token already means a
 * row break inside aligned/matrix/cases.  Detect only a genuinely top-level
 * separator: braces and environments keep their own row semantics. */
bool has_top_level_row_break(const std::string& source) {
    int braces = 0;
    int environments = 0;
    for (size_t i = 0; i < source.size();) {
        if (source.compare(i, 7, "\\begin{") == 0 ||
            source.compare(i, 5, "\\end{") == 0) {
            const bool opening = source.compare(i, 7, "\\begin{") == 0;
            const size_t tokenSize = opening ? 7 : 5;
            const size_t close = source.find('}', i + tokenSize);
            if (close == std::string::npos) break;
            if (opening)
                ++environments;
            else
                environments = std::max(0, environments - 1);
            i = close + 1;
            continue;
        }
        if (source[i] == '{') {
            ++braces;
            ++i;
            continue;
        }
        if (source[i] == '}') {
            braces = std::max(0, braces - 1);
            ++i;
            continue;
        }
        if (source[i] == '\\' && i + 1 < source.size() &&
            source[i + 1] == '\\') {
            if (braces == 0 && environments == 0) return true;
            i += 2;
            continue;
        }
        if (source[i] == '\\' && i + 1 < source.size()) {
            /* An escaped brace is data, not group structure.  Other commands
             * can be skipped one character at a time without hiding `\\\\`. */
            if (source[i + 1] == '{' || source[i + 1] == '}') {
                i += 2;
                continue;
            }
        }
        ++i;
    }
    return false;
}

}  // namespace

std::unique_ptr<LineNode> parse_latex(const std::string& latex,
                                      bool* depthExceeded) {
    std::string source = strip_delimiters(latex);
    if (has_top_level_row_break(source))
        source = "\\begin{aligned}" + source + "\\end{aligned}";
    TexParser p(std::move(source));
    std::unique_ptr<LineNode> root = p.parse();
    if (depthExceeded) *depthExceeded = p.depth_exceeded();
    return root;
}

}  // namespace eqnedit
