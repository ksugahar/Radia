/*
 * eq_edit.cpp -- the equation editing model
 *
 * The caret is a path of (child, slot) steps from the root plus an index in
 * the slot it lands in.  Every operation is expressed on that pair, which is
 * what lets arrow keys walk *into* a fraction the way Equation Editor's do
 * instead of skipping over it as one opaque object.
 */
#include "eq_edit.h"
#include "tex_parser.h"
#include "latex_emitter.h"
#include "tex2mtef.h"
#include "math_layout.h"
#include "mtef_common.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace mtef {

/* ------------------------------------------------------------------ slots */

std::vector<NodeList*> node_slots(Node& n) {
    switch (n.tag()) {
        case Node::kLine: {
            auto& l = static_cast<LineNode&>(n);
            return {&l.children};
        }
        case Node::kFrac: {
            auto& f = static_cast<FracNode&>(n);
            return {&f.numer, &f.denom};
        }
        case Node::kSqrt: {
            auto& s = static_cast<SqrtNode&>(n);
            if (s.hasIndex) return {&s.content, &s.index};
            return {&s.content};
        }
        case Node::kScript: {
            auto& s = static_cast<ScriptNode&>(n);
            std::vector<NodeList*> v{&s.base};
            if (s.hasSub) v.push_back(&s.sub);
            if (s.hasSup) v.push_back(&s.sup);
            return v;
        }
        case Node::kFence: {
            auto& f = static_cast<FenceNode&>(n);
            return {&f.content};
        }
        case Node::kIntegral: {
            auto& i = static_cast<IntegralNode&>(n);
            std::vector<NodeList*> v;
            if (i.hasLower) v.push_back(&i.lower);
            if (i.hasUpper) v.push_back(&i.upper);
            v.push_back(&i.body);
            return v;
        }
        case Node::kBigOp: {
            auto& b = static_cast<BigOpNode&>(n);
            std::vector<NodeList*> v;
            if (b.hasLower) v.push_back(&b.lower);
            if (b.hasUpper) v.push_back(&b.upper);
            v.push_back(&b.body);
            return v;
        }
        case Node::kEmbell: {
            auto& e = static_cast<EmbellNode&>(n);
            return {&e.content};
        }
        case Node::kDecoration: {
            auto& d = static_cast<DecorationNode&>(n);
            return {&d.content};
        }
        case Node::kBraceDeco: {
            auto& b = static_cast<BraceDecoNode&>(n);
            return {&b.content, &b.label};
        }
        case Node::kMatrix: {
            auto& m = static_cast<MatrixNode&>(n);
            std::vector<NodeList*> v;
            for (auto& e : m.elements)
                if (e && e->tag() == Node::kLine)
                    v.push_back(&static_cast<LineNode&>(*e).children);
            return v;
        }
        case Node::kPile: {
            auto& p = static_cast<PileNode&>(n);
            std::vector<NodeList*> v;
            for (auto& l : p.lines)
                if (l && l->tag() == Node::kLine)
                    v.push_back(&static_cast<LineNode&>(*l).children);
            return v;
        }
        default:
            return {};
    }
}

namespace {

std::string utf8_of(uint32_t cp) {
    std::string s;
    if (cp < 0x80) s += char(cp);
    else if (cp < 0x800) { s += char(0xC0 | (cp >> 6)); s += char(0x80 | (cp & 0x3F)); }
    else { s += char(0xE0 | (cp >> 12)); s += char(0x80 | ((cp >> 6) & 0x3F));
           s += char(0x80 | (cp & 0x3F)); }
    return s;
}

bool is_letter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }
bool is_digit(char c)  { return c >= '0' && c <= '9'; }

std::unique_ptr<CharNode> make_char(int tf, uint16_t code, char ascii = 0) {
    auto c = std::make_unique<CharNode>();
    c->typeface = tf; c->charCode = code; c->ch = ascii;
    return c;
}

int typeface_for_code(uint16_t code) {
    if (code >= 0x0391 && code <= 0x03A9) return TF_UCGREEK;
    if (code >= 0x03B1 && code <= 0x03C9) return TF_LCGREEK;
    return TF_SYMBOL;
}

/* Empty slots are what makes a template usable: you insert one and tab through
 * the holes.  The tree has no placeholder node, so "empty" is literally an
 * empty NodeList and the caret is put in the first one. */
bool first_empty_slot(Node& n, int& slot_out) {
    auto slots = node_slots(n);
    for (size_t i = 0; i < slots.size(); ++i)
        if (slots[i]->empty()) { slot_out = int(i); return true; }
    return false;
}

}  // namespace

/* ------------------------------------------------------------- lifecycle */

Equation::Equation() : root_(std::make_unique<LineNode>()) {}

void Equation::load_latex(const std::string& latex) {
    root_ = parse_latex(latex);
    if (!root_) root_ = std::make_unique<LineNode>();
    path_.clear();
    index_ = int(root_->children.size());
    undo_.clear();
    redo_.clear();
}

std::string Equation::latex() const { return tree_to_latex(*root_); }

std::string Equation::omml(const OmmlOptions& opt) const {
    return render_omml(*root_, opt, /*run_passes=*/false);
}

std::string Equation::svg(const SvgStyle& style) const {
    return render_svg(*root_, style);
}

/* ----------------------------------------------------------- addressing */

NodeList* Equation::slot_at(const std::vector<CaretStep>& path) const {
    NodeList* list = const_cast<NodeList*>(&root_->children);
    for (const CaretStep& st : path) {
        if (st.child < 0 || size_t(st.child) >= list->size()) return nullptr;
        Node* n = (*list)[st.child].get();
        if (!n) return nullptr;
        auto slots = node_slots(*n);
        if (st.slot < 0 || size_t(st.slot) >= slots.size()) return nullptr;
        list = slots[st.slot];
    }
    return list;
}

NodeList& Equation::here() {
    NodeList* l = slot_at(path_);
    if (!l) { path_.clear(); l = &root_->children; }
    return *l;
}

Node* Equation::parent_node(const std::vector<CaretStep>& path) const {
    if (path.empty()) return nullptr;
    std::vector<CaretStep> up(path.begin(), path.end() - 1);
    NodeList* l = slot_at(up);
    if (!l) return nullptr;
    int child = path.back().child;
    if (child < 0 || size_t(child) >= l->size()) return nullptr;
    return (*l)[child].get();
}

void Equation::clamp() {
    NodeList& l = here();
    if (index_ < 0) index_ = 0;
    if (size_t(index_) > l.size()) index_ = int(l.size());
}

/* --------------------------------------------------------------- caret */

std::string Equation::caret_text() const {
    std::string s;
    for (const CaretStep& st : path_) {
        if (!s.empty()) s += '/';
        s += std::to_string(st.child) + "." + std::to_string(st.slot);
    }
    s += ":" + std::to_string(index_);
    return s;
}

void Equation::set_caret_text(const std::string& s) {
    path_.clear();
    index_ = 0;
    size_t colon = s.rfind(':');
    std::string p = (colon == std::string::npos) ? s : s.substr(0, colon);
    if (colon != std::string::npos) index_ = atoi(s.c_str() + colon + 1);
    size_t i = 0;
    while (i < p.size()) {
        size_t slash = p.find('/', i);
        std::string step = p.substr(i, slash == std::string::npos ? std::string::npos
                                                                 : slash - i);
        size_t dot = step.find('.');
        if (dot != std::string::npos)
            path_.push_back({atoi(step.c_str()), atoi(step.c_str() + dot + 1)});
        if (slash == std::string::npos) break;
        i = slash + 1;
    }
    if (!slot_at(path_)) path_.clear();
    clamp();
}

std::string Equation::caret() const { return caret_text(); }

void Equation::move_home() { index_ = 0; }
void Equation::move_end()  { index_ = int(here().size()); }

bool Equation::move_out() {
    if (path_.empty()) return false;
    int child = path_.back().child;
    path_.pop_back();
    index_ = child + 1;
    clamp();
    return true;
}

bool Equation::move_left() {
    NodeList& l = here();
    if (index_ > 0) {
        Node* prev = l[index_ - 1].get();
        auto slots = prev ? node_slots(*prev) : std::vector<NodeList*>{};
        if (!slots.empty()) {                      /* step into its last slot */
            path_.push_back({index_ - 1, int(slots.size()) - 1});
            index_ = int(slots.back()->size());
            return true;
        }
        --index_;
        return true;
    }
    if (path_.empty()) return false;
    /* At the start of a slot: try the previous slot of the same node. */
    CaretStep st = path_.back();
    if (st.slot > 0) {
        path_.back().slot = st.slot - 1;
        NodeList* l2 = slot_at(path_);
        index_ = l2 ? int(l2->size()) : 0;
        return true;
    }
    path_.pop_back();
    index_ = st.child;
    clamp();
    return true;
}

bool Equation::move_right() {
    NodeList& l = here();
    if (size_t(index_) < l.size()) {
        Node* next = l[index_].get();
        auto slots = next ? node_slots(*next) : std::vector<NodeList*>{};
        if (!slots.empty()) {
            path_.push_back({index_, 0});
            index_ = 0;
            return true;
        }
        ++index_;
        return true;
    }
    if (path_.empty()) return false;
    CaretStep st = path_.back();
    Node* parent = parent_node(path_);
    int nslots = parent ? int(node_slots(*parent).size()) : 0;
    if (st.slot + 1 < nslots) {
        path_.back().slot = st.slot + 1;
        index_ = 0;
        return true;
    }
    path_.pop_back();
    index_ = st.child + 1;
    clamp();
    return true;
}

bool Equation::next_slot() {
    if (path_.empty()) return false;
    for (;;) {
        CaretStep st = path_.back();
        Node* parent = parent_node(path_);
        int nslots = parent ? int(node_slots(*parent).size()) : 0;
        if (st.slot + 1 < nslots) {
            path_.back().slot = st.slot + 1;
            index_ = 0;
            return true;
        }
        path_.pop_back();
        index_ = st.child + 1;
        if (path_.empty()) { clamp(); return true; }
    }
}

bool Equation::prev_slot() {
    if (path_.empty()) return false;
    CaretStep st = path_.back();
    if (st.slot > 0) {
        path_.back().slot = st.slot - 1;
        index_ = 0;
        return true;
    }
    path_.pop_back();
    index_ = st.child;
    clamp();
    return true;
}

/* ------------------------------------------------------------ geometry */

Equation::CaretGeometry Equation::caret_geometry(const SvgStyle& style) const {
    CaretGeometry g;
    Layout L = layout_math(*root_, style);
    /* caret_text() is "path:index"; the layout keys stops on the path alone. */
    std::string ct = caret_text();
    size_t colon = ct.rfind(':');
    std::string path = (colon == std::string::npos) ? ct : ct.substr(0, colon);
    int index = (colon == std::string::npos) ? 0 : atoi(ct.c_str() + colon + 1);

    const CaretStop* s = find_stop(L, path, index);
    if (!s) return g;
    g.found = true;
    g.x = s->x;
    g.top = s->top;
    g.bottom = s->bottom;
    return g;
}

bool Equation::move_to_point(double x, double y, const SvgStyle& style) {
    Layout L = layout_math(*root_, style);
    const CaretStop* s = nearest_stop(L, x, y);
    if (!s) return false;
    set_caret_text(s->path + ":" + std::to_string(s->index));
    return true;
}

void Equation::extents(double& w, double& asc, double& desc,
                       const SvgStyle& style) const {
    Layout L = layout_math(*root_, style);
    w = L.w;
    asc = L.asc;
    desc = L.desc;
}

Layout Equation::layout(const SvgStyle& style) const {
    return layout_math(*root_, style);
}

/* ------------------------------------------------------------- history */

void Equation::checkpoint() {
    undo_.push_back({latex(), caret_text()});
    redo_.clear();
    if (undo_.size() > 200) undo_.erase(undo_.begin());
}

void Equation::restore(const Snapshot& s) {
    root_ = parse_latex(s.latex);
    if (!root_) root_ = std::make_unique<LineNode>();
    set_caret_text(s.caret);
}

bool Equation::undo() {
    if (undo_.empty()) return false;
    redo_.push_back({latex(), caret_text()});
    Snapshot s = undo_.back();
    undo_.pop_back();
    restore(s);
    return true;
}

bool Equation::redo() {
    if (redo_.empty()) return false;
    undo_.push_back({latex(), caret_text()});
    Snapshot s = redo_.back();
    redo_.pop_back();
    restore(s);
    return true;
}

/* ------------------------------------------------------------- editing */

void Equation::insert_text(const std::string& utf8) {
    checkpoint();
    NodeList& l = here();
    clamp();
    /* One character, one node -- and a character may be several bytes.  What
     * an IME hands the window is a composed character, so walking bytes here
     * would turn every Japanese keystroke into three nodes the same way the
     * parser once did. */
    for (size_t i = 0; i < utf8.size(); ) {
        const char c = utf8[i];
        std::unique_ptr<CharNode> node;
        if ((unsigned char)c >= 0x80) {
            const int extra = ((unsigned char)c >= 0xF0) ? 3 :
                              ((unsigned char)c >= 0xE0) ? 2 : 1;
            uint32_t cp = uint32_t((unsigned char)c) & uint32_t(0x3F >> extra);
            ++i;
            for (int k = 0; k < extra && i < utf8.size(); ++k, ++i)
                cp = (cp << 6) | (uint32_t((unsigned char)utf8[i]) & 0x3Fu);
            node = make_char(TF_SYMBOL, uint16_t(cp > 0xFFFF ? 0xFFFD : cp));
        } else {
            ++i;
            if (is_digit(c))       node = make_char(TF_NUMBER,   uint16_t(c), c);
            else if (is_letter(c)) node = make_char(TF_VARIABLE, uint16_t(c), c);
            else if (c == '-')     node = make_char(TF_SYMBOL,   0x2212, c);
            else                   node = make_char(TF_SYMBOL,
                                                    uint16_t((unsigned char)c), c);
        }
        l.insert(l.begin() + index_, std::move(node));
        ++index_;
    }
}

bool Equation::insert_symbol(const std::string& cmd) {
    int code = tex_command_to_unicode(cmd.c_str());
    if (code < 0) return false;
    checkpoint();
    NodeList& l = here();
    clamp();
    auto ch = make_char(typeface_for_code(uint16_t(code)), uint16_t(code));
    ch->latex = cmd;
    l.insert(l.begin() + index_, std::move(ch));
    ++index_;
    return true;
}

/* Descend to the first hole a person would start typing in.  It is not always
 * a direct slot: `cases` is a fence wrapping a pile, and its first hole is a
 * row of that pile, two levels down. */
void Equation::enter_first_empty_slot(Node& n) {
    Node* cur = &n;
    int child = index_;
    for (;;) {
        int slot = 0;
        auto slots = node_slots(*cur);
        if (slots.empty()) return;
        if (!first_empty_slot(*cur, slot)) {
            /* No empty slot here -- follow the first slot that holds exactly
             * one node and look inside that. */
            slot = 0;
            path_.push_back({child, slot});
            index_ = 0;
            if (slots[0]->size() != 1) return;
            cur = (*slots[0])[0].get();
            child = 0;
            if (!cur) return;
            continue;
        }
        path_.push_back({child, slot});
        index_ = 0;
        return;
    }
}

bool Equation::insert_template(const std::string& kind) {
    NodePtr node;

    auto fence = [](int sel) {
        auto f = std::make_unique<FenceNode>();
        f->selector = sel;
        return f;
    };

    if (kind == "frac" || kind == "slashfrac") {
        auto f = std::make_unique<FracNode>();
        f->slashed = (kind == "slashfrac");
        node = std::move(f);
    } else if (kind == "sqrt") {
        node = std::make_unique<SqrtNode>();
    } else if (kind == "nthroot") {
        auto s = std::make_unique<SqrtNode>();
        s->hasIndex = true;
        node = std::move(s);
    } else if (kind == "sub" || kind == "sup" || kind == "subsup") {
        auto s = std::make_unique<ScriptNode>();
        s->hasSub = (kind != "sup");
        s->hasSup = (kind != "sub");
        /* The character just typed is the natural base, exactly as pressing
         * Ctrl+L after `x` gives `x` a subscript rather than an empty box. */
        NodeList& l = here();
        if (index_ > 0) {
            s->base.push_back(std::move(l[index_ - 1]));
            l.erase(l.begin() + index_ - 1);
            --index_;
        }
        node = std::move(s);
    } else if (kind == "paren")   { node = fence(tmPAREN); }
    else if (kind == "bracket")   { node = fence(tmBRACK); }
    else if (kind == "brace")     { node = fence(tmBRACE); }
    else if (kind == "abs")       { node = fence(tmBAR); }
    else if (kind == "angle")     { node = fence(tmANGLE); }
    else if (kind == "int" || kind == "oint" || kind == "iint") {
        auto i = std::make_unique<IntegralNode>();
        i->selector = (kind == "oint") ? tmSSINT
                    : (kind == "iint") ? tmDINT : tmSINT;
        i->hasLower = i->hasUpper = true;
        node = std::move(i);
    } else if (kind == "sum" || kind == "prod") {
        auto b = std::make_unique<BigOpNode>();
        b->selector = (kind == "sum") ? tmSUM : tmPROD;
        b->hasLower = b->hasUpper = true;
        b->hasLimits = true;
        node = std::move(b);
    } else if (kind == "overline" || kind == "underline") {
        auto d = std::make_unique<DecorationNode>();
        d->selector = (kind == "overline") ? tmOBAR : tmUBAR;
        node = std::move(d);
    } else if (kind == "hat" || kind == "vec" || kind == "bar" ||
               kind == "tilde" || kind == "dot") {
        auto e = std::make_unique<EmbellNode>();
        e->embellType = (kind == "hat")   ? 9
                      : (kind == "vec")   ? 11
                      : (kind == "bar")   ? 17
                      : (kind == "tilde") ? 8 : 2;
        node = std::move(e);
    } else if (kind == "matrix2x2" || kind == "matrix3x3") {
        const int n = (kind == "matrix2x2") ? 2 : 3;
        auto m = std::make_unique<MatrixNode>();
        m->rows = m->cols = n;
        for (int i = 0; i < n * n; ++i)
            m->elements.push_back(std::make_unique<LineNode>());
        node = std::move(m);
    } else if (kind == "cases") {
        auto p = std::make_unique<PileNode>();
        p->lines.push_back(std::make_unique<LineNode>());
        p->lines.push_back(std::make_unique<LineNode>());
        auto f = fence(tmBRACE);
        f->variation = 1;
        f->content.push_back(std::move(p));
        node = std::move(f);
    } else {
        return false;
    }

    checkpoint();
    NodeList& l = here();
    clamp();
    Node* raw = node.get();
    l.insert(l.begin() + index_, std::move(node));
    enter_first_empty_slot(*raw);
    return true;
}

bool Equation::backspace() {
    NodeList& l = here();
    clamp();
    if (index_ > 0) {
        checkpoint();
        NodeList& cur = here();
        cur.erase(cur.begin() + index_ - 1);
        --index_;
        return true;
    }
    /* At the start of an empty slot, backspace unwraps the template: its
     * contents are spliced into the parent, so nothing is silently lost. */
    if (path_.empty()) return false;
    CaretStep st = path_.back();
    std::vector<CaretStep> up(path_.begin(), path_.end() - 1);
    NodeList* parent_list = slot_at(up);
    if (!parent_list || size_t(st.child) >= parent_list->size()) return false;

    checkpoint();
    parent_list = slot_at(up);
    NodePtr owned = std::move((*parent_list)[st.child]);
    parent_list->erase(parent_list->begin() + st.child);

    NodeList salvage;
    for (NodeList* s : node_slots(*owned))
        for (auto& n : *s)
            if (n) salvage.push_back(std::move(n));

    int at = st.child;
    for (auto& n : salvage)
        parent_list->insert(parent_list->begin() + at++, std::move(n));

    path_ = up;
    index_ = st.child;
    clamp();
    return true;
}

bool Equation::erase() {
    NodeList& l = here();
    clamp();
    if (size_t(index_) >= l.size()) return false;
    checkpoint();
    NodeList& cur = here();
    cur.erase(cur.begin() + index_);
    return true;
}

/* ------------------------------------------------------------ commands */

const std::vector<std::string>& Equation::templates() {
    /* "slashfrac" is deliberately absent.  The LaTeX emitter writes a slashed
     * fraction as MTEF's {}^{a}/{}_{b}, which does not read back as one
     * fraction -- so offering it in the editor would let an equation change
     * shape when it is saved and reopened.  FracNode::slashed still exists for
     * reading legacy MTEF; it is only unavailable as a template. */
    static const std::vector<std::string> kAll = {
        "frac", "sqrt", "nthroot",
        "sub", "sup", "subsup",
        "paren", "bracket", "brace", "abs", "angle",
        "int", "iint", "oint", "sum", "prod",
        "overline", "underline", "hat", "vec", "bar", "tilde", "dot",
        "matrix2x2", "matrix3x3", "cases",
    };
    return kAll;
}

/* Equation Editor 3.0's chords.  Only the ones that are unambiguous in the
 * editor's own Help are listed; everything else in templates() is reachable
 * from a palette and is deliberately left unbound rather than invented. */
/* Which palette a symbol belongs to, decided from its code point.  Equation
 * Editor 3.0's arrangement, reimplemented -- the grouping is observable from
 * using it, which is all we take. */
static int palette_of(uint32_t c) {
    enum { kRelation, kEllipsis, kOperator, kArrow, kLogic, kSet, kMisc,
           kGreekLower, kGreekUpper };
    if (c >= 0x03B1 && c <= 0x03C9) return kGreekLower;
    if (c >= 0x0391 && c <= 0x03A9) return kGreekUpper;
    if ((c >= 0x2190 && c <= 0x21FF) || (c >= 0x27F8 && c <= 0x27FA))
        return kArrow;

    switch (c) {
        case 0x2026: case 0x22EF: case 0x22EE: case 0x22F1: case 0x22F0:
            return kEllipsis;
        default:
            break;
    }

    switch (c) {
        case 0x2264: case 0x2265: case 0x2260: case 0x2248: case 0x2261:
        case 0x221D: case 0x223C: case 0x2245: case 0x226A: case 0x226B:
        case 0x2250: case 0x2243:
            return kRelation;

        case 0x00B1: case 0x2213: case 0x00D7: case 0x00F7: case 0x2217:
        case 0x2218: case 0x2219: case 0x2295: case 0x2297: case 0x2299:
        case 0x22C5: case 0x2216:
            return kOperator;

        case 0x2200: case 0x2203: case 0x2234: case 0x2235: case 0x00AC:
        case 0x2227: case 0x2228: case 0x22A2: case 0x22A8:
            return kLogic;

        case 0x2208: case 0x2209: case 0x220B: case 0x2282: case 0x2283:
        case 0x2286: case 0x2287: case 0x222A: case 0x2229: case 0x2205:
        case 0x2284:
            return kSet;

        default:
            break;
    }
    /* Miscellaneous is the catch-all, and it is TOTAL: every entry in the
     * table lands somewhere.  Classifying only the ranges I thought of left
     * \ldots, \cdots, \dagger and \bullet insertable from the keyboard but
     * absent from the bar -- the palette was quietly hiding part of the
     * editor, which is the one thing a generated palette exists to prevent. */
    return kMisc;
}

const std::vector<Equation::PaletteGroup>& Equation::symbol_palettes() {
    static const std::vector<PaletteGroup> kGroups = [] {
        /* This array is indexed by palette_of()'s enum.  Keep the two in step:
         * adding a group to one and not the other silently relabels every
         * palette after it, which is how "Greek" came to hold arrows. */
        static const char* kNames[] = {
            "Relations", "Ellipses", "Operators", "Arrows", "Logic",
            "Set theory", "Miscellaneous", "Greek", "Greek capitals",
        };
        const int n = int(sizeof(kNames) / sizeof(kNames[0]));
        std::vector<PaletteGroup> groups;
        for (int i = 0; i < n; ++i) groups.push_back({kNames[i], {}});

        for (int i = 0; i < tex_command_count(); ++i) {
            const int code = tex_command_code_at(i);
            if (code < 0) continue;
            const int g = palette_of(uint32_t(code));
            if (g < 0 || g >= n) continue;
            PaletteItem it;
            it.command = tex_command_name(i);
            it.code = uint32_t(code);
            groups[size_t(g)].items.push_back(it);
        }

        /* Embellishments are templates, not table symbols, but Equation
         * Editor keeps them on the symbol row and so do we. */
        PaletteGroup emb{"Embellishments", {}};
        for (const char* t : {"hat", "vec", "bar", "tilde", "dot",
                              "overline", "underline"}) {
            PaletteItem it;
            it.command = t;
            it.is_template = true;
            emb.items.push_back(it);
        }
        groups.insert(groups.begin() + 2, emb);

        groups.erase(std::remove_if(groups.begin(), groups.end(),
                         [](const PaletteGroup& g) { return g.items.empty(); }),
                     groups.end());
        return groups;
    }();
    return kGroups;
}

const std::vector<Equation::PaletteGroup>& Equation::template_palettes() {
    static const std::vector<PaletteGroup> kGroups = [] {
        struct Def { const char* name; std::vector<const char*> items; };
        static const Def kDefs[] = {
            {"Fences",        {"paren", "bracket", "brace", "abs", "angle"}},
            {"Fractions",     {"frac", "sqrt", "nthroot"}},
            {"Scripts",       {"sub", "sup", "subsup"}},
            {"Summation",     {"sum"}},
            {"Integrals",     {"int", "iint", "oint"}},
            {"Products",      {"prod"}},
            {"Matrices",      {"matrix2x2", "matrix3x3", "cases"}},
        };
        /* Equation Editor also has a labelled-arrow palette.  This build has
         * no such template, so the button is absent rather than empty. */
        const std::vector<std::string>& known = templates();
        std::vector<PaletteGroup> groups;
        for (const Def& d : kDefs) {
            PaletteGroup g{d.name, {}};
            for (const char* t : d.items) {
                if (std::find(known.begin(), known.end(), t) == known.end())
                    continue;          /* never offer what cannot be inserted */
                PaletteItem it;
                it.command = t;
                it.is_template = true;
                g.items.push_back(it);
            }
            if (!g.items.empty()) groups.push_back(g);
        }
        return groups;
    }();
    return kGroups;
}

const std::vector<Equation::Binding>& Equation::shortcuts() {
    static const std::vector<Binding> kBindings = {
        {"Ctrl+F",       "template.frac",      "fraction"},
        {"Ctrl+R",       "template.sqrt",      "square root"},
        {"Ctrl+H",       "template.sup",       "superscript (high)"},
        {"Ctrl+L",       "template.sub",       "subscript (low)"},
        {"Ctrl+J",       "template.subsup",    "subscript and superscript"},
        {"Ctrl+9",       "template.paren",     "parentheses"},
        {"Ctrl+[",       "template.bracket",   "brackets"},
        {"Ctrl+{",       "template.brace",     "braces"},
        {"Ctrl+I",       "template.int",       "integral"},
        {"Ctrl+T, S",    "template.sum",       "summation"},
        {"Ctrl+T, P",    "template.prod",      "product"},
        {"Ctrl+T, M",    "template.matrix2x2", "matrix"},
        {"Tab",          "caret.next_slot",    "next slot"},
        {"Shift+Tab",    "caret.prev_slot",    "previous slot"},
        {"Left",         "caret.left",         "left"},
        {"Right",        "caret.right",        "right"},
        {"Home",         "caret.home",         "start of slot"},
        {"End",          "caret.end",          "end of slot"},
        {"Backspace",    "edit.backspace",     "delete backwards / unwrap"},
        {"Delete",       "edit.delete",        "delete forwards"},
        {"Ctrl+Z",       "edit.undo",          "undo"},
        {"Ctrl+Y",       "edit.redo",          "redo"},
    };
    return kBindings;
}

bool Equation::command(const std::string& name) {
    if (name.compare(0, 9, "template.") == 0)
        return insert_template(name.substr(9));
    if (name == "caret.left")       return move_left();
    if (name == "caret.right")      return move_right();
    if (name == "caret.next_slot")  return next_slot();
    if (name == "caret.prev_slot")  return prev_slot();
    if (name == "caret.out")        return move_out();
    if (name == "caret.home")       { move_home(); return true; }
    if (name == "caret.end")        { move_end();  return true; }
    if (name == "edit.backspace")   return backspace();
    if (name == "edit.delete")      return erase();
    if (name == "edit.undo")        return undo();
    if (name == "edit.redo")        return redo();
    return false;
}

}  // namespace mtef
