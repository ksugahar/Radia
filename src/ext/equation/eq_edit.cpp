/*
 * eq_edit.cpp -- the equation editing model
 *
 * The caret is a path of (child, slot) steps from the root plus an index in
 * the slot it lands in.  Every operation is expressed on that pair, which is
 * what lets arrow keys walk *into* a fraction the way Equation Editor's do
 * instead of skipping over it as one opaque object.
 */
#include <cmath>
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

std::string utf8_of(uint32_t cp) { return mtef_utf8_of(cp); }

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
    anchor_ = -1;
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

/* A plain move drops the selection.  Shift is what keeps it, and a caret that
 * wandered off while a range stayed highlighted would delete something other
 * than what is shown. */
void Equation::move_home() { clear_selection(); index_ = 0; }
void Equation::move_end()  { clear_selection(); index_ = int(here().size()); }

bool Equation::move_out() {
    clear_selection();
    if (path_.empty()) return false;
    int child = path_.back().child;
    path_.pop_back();
    index_ = child + 1;
    clamp();
    return true;
}

/* One whole item left or right, without stepping inside it.
 *
 * move_left/move_right descend into a template because that is how you get at
 * the letters in it.  This is the other motion: past the template.  At the
 * edge of a slot it does NOT jump out -- the same reasoning as the selection,
 * that a motion which silently changes which slot you are in is worse than one
 * that stops. */
/* Enter: break the line here.
 *
 * The slot the caret is in becomes two lines of a stack, split at the caret --
 * or, if it already IS one line of a stack, a new line is opened after it.
 * That is the operation Equation Editor binds to Enter, and it is only
 * possible now that a pile has a layout: before, a second line would have been
 * stored and drawn as nothing.
 *
 * It splits rather than merely appending because that is what Enter does
 * everywhere else -- pressing it in the middle of a line takes the rest of the
 * line with it. */
bool Equation::newline() {
    checkpoint();
    delete_selection();
    clear_selection();

    Node* parent = parent_node(path_);
    if (parent && parent->tag() == Node::kPile) {
        auto& pile = static_cast<PileNode&>(*parent);
        const int which = path_.back().slot;
        NodeList& cur = here();
        auto fresh = std::make_unique<LineNode>();
        for (size_t i = size_t(index_); i < cur.size(); ++i)
            fresh->children.push_back(std::move(cur[i]));
        cur.resize(size_t(index_));
        pile.lines.insert(pile.lines.begin() + which + 1, std::move(fresh));
        path_.back().slot = which + 1;
        index_ = 0;
        return true;
    }

    /* Not in a stack yet: wrap what is here in one. */
    NodeList& cur = here();
    auto first = std::make_unique<LineNode>();
    auto second = std::make_unique<LineNode>();
    for (size_t i = 0; i < cur.size(); ++i) {
        if (int(i) < index_) first->children.push_back(std::move(cur[i]));
        else                 second->children.push_back(std::move(cur[i]));
    }
    cur.clear();
    auto pile = std::make_unique<PileNode>();
    pile->halign = 0;               /* gathered, which is what Enter means */
    pile->ncols = 1;
    pile->lines.push_back(std::move(first));
    pile->lines.push_back(std::move(second));
    cur.push_back(std::move(pile));
    path_.push_back({0, 1});        /* into the pile, its second line */
    index_ = 0;
    return true;
}

bool Equation::move_item(int dir, bool extend) {
    if (extend) {
        if (anchor_ < 0) anchor_ = index_;
    } else {
        clear_selection();
    }
    NodeList& l = here();
    const int want = index_ + (dir < 0 ? -1 : +1);
    if (want < 0 || want > int(l.size())) return false;
    index_ = want;
    if (extend && anchor_ == index_) clear_selection();
    return true;
}

bool Equation::move_left() {
    clear_selection();
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
    clear_selection();
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
    clear_selection();
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
    clear_selection();
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

/* --------------------------------------------------------------- style */

/* The Adobe Symbol keyboard, which is the one Equation Editor's Greek style
 * types on: `a` gives alpha, `q` gives theta.  Greek is not a font here -- the
 * character itself changes -- because restyling a Latin `a` would leave it a
 * Latin `a` in the OMML, and Word would set it in a Latin font. */
static uint32_t greek_of(uint32_t c) {
    switch (c) {
        case 'a': return 0x3B1; case 'b': return 0x3B2; case 'c': return 0x3C7;
        case 'd': return 0x3B4; case 'e': return 0x3B5; case 'f': return 0x3C6;
        case 'g': return 0x3B3; case 'h': return 0x3B7; case 'i': return 0x3B9;
        case 'j': return 0x3D5; case 'k': return 0x3BA; case 'l': return 0x3BB;
        case 'm': return 0x3BC; case 'n': return 0x3BD; case 'o': return 0x3BF;
        case 'p': return 0x3C0; case 'q': return 0x3B8; case 'r': return 0x3C1;
        case 's': return 0x3C3; case 't': return 0x3C4; case 'u': return 0x3C5;
        case 'v': return 0x3D6; case 'w': return 0x3C9; case 'x': return 0x3BE;
        case 'y': return 0x3C8; case 'z': return 0x3B6;
        case 'A': return 0x391; case 'B': return 0x392; case 'C': return 0x3A7;
        case 'D': return 0x394; case 'E': return 0x395; case 'F': return 0x3A6;
        case 'G': return 0x393; case 'H': return 0x397; case 'I': return 0x399;
        case 'J': return 0x3D1; case 'K': return 0x39A; case 'L': return 0x39B;
        case 'M': return 0x39C; case 'N': return 0x39D; case 'O': return 0x39F;
        case 'P': return 0x3A0; case 'Q': return 0x398; case 'R': return 0x3A1;
        case 'S': return 0x3A3; case 'T': return 0x3A4; case 'U': return 0x3A5;
        case 'V': return 0x3C2; case 'W': return 0x3A9; case 'X': return 0x39E;
        case 'Y': return 0x3A8; case 'Z': return 0x396;
        default:  return 0;
    }
}

/* Walk everything, not just the immediate children: a selection may hold a
 * whole fraction, and restyling only its top level would leave the numerator
 * in the old face. */
static void apply_style(NodeList& list, const std::string& name);

static void style_one(Node& n, const std::string& name) {
    if (n.tag() == Node::kChar) {
        auto& c = static_cast<CharNode&>(n);
        uint32_t cp = c.charCode ? c.charCode : uint32_t(uint8_t(c.ch));

        if (name == "greek") {
            if (uint32_t g = greek_of(cp)) {
                c.charCode = uint16_t(g);
                c.ch = 0;
                c.latex.clear();
                c.typeface = (g >= 0x3B1) ? TF_LCGREEK : TF_UCGREEK;
            }
            return;
        }
        if (name == "text")     { c.typeface = TF_TEXT;     return; }
        if (name == "function") { c.typeface = TF_FUNCTION; return; }
        if (name == "variable") { c.typeface = TF_VARIABLE; return; }
        if (name == "vector")   { c.typeface = TF_VECTOR;   return; }
        if (name == "math") {
            /* Back to what the parser would have chosen: a letter is a
             * variable, a digit is a number, anything else is a symbol. */
            if (cp >= '0' && cp <= '9')                          c.typeface = TF_NUMBER;
            else if ((cp >= 'a' && cp <= 'z') ||
                     (cp >= 'A' && cp <= 'Z'))                   c.typeface = TF_VARIABLE;
            else if (cp >= 0x3B1 && cp <= 0x3C9)                 c.typeface = TF_LCGREEK;
            else if (cp >= 0x391 && cp <= 0x3A9)                 c.typeface = TF_UCGREEK;
            else                                                 c.typeface = TF_SYMBOL;
            return;
        }
        return;
    }
    for (NodeList* slot : node_slots(n))
        if (slot) apply_style(*slot, name);
}

static void apply_style(NodeList& list, const std::string& name) {
    for (auto& n : list) if (n) style_one(*n, name);
}

const std::vector<std::string>& Equation::styles() {
    /* Equation Editor's Style menu, minus Other... and Define..., which are
     * dialogs over a font list rather than styles in their own right. */
    static const std::vector<std::string> kAll = {
        "math", "text", "function", "variable", "vector", "greek",
    };
    return kAll;
}

bool Equation::set_style(const std::string& name) {
    const std::vector<std::string>& all = styles();
    if (std::find(all.begin(), all.end(), name) == all.end()) return false;

    /* The mode changes either way.  Requiring a selection made one Greek
     * letter cost three operations -- type it, select it, restyle it -- where
     * the editor being imitated costs one, because there the style is
     * something you are IN. */
    style_ = name;
    if (!has_selection()) return true;

    checkpoint();
    NodeList& l = here();
    const int lo = std::min(anchor_, index_);
    const int hi = std::max(anchor_, index_);
    for (int i = lo; i < hi; ++i)
        if (l[size_t(i)]) style_one(*l[size_t(i)], name);
    return true;
}

/* ----------------------------------------------------------- selection */

bool Equation::has_selection() const {
    return anchor_ >= 0 && anchor_ != index_;
}

void Equation::clear_selection() { anchor_ = -1; }

void Equation::select_all() {
    anchor_ = 0;
    index_ = int(here().size());
}

/* Shift+move starts a selection at wherever the caret was. */
bool Equation::extend_left() {
    if (anchor_ < 0) anchor_ = index_;
    if (index_ <= 0) return false;      /* stop at the slot edge */
    --index_;
    return true;
}

bool Equation::extend_right() {
    if (anchor_ < 0) anchor_ = index_;
    if (size_t(index_) >= here().size()) return false;
    ++index_;
    return true;
}

void Equation::extend_home() {
    if (anchor_ < 0) anchor_ = index_;
    index_ = 0;
}

void Equation::extend_end() {
    if (anchor_ < 0) anchor_ = index_;
    index_ = int(here().size());
}

bool Equation::extend_to_point(double x, double y, const SvgStyle& style) {
    const int keep = (anchor_ < 0) ? index_ : anchor_;
    const std::vector<CaretStep> keep_path = path_;
    if (!move_to_point(x, y, style)) return false;
    /* Dragging out of the slot the selection started in would silently change
     * what is anchored, so the anchor holds and the caret is clamped back. */
    if (path_ != keep_path) {
        path_ = keep_path;
        clamp();
        index_ = (index_ < keep) ? 0 : int(here().size());
    }
    anchor_ = keep;
    return true;
}

std::string Equation::selected_latex() {
    if (!has_selection()) return std::string();
    NodeList& l = here();
    const int lo = std::min(anchor_, index_);
    const int hi = std::max(anchor_, index_);

    /* Borrow the nodes into a line, write it, and put them back.  Teaching
     * every node type to clone itself would be a lot of code for one caller,
     * and a clone that fell behind a new node type would silently drop it. */
    LineNode tmp;
    for (int i = lo; i < hi; ++i) tmp.children.push_back(std::move(l[size_t(i)]));
    std::string out = tree_to_latex(tmp);
    for (int i = lo; i < hi; ++i)
        l[size_t(i)] = std::move(tmp.children[size_t(i - lo)]);
    return out;
}

NodeList Equation::take_selection() {
    NodeList taken;
    if (!has_selection()) return taken;
    NodeList& l = here();
    const int lo = std::min(anchor_, index_);
    const int hi = std::max(anchor_, index_);
    for (int i = lo; i < hi; ++i) taken.push_back(std::move(l[size_t(i)]));
    l.erase(l.begin() + lo, l.begin() + hi);
    index_ = lo;
    anchor_ = -1;
    clamp();
    return taken;
}

bool Equation::delete_selection() {
    if (!has_selection()) return false;
    checkpoint();
    take_selection();
    return true;
}

Equation::SelectionBox Equation::selection_geometry(const SvgStyle& style) const {
    SelectionBox box;
    if (!has_selection()) return box;

    Layout L = layout_math(*root_, style);
    std::string ct = caret_text();
    size_t colon = ct.rfind(':');
    std::string path = (colon == std::string::npos) ? ct : ct.substr(0, colon);

    const CaretStop* a = find_stop(L, path, std::min(anchor_, index_));
    const CaretStop* b = find_stop(L, path, std::max(anchor_, index_));
    if (!a || !b) return box;

    box.found = true;
    box.x0 = a->x;
    box.x1 = b->x;
    box.top = std::min(a->top, b->top);
    box.bottom = std::max(a->bottom, b->bottom);
    return box;
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

/* Up and down, which the arrow keys had no answer for at all.
 *
 * Equation Editor's own key table says these are ordinary caret motion -- the
 * same command group as left and right -- so a person expects Up inside a
 * denominator to reach the numerator.  There is no way to know that from the
 * tree: a numerator is not "above" a denominator in any structural sense, only
 * on the page.  So the page is what is asked.
 *
 * Among the stops whose band lies clear of this one in the wanted direction,
 * take the nearest band; within it, the nearest column.  That is what every
 * text editor does with a wrapped line, and it gives the right answer here for
 * fractions, scripts, limits and matrix rows without any of them being named. */
bool Equation::move_vertical(int dir, const SvgStyle& style) {
    Layout L = layout_math(*root_, style);
    std::string ct = caret_text();
    const size_t colon = ct.rfind(':');
    const std::string path = (colon == std::string::npos) ? ct : ct.substr(0, colon);
    const int index = (colon == std::string::npos) ? 0
                                                  : atoi(ct.c_str() + colon + 1);
    const CaretStop* cur = find_stop(L, path, index);
    if (!cur) return false;

    const double eps = 0.5;
    const CaretStop* best = nullptr;
    double bestGap = 0, bestDx = 0;
    for (const CaretStop& s : L.stops) {
        if (s.path == path) continue;              /* same slot is sideways */
        const double gap = (dir < 0) ? (cur->top - s.bottom)
                                     : (s.top - cur->bottom);
        if (gap < -eps) continue;                  /* not clear of us */
        const double dx = std::fabs(s.x - cur->x);
        if (!best || gap < bestGap - eps ||
            (std::fabs(gap - bestGap) <= eps && dx < bestDx)) {
            best = &s; bestGap = gap; bestDx = dx;
        }
    }
    if (!best) return false;
    clear_selection();
    set_caret_text(best->path + ":" + std::to_string(best->index));
    return true;
}

bool Equation::move_to_point(double x, double y, const SvgStyle& style) {
    clear_selection();
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
    take_selection();          /* typing replaces what is highlighted */
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
        if (style_ != "math" && node) style_one(*node, style_);
        l.insert(l.begin() + index_, std::move(node));
        ++index_;
    }
}

bool Equation::insert_symbol(const std::string& cmd) {
    int code = tex_command_to_unicode(cmd.c_str());
    if (code < 0) return false;
    checkpoint();
    take_selection();
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

bool Equation::insert_latex(const std::string& latex) {
    std::unique_ptr<LineNode> parsed = parse_latex(latex);
    if (!parsed || parsed->children.empty()) return false;

    checkpoint();
    take_selection();
    NodeList& l = here();
    clamp();
    const int n = int(parsed->children.size());
    for (int i = 0; i < n; ++i)
        l.insert(l.begin() + index_ + i, std::move(parsed->children[size_t(i)]));
    index_ += n;
    return true;
}

bool Equation::insert_template(const std::string& kind) {
    /* A selection is what the template WRAPS.  Select B, press the vector
     * chord, get \vec{B} -- which is how a vector actually gets written, and
     * the reason selection had to come before the Style menu. */
    NodeList wrapped = take_selection();
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
         * Ctrl+L after `x` gives `x` a subscript rather than an empty box --
         * unless a selection was made, which is a more explicit statement of
         * the same intent and takes precedence. */
        NodeList& l = here();
        if (wrapped.empty() && index_ > 0) {
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

    /* What was selected goes into the template's first slot -- the numerator
     * of a fraction, the body of a root, the base of a script, the thing an
     * accent sits over.  That is the slot node_slots() reports first, so no
     * table of "which slot wraps" is needed and none can fall behind a new
     * node type. */
    if (!wrapped.empty()) {
        auto slots = node_slots(*raw);
        if (!slots.empty()) {
            for (auto& n : wrapped) slots[0]->push_back(std::move(n));
        }
    }

    l.insert(l.begin() + index_, std::move(node));
    enter_first_empty_slot(*raw);
    return true;
}

bool Equation::backspace() {
    if (has_selection()) return delete_selection();
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
    if (has_selection()) return delete_selection();
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
        {"Ctrl+/",       "template.frac",      "fraction (skewed)"},
        {"Ctrl+T, S",    "template.sum",       "summation"},
        {"Ctrl+T, P",    "template.prod",      "product"},
        {"Ctrl+T, M",    "template.matrix2x2", "matrix"},
        /* Equation Editor keeps its bindings in a resource, not in code, and
         * the table below is what that resource says -- read out rather than
         * remembered.  Insert really is a second Tab there, and up and down
         * really are ordinary caret motion in the same command group as left
         * and right, which is why they belong here and not in some other
         * mode. */
        {"Tab",          "caret.next_slot",    "next slot"},
        {"Insert",       "caret.next_slot",    "next slot"},
        {"Shift+Tab",    "caret.prev_slot",    "previous slot"},
        {"Left",         "caret.left",         "left"},
        {"Right",        "caret.right",        "right"},
        {"Up",           "caret.up",           "up"},
        {"Down",         "caret.down",         "down"},
        {"Home",         "caret.home",         "start of slot"},
        {"End",          "caret.end",          "end of slot"},
        {"Backspace",    "edit.backspace",     "delete backwards / unwrap"},
        {"Delete",       "edit.delete",        "delete forwards"},
        {"Ctrl+Z",       "edit.undo",          "undo"},
        {"Ctrl+Y",       "edit.redo",          "redo"},
        /* Selection.  Equation Editor's Edit menu is Select All, Cut, Copy,
         * Paste and Clear; there is no menu bar here, so these are the whole
         * of it.  Shift stops at the slot edge rather than jumping out, so a
         * selection never silently changes what it is anchored to. */
        {"Shift+Left",   "select.left",        "extend selection left"},
        {"Shift+Right",  "select.right",       "extend selection right"},
        {"Shift+Home",   "select.home",        "extend to start of slot"},
        {"Shift+End",    "select.end",         "extend to end of slot"},
        {"Ctrl+A",       "select.all",         "select all"},
        /* Style, applied to the selection.  The letters are Equation Editor's
         * own menu mnemonics -- &Math, &Text, &Function, &Variable, &Greek,
         * Matri&x-Vector -- so the muscle memory carries over even though
         * there is no menu bar to press Alt+S in. */
        /* Equation Editor's, from the same resource.  Four of these were
         * guessed and four were wrong: it uses the letter that names the
         * EFFECT -- I for italic gives Variable, B for bold gives
         * Matrix-Vector -- and "=" for going back to plain Math.  Only Greek
         * and Function happened to match. */
        {"Ctrl+Shift+=", "style.math",         "math style"},
        {"Ctrl+Shift+E", "style.text",         "text style"},
        {"Ctrl+Shift+F", "style.function",     "function style"},
        {"Ctrl+Shift+I", "style.variable",     "variable style"},
        {"Ctrl+Shift+B", "style.vector",       "matrix-vector (bold)"},
        {"Ctrl+Shift+G", "style.greek",        "greek"},

        /* Move and select by WHOLE ITEM.  Left and Right walk into a
         * template, which is what you want for correcting a letter and not
         * what you want for stepping past a fraction to get at what follows
         * it: from the left of \frac{a}{b} it takes four presses to get past
         * it and Ctrl+Right takes one.  Equation Editor binds the same keys
         * to the same idea. */
        /* Enter breaks the line, which a stack can now be drawn for. */
        {"Enter",        "edit.newline",       "new line"},
        /* Format > Align.  Equation Editor's own three chords, read off the
         * same key table: Ctrl+Shift+L is command 3,0, Ctrl+Shift+C is 3,1 and
         * Ctrl+Shift+R is 3,2, and group 3 is the Format menu -- confirmed by
         * the Style group matching all six of its chords. */
        {"Ctrl+Shift+L", "format.left",        "align left"},
        {"Ctrl+Shift+C", "format.center",      "align centre"},
        {"Ctrl+Shift+R", "format.right",       "align right"},
        {"Ctrl+Left",    "caret.left_item",    "left one item"},
        {"Ctrl+Right",   "caret.right_item",   "right one item"},
        {"Ctrl+Shift+Left",  "select.left_item",  "select left one item"},
        {"Ctrl+Shift+Right", "select.right_item", "select right one item"},
    };

    /* Ctrl+K and a letter: one SYMBOL.
     *
     * Read out of Equation Editor's own key table rather than remembered.
     * The table stores these as "kind 4" records whose command number is the
     * code point itself -- 8706 is U+2202, 8734 is U+221E -- which is how the
     * pairing can be recovered at all: the label column in that resource is
     * unusable, but the numbers are not. */
    static const struct { const char* key; uint32_t cp; const char* name; } kSymbols[] = {
        {"T",       0x00D7, "times"},
        {"A",       0x2192, "right arrow"},
        {"D",       0x2202, "partial"},
        {"E",       0x2208, "element of"},
        {"Shift+E", 0x2209, "not an element of"},
        {"I",       0x221E, "infinity"},
        {"<",       0x2264, "less than or equal"},
        {">",       0x2265, "greater than or equal"},
        {"C",       0x2282, "subset of"},
        {"Shift+C", 0x2284, "not a subset of"},
    };

    /* Ctrl+G and a letter: ONE Greek letter, without changing the style.
     *
     * Equation Editor has both, and they are not the same thing.  Ctrl+Shift+G
     * is the STYLE -- everything typed after it is Greek, which is what you
     * want for a run of them.  Ctrl+G is a prefix for a single letter, which
     * is what you want for the mu in a line of otherwise Latin algebra, and is
     * the one people actually use.  Only the style existed here.
     *
     * The 52 chords are generated rather than written out: a hand-written
     * table of them would be 52 chances to typo a letter, and the mapping is
     * already in greek_of, where the style gets it. */
    static const std::vector<Binding>& kAll = [&]() -> const std::vector<Binding>& {
        static std::vector<std::string> pool;
        static std::vector<Binding> all;
        pool.reserve(62 * 3);
        all = kBindings;
        for (const auto& sym : kSymbols) {
            pool.push_back(std::string("Ctrl+K, ") + sym.key);
            pool.push_back("symbol." + std::to_string(sym.cp));
            pool.push_back(sym.name);
            const size_t i = pool.size() - 3;
            all.push_back({pool[i].c_str(), pool[i + 1].c_str(),
                           pool[i + 2].c_str()});
        }
        for (int upper = 0; upper < 2; ++upper) {
            for (char c = 'A'; c <= 'Z'; ++c) {
                const char latin = upper ? c : char(c - 'A' + 'a');
                if (!greek_of(uint32_t((unsigned char)latin))) continue;
                pool.push_back(std::string("Ctrl+G, ") +
                               (upper ? "Shift+" : "") + c);
                pool.push_back(std::string("greek.") + latin);
                pool.push_back(std::string("greek ") + latin);
                const size_t i = pool.size() - 3;
                all.push_back({pool[i].c_str(), pool[i + 1].c_str(),
                               pool[i + 2].c_str()});
            }
        }
        return all;
    }();
    return kAll;
}

bool Equation::command(const std::string& name) {
    if (name.compare(0, 9, "template.") == 0)
        return insert_template(name.substr(9));
    if (name == "caret.left")       return move_left();
    if (name == "caret.right")      return move_right();
    if (name == "caret.up")         return move_vertical(-1);
    if (name == "caret.down")       return move_vertical(+1);
    if (name == "caret.next_slot")  return next_slot();
    if (name == "caret.prev_slot")  return prev_slot();
    if (name == "caret.out")        return move_out();
    if (name == "caret.home")       { move_home(); return true; }
    if (name == "caret.end")        { move_end();  return true; }
    if (name == "edit.backspace")   return backspace();
    if (name == "edit.delete")      return erase();
    if (name == "select.left")      return extend_left();
    if (name == "select.right")     return extend_right();
    if (name == "select.home")      { extend_home(); return true; }
    if (name == "select.end")       { extend_end();  return true; }
    if (name == "select.all")       { select_all();  return true; }
    if (name.compare(0, 6, "style.") == 0) return set_style(name.substr(6));
    if (name.compare(0, 7, "format.") == 0) {
        const std::string how = name.substr(7);
        const int want = (how == "left")   ? 2
                       : (how == "right")  ? 3
                       : (how == "center") ? 0 : -1;
        if (want < 0) return false;
        /* The nearest enclosing stack, which is what the Format menu acts on.
         * Nothing to align outside one, and saying so is better than silently
         * doing nothing to the whole equation. */
        for (size_t depth = path_.size(); depth-- > 0; ) {
            std::vector<CaretStep> upto(path_.begin(), path_.begin() + depth);
            Node* n = parent_node(std::vector<CaretStep>(
                path_.begin(), path_.begin() + depth + 1));
            if (n && n->tag() == Node::kPile) {
                checkpoint();
                static_cast<PileNode&>(*n).halign = want;
                return true;
            }
        }
        return false;
    }
    if (name == "edit.newline")     return newline();
    if (name == "caret.left_item")  return move_item(-1, false);
    if (name == "caret.right_item") return move_item(+1, false);
    if (name == "select.left_item") return move_item(-1, true);
    if (name == "select.right_item")return move_item(+1, true);
    if (name.compare(0, 7, "symbol.") == 0) {
        const uint32_t cp = uint32_t(std::strtoul(name.c_str() + 7, nullptr, 10));
        if (!cp) return false;
        checkpoint();
        take_selection();
        NodeList& l = here();
        clamp();
        l.insert(l.begin() + index_,
                 make_char(typeface_for_code(uint16_t(cp)), uint16_t(cp)));
        ++index_;
        return true;
    }
    if (name.compare(0, 6, "greek.") == 0 && name.size() == 7) {
        const uint32_t g = greek_of(uint32_t((unsigned char)name[6]));
        if (!g) return false;
        checkpoint();
        take_selection();
        NodeList& l = here();
        clamp();
        l.insert(l.begin() + index_,
                 make_char(typeface_for_code(uint16_t(g)), uint16_t(g)));
        ++index_;
        return true;
    }
    if (name == "edit.undo")        return undo();
    if (name == "edit.redo")        return redo();
    return false;
}

}  // namespace mtef
