/*
 * equation_edit.cpp -- the equation editing model
 *
 * The caret is a path of (child, slot) steps from the root plus an index in
 * the slot it lands in.  Every operation is expressed on that pair, which is
 * what lets arrow keys walk *into* a fraction the way Equation Editor's do
 * instead of skipping over it as one opaque object.
 */
#include "equation_edit.h"
#include "tex_parser.h"
#include "latex_emitter.h"
#include "math_symbols.h"

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>

namespace eqnedit {

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
        case Node::kDirac: {
            auto& d = static_cast<DiracNode&>(n);
            return {&d.bra, &d.ket};
        }
        case Node::kLim: {
            auto& l = static_cast<LimNode&>(n);
            return {&l.content};
        }
        case Node::kMathbf: {
            auto& b = static_cast<MathbfNode&>(n);
            return {&b.content};
        }
        case Node::kGroup: {
            auto& g = static_cast<GroupNode&>(n);
            return {&g.children};
        }
        case Node::kOverset: {
            auto& o = static_cast<OversetNode&>(n);
            return {&o.base, &o.over};
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

constexpr const char* kDepthError = "maximum-nesting-depth";

int node_nesting(Node& node);

int list_nesting(NodeList& list) {
    int deepest = 0;
    for (auto& child : list)
        if (child) deepest = std::max(deepest, node_nesting(*child));
    return deepest;
}

/* Count structural containers below a slot.  LineNode is only a list wrapper;
 * every other node with slots consumes one level of the parser/editor limit. */
int node_nesting(Node& node) {
    const std::vector<NodeList*> slots = node_slots(node);
    if (slots.empty()) return 0;
    int deepest = 0;
    for (NodeList* slot : slots)
        if (slot) deepest = std::max(deepest, list_nesting(*slot));
    return (node.tag() == Node::kLine) ? deepest : 1 + deepest;
}

bool matrix_template_dimensions(const std::string& kind,
                                int* rows, int* cols) {
    if (kind.rfind("matrix", 0) != 0) return false;
    const size_t x = kind.find('x', 6);
    if (x == std::string::npos || x == 6 || x + 1 >= kind.size()) return false;
    char* rowEnd = nullptr;
    char* colEnd = nullptr;
    const long r = std::strtol(kind.c_str() + 6, &rowEnd, 10);
    const long c = std::strtol(kind.c_str() + x + 1, &colEnd, 10);
    if (rowEnd != kind.c_str() + x || colEnd != kind.c_str() + kind.size() ||
        r < 1 || r > 99 || c < 1 || c > 99) return false;
    if (rows) *rows = int(r);
    if (cols) *cols = int(c);
    return true;
}

std::string utf8_of(uint32_t cp) {
    std::string s;
    if (cp < 0x80) s += char(cp);
    else if (cp < 0x800) { s += char(0xC0 | (cp >> 6)); s += char(0x80 | (cp & 0x3F)); }
    else if (cp < 0x10000) {
        s += char(0xE0 | (cp >> 12)); s += char(0x80 | ((cp >> 6) & 0x3F));
        s += char(0x80 | (cp & 0x3F));
    } else {
        s += char(0xF0 | (cp >> 18)); s += char(0x80 | ((cp >> 12) & 0x3F));
        s += char(0x80 | ((cp >> 6) & 0x3F)); s += char(0x80 | (cp & 0x3F));
    }
    return s;
}

bool is_letter(char c) { return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z'); }
bool is_digit(char c)  { return c >= '0' && c <= '9'; }

bool is_automatic_function_name(const std::string& name) {
    /* Eqnedit32 string resources 11000--11038, plus common engineering
     * operators already understood by Eqnedit64's TeX emitter/parser. */
    static const char* kNames[] = {
        "Im", "Pr", "Re", "arg", "arcsin", "arccos", "arctan",
        "cosh", "cos", "coth", "cot", "cov", "csc", "deg", "det",
        "dim", "exp", "gcd", "glb", "hom", "inf", "int", "ker",
        "ln", "lg", "lim", "log", "lub", "max", "min", "mod", "sec",
        "sgn", "sinh", "sin", "sup", "tanh", "tan", "var",
        "curl", "div", "grad", "rot", "tr", "diag", "Res", "const",
        nullptr
    };
    for (int i = 0; kNames[i]; ++i)
        if (name == kNames[i]) return true;
    return false;
}

bool automatic_letter(const NodePtr& node) {
    if (!node || node->tag() != Node::kChar) return false;
    const auto* ch = static_cast<const CharNode*>(node.get());
    return ch->ch && is_letter(ch->ch) &&
           (ch->typeface == TF_VARIABLE || ch->automaticFunction);
}

/* End of the already-recognised function immediately before the caret,
 * allowing an ordinary variable suffix after it.  This is what makes typing
 * `sinx` behave as `\sin x`: once `sin` has become a function, the next
 * variable is its argument rather than a reason to turn `sin` italic again. */
int automatic_function_prefix_end(const NodeList& list, int caret) {
    int end = std::clamp(caret, 0, int(list.size()));
    while (end > 0 && automatic_letter(list[size_t(end - 1)]) &&
           !static_cast<const CharNode*>(list[size_t(end - 1)].get())
                ->automaticFunction)
        --end;
    if (end <= 0 || !automatic_letter(list[size_t(end - 1)]) ||
        !static_cast<const CharNode*>(list[size_t(end - 1)].get())
             ->automaticFunction)
        return -1;
    int first = end;
    while (first > 0 && automatic_letter(list[size_t(first - 1)]) &&
           static_cast<const CharNode*>(list[size_t(first - 1)].get())
               ->automaticFunction)
        --first;
    std::string name;
    for (int i = first; i < end; ++i)
        name.push_back(static_cast<const CharNode*>(list[size_t(i)].get())->ch);
    return is_automatic_function_name(name) ? end : -1;
}

std::vector<uint32_t> utf8_codes(const std::string& s) {
    std::vector<uint32_t> out;
    for (size_t i = 0; i < s.size();) {
        unsigned char c = static_cast<unsigned char>(s[i]);
        uint32_t cp = 0;
        size_t n = 1;
        if (c < 0x80) cp = c;
        else if ((c & 0xE0) == 0xC0 && i + 1 < s.size()) { cp = c & 0x1F; n = 2; }
        else if ((c & 0xF0) == 0xE0 && i + 2 < s.size()) { cp = c & 0x0F; n = 3; }
        else if ((c & 0xF8) == 0xF0 && i + 3 < s.size()) { cp = c & 0x07; n = 4; }
        else cp = 0xFFFD;
        if (n > 1) {
            bool valid = true;
            for (size_t j = 1; j < n; ++j) {
                unsigned char d = static_cast<unsigned char>(s[i + j]);
                if ((d & 0xC0) != 0x80) { valid = false; break; }
                cp = (cp << 6) | (d & 0x3F);
            }
            if (!valid) { cp = 0xFFFD; n = 1; }
        }
        out.push_back(cp);
        i += n;
    }
    return out;
}

bool same_path(const std::vector<CaretStep>& a,
               const std::vector<CaretStep>& b) {
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i)
        if (a[i].child != b[i].child || a[i].slot != b[i].slot) return false;
    return true;
}

std::unique_ptr<CharNode> make_char(int tf, uint32_t code, char ascii = 0) {
    auto c = std::make_unique<CharNode>();
    c->typeface = tf; c->charCode = code; c->ch = ascii;
    return c;
}

/* The slot holding the last item a reader can see inside `node`, descending
 * through nested templates.  Backspace uses it so one keystroke removes one
 * visible thing: without the descent, `\sqrt{x^{2}}` lost its x and its 2
 * together, and `c^{2}` lost the c along with the 2.  Null when the node has
 * no slots at all; the slot itself when every slot below is empty, so the
 * caller removes the empty shell. */
NodeList* deepest_nonempty_slot(Node& node) {
    NodeList* found = nullptr;
    Node* current = &node;
    for (;;) {
        auto slots = node_slots(*current);
        int last = -1;
        for (int i = int(slots.size()) - 1; i >= 0; --i)
            if (!slots[size_t(i)]->empty()) { last = i; break; }
        if (last < 0) return found;
        found = slots[size_t(last)];
        Node* next = found->back().get();
        if (!next || node_slots(*next).empty()) return found;
        current = next;
    }
}

/* The same descent as deepest_nonempty_slot, but recording how it got there.
 *
 * Backspace used to remove the last visible item and leave the caret outside
 * the structure it had just reached into, so the next character typed landed
 * after the whole thing: deleting the 2 of x_1^2 and typing 9 gave x_1^{}9
 * rather than putting the 9 back in the exponent.  The caret has to follow
 * the deletion into the slot, which means knowing the path and not just the
 * slot. */
NodeList* deepest_nonempty_path(NodeList& list, int childIndex,
                                std::vector<CaretStep>& steps) {
    steps.clear();
    if (childIndex < 0 || size_t(childIndex) >= list.size() ||
        !list[size_t(childIndex)])
        return nullptr;
    NodeList* found = nullptr;
    Node* current = list[size_t(childIndex)].get();
    int child = childIndex;
    for (;;) {
        auto slots = node_slots(*current);
        int last = -1;
        for (int i = int(slots.size()) - 1; i >= 0; --i)
            if (!slots[size_t(i)]->empty()) { last = i; break; }
        if (last < 0) return found;
        steps.push_back(CaretStep{child, last});
        found = slots[size_t(last)];
        Node* next = found->back().get();
        if (!next || node_slots(*next).empty()) return found;
        child = int(found->size()) - 1;
        current = next;
    }
}

/* The mirror of deepest_nonempty_path, for forward delete.  Written at the
 * same time as its Backspace twin on purpose: the last two defects in this
 * area were each fixed for one key and left standing for the other. */
NodeList* first_nonempty_path(NodeList& list, int childIndex,
                              std::vector<CaretStep>& steps) {
    steps.clear();
    if (childIndex < 0 || size_t(childIndex) >= list.size() ||
        !list[size_t(childIndex)])
        return nullptr;
    NodeList* found = nullptr;
    Node* current = list[size_t(childIndex)].get();
    int child = childIndex;
    for (;;) {
        auto slots = node_slots(*current);
        int first = -1;
        for (int i = 0; i < int(slots.size()); ++i)
            if (!slots[size_t(i)]->empty()) { first = i; break; }
        if (first < 0) return found;
        steps.push_back(CaretStep{child, first});
        found = slots[size_t(first)];
        Node* next = found->front().get();
        if (!next || node_slots(*next).empty()) return found;
        child = 0;
        current = next;
    }
}

/* The mirror of deepest_nonempty_slot, for forward delete: the slot holding
 * the first item a reader can see inside `node`.  Delete had the same defect
 * Backspace did -- `c^{2}` lost its c and its 2 to one press -- and it went
 * unnoticed because the rule was written for one key and the test exercised
 * only that key. */
NodeList* first_nonempty_slot(Node& node) {
    NodeList* found = nullptr;
    Node* current = &node;
    for (;;) {
        auto slots = node_slots(*current);
        int first = -1;
        for (int i = 0; i < int(slots.size()); ++i)
            if (!slots[size_t(i)]->empty()) { first = i; break; }
        if (first < 0) return found;
        found = slots[size_t(first)];
        Node* next = found->front().get();
        if (!next || node_slots(*next).empty()) return found;
        current = next;
    }
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

bool Equation::load_latex(const std::string& latex) {
    lastError_.clear();
    bool depthExceeded = false;
    std::unique_ptr<LineNode> parsed = parse_latex(latex, &depthExceeded);
    if (depthExceeded) {
        lastError_ = kDepthError;
        return false;
    }
    root_ = std::move(parsed);
    if (!root_) root_ = std::make_unique<LineNode>();
    path_.clear();
    index_ = int(root_->children.size());
    selecting_ = false;
    anchorPath_.clear();
    anchorIndex_ = 0;
    undo_.clear();
    redo_.clear();
    return true;
}

bool Equation::replace_latex(const std::string& latex, bool checkpointFirst) {
    lastError_.clear();
    bool depthExceeded = false;
    std::unique_ptr<LineNode> parsed = parse_latex(latex, &depthExceeded);
    if (depthExceeded) {
        lastError_ = kDepthError;
        return false;
    }
    if (checkpointFirst) checkpoint("TeX Edit");
    root_ = std::move(parsed);
    if (!root_) root_ = std::make_unique<LineNode>();
    path_.clear();
    index_ = int(root_->children.size());
    clear_selection();
    return true;
}

std::string Equation::latex() const { return tree_to_latex(*root_); }

std::string Equation::svg(const SvgStyle& style) const {
    return render_svg(*root_, style);
}

RenderMetrics Equation::metrics(const SvgStyle& style) const {
    return measure_equation(*root_, style);
}

bool Equation::caret_geometry(CaretGeometry* geometry,
                              const SvgStyle& style) const {
    return caret_geometry_equation(*root_, slot_at(path_), index_, style,
                                   geometry);
}

#ifdef _WIN32
void Equation::draw_gdi(HDC hdc, double left, double top, double scale,
                        const SvgStyle& style, bool show_placeholders,
                        bool show_caret, bool text_as_outlines) const {
    NodeList* caretSlot = slot_at(path_);
    NodeList* selectionSlot = nullptr;
    int first = -1, last = -1;
    selection_range(&selectionSlot, &first, &last);
    draw_equation_gdi(*root_, hdc, left, top, scale, style,
                      caretSlot, index_, selectionSlot, first, last,
                      show_placeholders, show_caret, text_as_outlines);
}
#endif

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

bool Equation::find_slot_path(const NodeList* target, NodeList& list,
                              std::vector<CaretStep>& path) const {
    if (&list == target) return true;
    for (size_t child = 0; child < list.size(); ++child) {
        if (!list[child]) continue;
        auto slots = node_slots(*list[child]);
        for (size_t slot = 0; slot < slots.size(); ++slot) {
            path.push_back({int(child), int(slot)});
            if (find_slot_path(target, *slots[slot], path)) return true;
            path.pop_back();
        }
    }
    return false;
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

bool Equation::selection_range(NodeList** slot, int* first, int* last) const {
    if (!selecting_ || !same_path(anchorPath_, path_) ||
        anchorIndex_ == index_) return false;
    NodeList* l = slot_at(path_);
    if (!l) return false;
    /* Clamp to the slot as it is NOW.  The anchor was recorded when the
     * selection began, and the slot can have shrunk since -- Undo rebuilds
     * the tree under an unchanged path -- so consumers erasing [first,last)
     * read past the end of the vector.  That is the heap corruption the
     * fuzzer reproduced from the user's crash (seed 14, op 625): ASan
     * placed the overflow in insert_text_typeface's replacing erase, fed by
     * exactly this stale range.  Clamping at the source fixes every
     * consumer, including the ones not written yet. */
    const int size = int(l->size());
    const int lo = std::min({anchorIndex_, index_, size});
    const int hi = std::min(std::max(anchorIndex_, index_), size);
    if (lo >= hi) return false;
    if (slot) *slot = l;
    if (first) *first = lo;
    if (last) *last = hi;
    return true;
}

bool Equation::has_selection() const {
    return selection_range(nullptr, nullptr, nullptr);
}

void Equation::begin_selection() {
    if (!has_selection()) {
        selecting_ = true;
        anchorPath_ = path_;
        anchorIndex_ = index_;
    }
}

bool Equation::select_step_left() {
    begin_selection();
    clamp();
    if (index_ <= 0) return false;      /* at the slot start: nothing more here */
    --index_;
    return true;
}

bool Equation::select_step_right() {
    begin_selection();
    clamp();
    if (size_t(index_) >= here().size()) return false;
    ++index_;
    return true;
}

void Equation::clear_selection() {
    selecting_ = false;
    anchorPath_.clear();
    anchorIndex_ = index_;
}

bool Equation::select_current_slot() {
    NodeList* slot = slot_at(path_);
    if (!slot || slot->empty()) {
        clear_selection();
        return false;
    }
    anchorPath_ = path_;
    anchorIndex_ = 0;
    index_ = int(slot->size());
    selecting_ = true;
    return true;
}

bool Equation::select_containing_structure() {
    if (path_.empty()) {
        clear_selection();
        return false;
    }
    const CaretStep container = path_.back();
    std::vector<CaretStep> parentPath(path_.begin(), path_.end() - 1);
    NodeList* parent = slot_at(parentPath);
    if (!parent || container.child < 0 ||
        size_t(container.child) >= parent->size()) {
        clear_selection();
        return false;
    }
    path_ = std::move(parentPath);
    anchorPath_ = path_;
    anchorIndex_ = container.child;
    index_ = container.child + 1;
    selecting_ = true;
    return true;
}

bool Equation::select_all() {
    path_.clear();
    return select_current_slot();
}

std::string Equation::selection_latex() const {
    NodeList* l = nullptr;
    int first = 0, last = 0;
    if (!selection_range(&l, &first, &last) || !l) return std::string();
    LaTeXEmitter em;
    return em.emit_range(*l, size_t(first), size_t(last));
}

bool Equation::delete_selection() {
    NodeList* l = nullptr;
    int first = 0, last = 0;
    if (!selection_range(&l, &first, &last) || !l) return false;
    checkpoint("Delete");
    l = slot_at(path_);
    if (!l) return false;
    first = std::max(0, std::min(first, int(l->size())));
    last = std::max(first, std::min(last, int(l->size())));
    l->erase(l->begin() + first, l->begin() + last);
    index_ = first;
    clear_selection();
    refresh_auto_function_word();
    return true;
}

bool Equation::hit_test(double x_points, double y_points,
                        const SvgStyle& style, bool extend) {
    const NodeList* target = nullptr;
    int targetIndex = 0;
    if (!hit_test_equation(*root_, x_points, y_points, style,
                           &target, &targetIndex)) return false;
    if (extend) begin_selection();
    else clear_selection();
    std::vector<CaretStep> p;
    if (!find_slot_path(target, root_->children, p)) return false;
    path_ = std::move(p);
    index_ = targetIndex;
    clamp();
    /* A drag that crosses structural slots cannot be represented as one
     * contiguous tree range.  Keep the caret but drop the ambiguous range. */
    if (selecting_ && !same_path(anchorPath_, path_)) clear_selection();
    return true;
}

bool Equation::matrix_context(size_t* depth, MatrixNode** matrix,
                              int* cell, int layoutKind) const {
    for (size_t j = path_.size(); j-- > 0;) {
        std::vector<CaretStep> up(path_.begin(), path_.begin() + j);
        NodeList* list = slot_at(up);
        if (!list) continue;
        const CaretStep& st = path_[j];
        if (st.child < 0 || size_t(st.child) >= list->size() || !(*list)[st.child])
            continue;
        if ((*list)[st.child]->tag() != Node::kMatrix) continue;
        auto* m = static_cast<MatrixNode*>((*list)[st.child].get());
        if (layoutKind >= 0 && int(m->layoutKind) != layoutKind) continue;
        if (depth) *depth = j;
        if (matrix) *matrix = m;
        if (cell) *cell = st.slot;
        return true;
    }
    return false;
}

bool Equation::is_multiline() const {
    std::vector<NodeList*> pending;
    pending.push_back(const_cast<NodeList*>(&root_->children));
    while (!pending.empty()) {
        NodeList* list = pending.back(); pending.pop_back();
        for (auto& n : *list) {
            if (!n) continue;
            if (n->tag() == Node::kMatrix) {
                auto* m = static_cast<MatrixNode*>(n.get());
                if (m->layoutKind == MatrixNode::kAlignedLayout && m->rows > 1)
                    return true;
            }
            for (NodeList* s : node_slots(*n)) pending.push_back(s);
        }
    }
    return false;
}

bool Equation::new_line() {
    size_t depth = 0;
    MatrixNode* matrix = nullptr;
    int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell,
                        MatrixNode::kAlignedLayout)) {
        checkpoint("Matrix");
        auto aligned = std::make_unique<MatrixNode>();
        aligned->layoutKind = MatrixNode::kAlignedLayout;
        aligned->rows = 2; aligned->cols = 1;
        auto first = std::make_unique<LineNode>();
        first->children = std::move(root_->children);
        aligned->elements.push_back(std::move(first));
        aligned->elements.push_back(std::make_unique<LineNode>());
        root_->children.push_back(std::move(aligned));
        path_.clear(); path_.push_back({0, 1});
        index_ = 0; clear_selection();
        return true;
    }

    checkpoint("Matrix");
    matrix_context(&depth, &matrix, &cell, MatrixNode::kAlignedLayout);
    if (!matrix || matrix->cols <= 0) return false;
    int row = std::clamp(cell / matrix->cols, 0, std::max(0, matrix->rows - 1));
    int col = std::clamp(cell % matrix->cols, 0, matrix->cols - 1);
    size_t at = size_t(row + 1) * size_t(matrix->cols);
    for (int c = 0; c < matrix->cols; ++c)
        matrix->elements.insert(matrix->elements.begin() + at + size_t(c),
                                std::make_unique<LineNode>());
    ++matrix->rows;
    path_.resize(depth + 1);
    path_[depth].slot = (row + 1) * matrix->cols + col;
    index_ = 0; clear_selection();
    return true;
}

bool Equation::alignment_tab() {
    size_t depth = 0;
    MatrixNode* matrix = nullptr;
    int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell,
                        MatrixNode::kAlignedLayout)) {
        if (!path_.empty()) return false;
        checkpoint("Alignment");
        auto aligned = std::make_unique<MatrixNode>();
        aligned->layoutKind = MatrixNode::kAlignedLayout;
        aligned->rows = 1; aligned->cols = 2;
        auto left = std::make_unique<LineNode>();
        auto right = std::make_unique<LineNode>();
        int split = std::clamp(index_, 0, int(root_->children.size()));
        for (int i = 0; i < split; ++i)
            left->children.push_back(std::move(root_->children[size_t(i)]));
        for (size_t i = size_t(split); i < root_->children.size(); ++i)
            right->children.push_back(std::move(root_->children[i]));
        root_->children.clear();
        aligned->elements.push_back(std::move(left));
        aligned->elements.push_back(std::move(right));
        root_->children.push_back(std::move(aligned));
        path_.clear(); path_.push_back({0, 1});
        index_ = 0; clear_selection();
        return true;
    }
    if (!matrix || matrix->cols <= 0) return false;
    int row = cell / matrix->cols, col = cell % matrix->cols;
    if (col + 1 < matrix->cols) {
        path_.resize(depth + 1);
        path_[depth].slot = cell + 1;
        index_ = 0; clear_selection();
        return true;
    }

    checkpoint("Alignment");
    matrix_context(&depth, &matrix, &cell, MatrixNode::kAlignedLayout);
    row = cell / matrix->cols; col = cell % matrix->cols;
    const int oldCols = matrix->cols;
    NodeList rebuilt;
    rebuilt.reserve(size_t(matrix->rows) * size_t(oldCols + 1));
    for (int r = 0; r < matrix->rows; ++r) {
        for (int c = 0; c < oldCols; ++c) {
            size_t oldIndex = size_t(r * oldCols + c);
            rebuilt.push_back(std::move(matrix->elements[oldIndex]));
            if (c == col) {
                auto added = std::make_unique<LineNode>();
                if (r == row && depth + 1 == path_.size() &&
                    rebuilt.back() && rebuilt.back()->tag() == Node::kLine) {
                    auto& cur = static_cast<LineNode&>(*rebuilt.back()).children;
                    int split = std::clamp(index_, 0, int(cur.size()));
                    for (size_t i = size_t(split); i < cur.size(); ++i)
                        added->children.push_back(std::move(cur[i]));
                    cur.erase(cur.begin() + split, cur.end());
                }
                rebuilt.push_back(std::move(added));
            }
        }
    }
    matrix->elements = std::move(rebuilt);
    matrix->cols = oldCols + 1;
    path_.resize(depth + 1);
    path_[depth].slot = row * matrix->cols + col + 1;
    index_ = 0; clear_selection();
    return true;
}

bool Equation::move_up() {
    size_t depth = 0; MatrixNode* m = nullptr; int cell = 0;
    if (!matrix_context(&depth, &m, &cell) || !m || m->cols <= 0) return false;
    int row = cell / m->cols, col = cell % m->cols;
    if (row <= 0) return false;
    path_.resize(depth + 1); path_[depth].slot = (row - 1) * m->cols + col;
    clamp(); return true;
}

bool Equation::move_down() {
    size_t depth = 0; MatrixNode* m = nullptr; int cell = 0;
    if (!matrix_context(&depth, &m, &cell) || !m || m->cols <= 0) return false;
    int row = cell / m->cols, col = cell % m->cols;
    if (row + 1 >= m->rows) return false;
    path_.resize(depth + 1); path_[depth].slot = (row + 1) * m->cols + col;
    clamp(); return true;
}

bool Equation::matrix_dimensions(int* rows, int* cols) const {
    MatrixNode* matrix = nullptr;
    if (!matrix_context(nullptr, &matrix, nullptr, MatrixNode::kMatrixLayout) ||
        !matrix) return false;
    if (rows) *rows = matrix->rows;
    if (cols) *cols = matrix->cols;
    return matrix->rows > 0 && matrix->cols > 0;
}

bool Equation::matrix_add_row() {
    size_t depth = 0; MatrixNode* matrix = nullptr; int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell, MatrixNode::kMatrixLayout) ||
        !matrix || matrix->rows <= 0 || matrix->cols <= 0 ||
        matrix->rows >= 99) return false;
    const int row = std::clamp(cell / matrix->cols, 0, matrix->rows - 1);
    const int col = std::clamp(cell % matrix->cols, 0, matrix->cols - 1);
    checkpoint("Matrix Row");
    const size_t at = size_t(row + 1) * size_t(matrix->cols);
    for (int c = 0; c < matrix->cols; ++c)
        matrix->elements.insert(matrix->elements.begin() + at + size_t(c),
                                std::make_unique<LineNode>());
    ++matrix->rows;
    path_.resize(depth + 1);
    path_[depth].slot = (row + 1) * matrix->cols + col;
    index_ = 0; clear_selection();
    return true;
}

bool Equation::matrix_remove_row() {
    size_t depth = 0; MatrixNode* matrix = nullptr; int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell, MatrixNode::kMatrixLayout) ||
        !matrix || matrix->rows <= 1 || matrix->cols <= 0) return false;
    const int row = std::clamp(cell / matrix->cols, 0, matrix->rows - 1);
    const int col = std::clamp(cell % matrix->cols, 0, matrix->cols - 1);
    checkpoint("Matrix Row");
    const size_t first = size_t(row) * size_t(matrix->cols);
    matrix->elements.erase(matrix->elements.begin() + first,
                           matrix->elements.begin() + first + matrix->cols);
    --matrix->rows;
    const int targetRow = std::min(row, matrix->rows - 1);
    path_.resize(depth + 1);
    path_[depth].slot = targetRow * matrix->cols + col;
    index_ = 0; clear_selection();
    return true;
}

bool Equation::matrix_add_column() {
    size_t depth = 0; MatrixNode* matrix = nullptr; int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell, MatrixNode::kMatrixLayout) ||
        !matrix || matrix->rows <= 0 || matrix->cols <= 0 ||
        matrix->cols >= 99) return false;
    const int row = std::clamp(cell / matrix->cols, 0, matrix->rows - 1);
    const int col = std::clamp(cell % matrix->cols, 0, matrix->cols - 1);
    checkpoint("Matrix Column");
    const int oldCols = matrix->cols;
    NodeList rebuilt;
    rebuilt.reserve(size_t(matrix->rows) * size_t(oldCols + 1));
    for (int r = 0; r < matrix->rows; ++r) {
        for (int c = 0; c < oldCols; ++c) {
            rebuilt.push_back(std::move(
                matrix->elements[size_t(r * oldCols + c)]));
            if (c == col) rebuilt.push_back(std::make_unique<LineNode>());
        }
    }
    matrix->elements = std::move(rebuilt);
    matrix->cols = oldCols + 1;
    path_.resize(depth + 1);
    path_[depth].slot = row * matrix->cols + col + 1;
    index_ = 0; clear_selection();
    return true;
}

bool Equation::matrix_remove_column() {
    size_t depth = 0; MatrixNode* matrix = nullptr; int cell = 0;
    if (!matrix_context(&depth, &matrix, &cell, MatrixNode::kMatrixLayout) ||
        !matrix || matrix->rows <= 0 || matrix->cols <= 1) return false;
    const int row = std::clamp(cell / matrix->cols, 0, matrix->rows - 1);
    const int col = std::clamp(cell % matrix->cols, 0, matrix->cols - 1);
    checkpoint("Matrix Column");
    const int oldCols = matrix->cols;
    NodeList rebuilt;
    rebuilt.reserve(size_t(matrix->rows) * size_t(oldCols - 1));
    for (int r = 0; r < matrix->rows; ++r)
        for (int c = 0; c < oldCols; ++c)
            if (c != col)
                rebuilt.push_back(std::move(
                    matrix->elements[size_t(r * oldCols + c)]));
    matrix->elements = std::move(rebuilt);
    matrix->cols = oldCols - 1;
    const int targetCol = std::min(col, matrix->cols - 1);
    path_.resize(depth + 1);
    path_[depth].slot = row * matrix->cols + targetCol;
    index_ = 0; clear_selection();
    return true;
}

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

/* ------------------------------------------------------------- history */

void Equation::checkpoint(const char* action) {
    undo_.push_back({latex(), caret_text(), action ? action : "Edit"});
    redo_.clear();
    if (undo_.size() > 200) undo_.erase(undo_.begin());
}

void Equation::restore(const Snapshot& s) {
    root_ = parse_latex(s.latex);
    if (!root_) root_ = std::make_unique<LineNode>();
    set_caret_text(s.caret);
    clear_selection();
}

bool Equation::can_undo() const { return !undo_.empty(); }
bool Equation::can_redo() const { return !redo_.empty(); }

std::string Equation::undo_name() const {
    return undo_.empty() ? std::string() : undo_.back().action;
}

std::string Equation::redo_name() const {
    return redo_.empty() ? std::string() : redo_.back().action;
}

bool Equation::undo() {
    if (undo_.empty()) return false;
    Snapshot s = undo_.back();
    undo_.pop_back();
    redo_.push_back({latex(), caret_text(), s.action});
    restore(s);
    return true;
}

bool Equation::redo() {
    if (redo_.empty()) return false;
    Snapshot s = redo_.back();
    redo_.pop_back();
    undo_.push_back({latex(), caret_text(), s.action});
    restore(s);
    return true;
}

/* ------------------------------------------------------------- editing */

void Equation::insert_text(const std::string& utf8) {
    insert_text_typeface(utf8, -1);
}

bool Equation::insert_styled_text(const std::string& utf8,
                                  const std::string& style) {
    int typeface = -1;
    if (style == "text") typeface = TF_TEXT;
    else if (style == "function") typeface = TF_FUNCTION;
    else if (style == "vector") typeface = TF_VECTOR;
    else if (style == "variable") typeface = TF_VARIABLE;
    else return false;
    insert_text_typeface(utf8, typeface);
    return true;
}

bool Equation::restyle_selection(const std::string& style) {
    int typeface = -1;
    if (style == "text") typeface = TF_TEXT;
    else if (style == "function") typeface = TF_FUNCTION;
    else if (style == "vector") typeface = TF_VECTOR;
    else if (style == "variable") typeface = TF_VARIABLE;
    else return false;

    NodeList* slot = nullptr;
    int first = 0, last = 0;
    if (!selection_range(&slot, &first, &last) || !slot) return false;
    checkpoint("Style Change");
    slot = slot_at(path_);
    if (!slot) return false;
    first = std::max(0, std::min(first, int(slot->size())));
    last = std::max(first, std::min(last, int(slot->size())));

    /* Apply to every character in the range, descending into any structure so
     * the whole selected sub-expression takes the style. */
    std::vector<NodeList*> pending;
    for (int i = first; i < last; ++i) {
        Node* n = (*slot)[size_t(i)].get();
        if (!n) continue;
        if (n->tag() == Node::kChar) {
            auto& ch = static_cast<CharNode&>(*n);
            ch.typeface = typeface;
            ch.automaticFunction = false;
        } else {
            for (NodeList* s : node_slots(*n)) pending.push_back(s);
        }
    }
    while (!pending.empty()) {
        NodeList* s = pending.back();
        pending.pop_back();
        if (!s) continue;
        for (auto& np : *s) {
            if (!np) continue;
            if (np->tag() == Node::kChar) {
                auto& ch = static_cast<CharNode&>(*np);
                ch.typeface = typeface;
                ch.automaticFunction = false;
            } else {
                for (NodeList* c : node_slots(*np)) pending.push_back(c);
            }
        }
    }
    return true;
}

void Equation::refresh_auto_function_word(int leftBoundary) {
    NodeList& l = here();
    clamp();

    int first = index_;
    const int floor = std::max(0, leftBoundary);
    while (first > floor && automatic_letter(l[size_t(first - 1)])) --first;
    int last = index_;
    while (last < int(l.size()) && automatic_letter(l[size_t(last)])) ++last;
    if (first == last) return;

    std::string word;
    word.reserve(size_t(last - first));
    for (int i = first; i < last; ++i)
        word.push_back(static_cast<CharNode*>(l[size_t(i)].get())->ch);
    const bool recognised = is_automatic_function_name(word);
    for (int i = first; i < last; ++i) {
        auto* ch = static_cast<CharNode*>(l[size_t(i)].get());
        if (recognised) {
            ch->typeface = TF_FUNCTION;
            ch->automaticFunction = true;
        } else if (ch->automaticFunction) {
            ch->typeface = TF_VARIABLE;
            ch->automaticFunction = false;
        }
    }
}

void Equation::insert_text_typeface(const std::string& utf8,
                                    int forcedTypeface) {
    NodeList* selected = nullptr;
    int first = 0, last = 0;
    bool replacing = selection_range(&selected, &first, &last);
    int functionPrefixEnd = -1;
    if (!replacing && forcedTypeface < 0) {
        NodeList& before = here();
        clamp();
        functionPrefixEnd = automatic_function_prefix_end(before, index_);
    }
    checkpoint("Typing");
    if (replacing) {
        selected = slot_at(path_);
        if (selected) {
            selected->erase(selected->begin() + first, selected->begin() + last);
            index_ = first;
        }
        clear_selection();
    }
    NodeList& l = here();
    clamp();
    for (uint32_t cp : utf8_codes(utf8)) {
        std::unique_ptr<CharNode> node;
        char ascii = cp < 0x80 ? char(cp) : 0;
        if (forcedTypeface >= 0)
            node = make_char(forcedTypeface, cp, ascii);
        else if (ascii && is_digit(ascii))
            node = make_char(TF_NUMBER, cp, ascii);
        else if (ascii && is_letter(ascii))
            node = make_char(TF_VARIABLE, cp, ascii);
        else if (cp == '-')
            node = make_char(TF_SYMBOL, 0x2212, ascii);
        else {
            int tf = typeface_for_code(cp);
            if (cp >= 0x80 && tf == TF_SYMBOL && !is_latex_symbol_codepoint(cp))
                tf = TF_TEXT;
            node = make_char(tf, cp, ascii);
        }
        l.insert(l.begin() + index_, std::move(node));
        ++index_;
    }
    if (forcedTypeface < 0) {
        /* `sin` + `h` must still grow into the longer recognised `sinh`.
         * Otherwise retain the complete function and classify only the newly
         * typed argument suffix. */
        if (functionPrefixEnd >= 0) {
            int firstWord = index_;
            while (firstWord > 0 && automatic_letter(l[size_t(firstWord - 1)]))
                --firstWord;
            int lastWord = index_;
            while (lastWord < int(l.size()) &&
                   automatic_letter(l[size_t(lastWord)]))
                ++lastWord;
            std::string whole;
            for (int i = firstWord; i < lastWord; ++i)
                whole.push_back(static_cast<CharNode*>(l[size_t(i)].get())->ch);
            if (is_automatic_function_name(whole)) functionPrefixEnd = -1;
        }
        refresh_auto_function_word(functionPrefixEnd);
    }
}

bool Equation::insert_latex(const std::string& latex) {
    lastError_.clear();
    bool depthExceeded = false;
    std::unique_ptr<LineNode> parsed = parse_latex(latex, &depthExceeded);
    if (!parsed) return false;
    if (depthExceeded ||
        int(path_.size()) + list_nesting(parsed->children) > kMaxNestingDepth) {
        lastError_ = kDepthError;
        return false;
    }
    NodeList* selected = nullptr;
    int first = 0, last = 0;
    bool replacing = selection_range(&selected, &first, &last);
    checkpoint("Paste");
    if (replacing) {
        selected = slot_at(path_);
        if (selected) {
            selected->erase(selected->begin() + first, selected->begin() + last);
            index_ = first;
        }
        clear_selection();
    }
    NodeList& l = here();
    clamp();
    int at = index_;
    for (auto& n : parsed->children)
        if (n) l.insert(l.begin() + at++, std::move(n));
    index_ = at;
    return true;
}

bool Equation::insert_symbol(const std::string& cmd) {
    lastError_.clear();
    int code = latex_symbol_codepoint(cmd);
    if (code < 0) {
        lastError_ = "unknown-symbol";
        return false;
    }
    NodeList* selected = nullptr;
    int first = 0, last = 0;
    bool replacing = selection_range(&selected, &first, &last);
    checkpoint("Symbol");
    if (replacing) {
        selected = slot_at(path_);
        if (selected) {
            selected->erase(selected->begin() + first, selected->begin() + last);
            index_ = first;
        }
        clear_selection();
    }
    NodeList& l = here();
    clamp();
    auto ch = make_char(typeface_for_code(uint32_t(code)), uint32_t(code));
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
    lastError_.clear();
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
        node = std::move(s);
    } else if (kind == "paren")   { node = fence(tmPAREN); }
    else if (kind == "bracket")   { node = fence(tmBRACK); }
    else if (kind == "brace")     { node = fence(tmBRACE); }
    else if (kind == "abs")       { node = fence(tmBAR); }
    else if (kind == "angle")     { node = fence(tmANGLE); }
    else if (kind == "norm")      { node = fence(tmDBAR); }
    else if (kind == "floor")     { node = fence(tmFLOOR); }
    else if (kind == "ceil")      { node = fence(tmCEIL); }
    else if (kind == "dirac")     { node = std::make_unique<DiracNode>(); }
    /* No oiint or oiiint here on purpose.  \oiint needs esint and \oiiint
     * exists only in packages that replace the whole math font, so offering
     * either would hand the user TeX their paper cannot typeset.  Both are
     * still read, so a pasted equation keeps them; we simply never originate
     * one.  tests/test_tex_compiles.py is what holds this line. */
    else if (kind == "int" || kind == "oint" || kind == "iint" ||
             kind == "iiint") {
        auto i = std::make_unique<IntegralNode>();
        i->selector = (kind == "oint")   ? tmSSINT
                    : (kind == "iint")   ? tmDINT
                    : (kind == "iiint")  ? tmTINT : tmSINT;
        i->hasLower = i->hasUpper = true;
        node = std::move(i);
    } else if (kind == "sum" || kind == "prod" || kind == "coprod" ||
               kind == "bigcup" || kind == "bigcap") {
        auto b = std::make_unique<BigOpNode>();
        b->selector = (kind == "sum")    ? tmSUM
                    : (kind == "prod")   ? tmPROD
                    : (kind == "coprod") ? tmCOPROD
                    : (kind == "bigcup") ? tmUNION : tmINTER;
        b->hasLower = b->hasUpper = true;
        b->hasLimits = true;
        node = std::move(b);
    } else if (kind == "under" || kind == "over") {
        auto o = std::make_unique<OversetNode>();
        o->under = (kind == "under");
        node = std::move(o);
    } else if (kind == "overbrace" || kind == "underbrace") {
        auto b = std::make_unique<BraceDecoNode>();
        b->selector = (kind == "overbrace") ? tmUHBRACE : tmLHBRACE;
        node = std::move(b);
    } else if (kind == "lim") {
        /* A limit sets its condition underneath in display style, so it is a
         * LimNode (which the renderer stacks) and round-trips as \lim_{...}.
         * There was no way to enter one before -- no template, no button, no
         * shortcut, only typing \lim into the source.  The caret lands in the
         * empty condition slot, ready for x \to 0. */
        node = std::make_unique<LimNode>();
    } else if (kind == "overline" || kind == "underline" ||
               kind == "overrightarrow" || kind == "overleftarrow" ||
               kind == "overleftrightarrow") {
        /* Decorations stretch over a whole slot; the accents below sit on a
         * single character.  TeX keeps the two apart (\overline vs \bar,
         * \overrightarrow vs \vec) and so must we, or one becomes the other
         * on the next save. */
        auto d = std::make_unique<DecorationNode>();
        d->selector = (kind == "overline")       ? tmOBAR
                    : (kind == "underline")      ? tmUBAR
                    : (kind == "overrightarrow") ? tmRARROW
                    : (kind == "overleftarrow")  ? tmLARROW : tmBARROW;
        node = std::move(d);
    } else if (kind == "hat" || kind == "vec" || kind == "bar" ||
               kind == "tilde" || kind == "dot" || kind == "ddot" ||
               kind == "dddot" || kind == "prime" || kind == "dprime" ||
               kind == "tprime" || kind == "strike" || kind == "frown" ||
               kind == "smile") {
        auto e = std::make_unique<EmbellNode>();
        e->embellType = (kind == "hat")    ? EM_HAT
                      : (kind == "vec")    ? EM_RARROW
                      : (kind == "bar")    ? EM_OBAR
                      : (kind == "tilde")  ? EM_TILDE
                      : (kind == "dot")    ? EM_DOT
                      : (kind == "ddot")   ? EM_DDOT
                      : (kind == "dddot")  ? EM_TDOT
                      : (kind == "prime")  ? EM_PRIME
                      : (kind == "dprime") ? EM_DPRIME
                      : (kind == "tprime") ? EM_TPRIME
                      : (kind == "strike") ? EM_MBAR
                      : (kind == "frown")  ? EM_FROWN : EM_SMILE;
        node = std::move(e);
    } else {
        int matrixRows = 0, matrixCols = 0;
        if (matrix_template_dimensions(kind, &matrixRows, &matrixCols)) {
            auto m = std::make_unique<MatrixNode>();
            m->rows = matrixRows; m->cols = matrixCols;
            for (int i = 0; i < matrixRows * matrixCols; ++i)
                m->elements.push_back(std::make_unique<LineNode>());
            node = std::move(m);
        } else if (kind == "cases") {
            /* Same shape the parser builds for \begin{cases}: a left brace
             * around a two-column table, value in the first column and
             * condition in the second.  Tab walks the four cells. */
            auto m = std::make_unique<MatrixNode>();
            m->layoutKind = MatrixNode::kCasesLayout;
            m->rows = 2; m->cols = 2;
            for (int i = 0; i < 4; ++i)
                m->elements.push_back(std::make_unique<LineNode>());
            auto f = fence(tmBRACE);
            f->variation = 1;
            f->content.push_back(std::move(m));
            node = std::move(f);
        } else {
            lastError_ = "unknown-template";
            return false;
        }
    }

    auto primary_slot = [](Node& n) -> NodeList* {
        switch (n.tag()) {
            case Node::kFrac:       return &static_cast<FracNode&>(n).numer;
            case Node::kSqrt:       return &static_cast<SqrtNode&>(n).content;
            case Node::kScript:     return &static_cast<ScriptNode&>(n).base;
            case Node::kFence:      return &static_cast<FenceNode&>(n).content;
            case Node::kIntegral:   return &static_cast<IntegralNode&>(n).body;
            case Node::kBigOp:      return &static_cast<BigOpNode&>(n).body;
            case Node::kEmbell:     return &static_cast<EmbellNode&>(n).content;
            case Node::kDecoration: return &static_cast<DecorationNode&>(n).content;
            case Node::kOverset:    return &static_cast<OversetNode&>(n).base;
            case Node::kBraceDeco:  return &static_cast<BraceDecoNode&>(n).content;
            case Node::kDirac:      return &static_cast<DiracNode&>(n).bra;
            case Node::kLim:        return &static_cast<LimNode&>(n).content;
            default: return nullptr;
        }
    };

    NodeList* selected = nullptr;
    int first = 0, last = 0;
    const bool hasRange = selection_range(&selected, &first, &last);
    bool wrapSelection = hasRange && primary_slot(*node) != nullptr;
    const bool wrapPrevious = !wrapSelection &&
        (node->tag() == Node::kScript || node->tag() == Node::kEmbell ||
         node->tag() == Node::kDecoration || node->tag() == Node::kOverset) &&
        index_ > 0;

    int wrappedDepth = 0;
    if (wrapSelection) {
        for (int i = first; i < last; ++i)
            if ((*selected)[size_t(i)])
                wrappedDepth = std::max(
                    wrappedDepth, node_nesting(*(*selected)[size_t(i)]));
    } else if (wrapPrevious) {
        NodeList& current = here();
        if (index_ <= int(current.size()) && current[size_t(index_ - 1)])
            wrappedDepth = node_nesting(*current[size_t(index_ - 1)]);
    }
    int prospectiveDepth = node_nesting(*node);
    if (wrapSelection || wrapPrevious)
        prospectiveDepth = std::max(prospectiveDepth, 1 + wrappedDepth);
    if (int(path_.size()) + prospectiveDepth > kMaxNestingDepth) {
        lastError_ = kDepthError;
        return false;
    }

    checkpoint("Template");
    NodeList& l = here();
    clamp();
    NodeList wrapped;
    if (wrapSelection) {
        for (int i = first; i < last; ++i)
            wrapped.push_back(std::move(l[size_t(i)]));
        l.erase(l.begin() + first, l.begin() + last);
        index_ = first;
        clear_selection();
    } else if (wrapPrevious) {
        /* Scripts and accents apply to the item just typed, matching the old
         * editor's Ctrl+L/Ctrl+- muscle memory.  With no preceding item they
         * remain ordinary empty templates ready for input.  Checkpoint first
         * so Undo restores the unwrapped base. */
        wrapped.push_back(std::move(l[size_t(index_ - 1)]));
        l.erase(l.begin() + index_ - 1);
        --index_;
        clear_selection();
    } else if (has_selection()) {
        /* Templates without a meaningful wrapper slot (matrix/cases) replace
         * the selected range, matching ordinary editor insertion. */
        for (int i = first; i < last; ++i)
            l[size_t(i)].reset();
        l.erase(l.begin() + first, l.begin() + last);
        index_ = first;
        clear_selection();
    }
    if (!wrapped.empty()) {
        NodeList* dst = primary_slot(*node);
        if (dst) *dst = std::move(wrapped);
    }
    Node* raw = node.get();
    l.insert(l.begin() + index_, std::move(node));
    enter_first_empty_slot(*raw);
    return true;
}

bool Equation::backspace() {
    if (has_selection()) return delete_selection();
    clamp();
    if (index_ > 0) {
        NodeList& current = here();
        const size_t target = size_t(index_ - 1);
        if (current[target] && current[target]->tag() == Node::kScript) {
            auto& script = static_cast<ScriptNode&>(*current[target]);
            const bool emptySub = script.hasSub && script.sub.empty();
            const bool emptySup = script.hasSup && script.sup.empty();
            if (emptySub || emptySup) {
                checkpoint("Delete");
                NodeList& cur = here();
                auto& editable = static_cast<ScriptNode&>(*cur[target]);
                if (editable.hasSub && editable.sub.empty()) editable.hasSub = false;
                if (editable.hasSup && editable.sup.empty()) editable.hasSup = false;
                if (!editable.hasSub && !editable.hasSup) {
                    NodeList base = std::move(editable.base);
                    cur.erase(cur.begin() + target);
                    size_t at = target;
                    for (auto& node : base)
                        cur.insert(cur.begin() + at++, std::move(node));
                    index_ = int(at);
                }
                refresh_auto_function_word();
                return true;
            }
        }
        /* A template is not one keystroke's worth of deletion.  Erasing the
         * node outright took the `c` and the `2` out of `c^{2}` together, so
         * one Backspace destroyed two things the reader could see.  Take only
         * the last item of the last slot that still holds something, and move
         * the caret into that slot: repeating Backspace peels the structure
         * one visible item at a time, and whatever is typed next goes back
         * where the deletion happened.  Leaving the caret outside meant that
         * after deleting the 2 of x_1^2, typing put the character after the
         * whole script instead of into the exponent it had just emptied. */
        if (current[target] && deepest_nonempty_slot(*current[target])) {
            checkpoint("Delete");
            NodeList& fresh = here();
            std::vector<CaretStep> steps;
            if (NodeList* slot = deepest_nonempty_path(fresh, int(target),
                                                       steps)) {
                slot->pop_back();
                for (const CaretStep& step : steps) path_.push_back(step);
                index_ = int(slot->size());
            }
            refresh_auto_function_word();
            return true;
        }

        int functionPrefixEnd = -1;
        if (current[target] && current[target]->tag() == Node::kChar &&
            !static_cast<const CharNode*>(current[target].get())
                 ->automaticFunction)
            functionPrefixEnd = automatic_function_prefix_end(current, index_);
        checkpoint("Delete");
        NodeList& cur = here();
        cur.erase(cur.begin() + index_ - 1);
        --index_;
        refresh_auto_function_word(functionPrefixEnd);
        return true;
    }
    /* At the start of an empty slot, backspace unwraps the template: its
     * contents are spliced into the parent, so nothing is silently lost. */
    if (path_.empty()) return false;
    CaretStep st = path_.back();
    std::vector<CaretStep> up(path_.begin(), path_.end() - 1);
    NodeList* parent_list = slot_at(up);
    if (!parent_list || size_t(st.child) >= parent_list->size()) return false;

    checkpoint("Delete");
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
    /* After the wrapper is gone the caret belongs at the END of what was
     * inside it, not before.  Sitting in front meant the next Backspace had
     * nothing to its left and did nothing at all: repeating the key stopped
     * partway through instead of peeling the rest of the equation. */
    index_ = at;
    clamp();
    refresh_auto_function_word();
    return true;
}

bool Equation::erase() {
    if (has_selection()) return delete_selection();
    NodeList& l = here();
    clamp();
    if (size_t(index_) >= l.size()) {
        /* Nothing ahead in this slot.  If the slot is the base of a script
         * whose sub/superscript are empty, those empty boxes are what sits
         * ahead of the caret, so forward-delete removes them -- the mirror of
         * Backspace unwrapping a template from the start of an empty slot.
         * Without this a fresh sub/sup applied to a whole selection could not
         * be taken off again while keeping the base: the caret sat at the end
         * of the base, Backspace ate the base content, and Delete did nothing.
         * Reported for `{...}_{}^{}` with the caret at the end of the base. */
        if (path_.empty()) return false;
        const CaretStep st = path_.back();
        Node* owner = parent_node(path_);
        if (st.slot != 0 || !owner || owner->tag() != Node::kScript) {
            /* Nothing ahead in this slot and no script rule applies.  Step
             * back out to just before the template and let the caller press
             * again: the next Delete re-enters it and takes the first item of
             * whatever slot still holds something.  Without this, Delete
             * stopped dead once the caret had followed a deletion into a slot
             * and emptied it -- the same early stop Backspace had. */
            path_.pop_back();
            index_ = st.child;
            clamp();
            return true;
        }
        Node* parent = owner;
        auto& script = static_cast<ScriptNode&>(*parent);
        const bool emptySub = script.hasSub && script.sub.empty();
        const bool emptySup = script.hasSup && script.sup.empty();
        if (!emptySub && !emptySup) {
            /* The base is exhausted but the scripts still hold something.
             * Step out so the next press re-enters and takes the first of
             * them, rather than reporting that there is nothing left. */
            path_.pop_back();
            index_ = st.child;
            clamp();
            return true;
        }

        checkpoint("Delete");
        parent = parent_node(path_);
        auto& editable = static_cast<ScriptNode&>(*parent);
        if (editable.hasSub && editable.sub.empty()) editable.hasSub = false;
        if (editable.hasSup && editable.sup.empty()) editable.hasSup = false;
        if (!editable.hasSub && !editable.hasSup) {
            /* No scripts left: unwrap, splicing the base into the parent and
             * leaving the caret after that content. */
            std::vector<CaretStep> up(path_.begin(), path_.end() - 1);
            NodeList* grand = slot_at(up);
            if (grand && size_t(st.child) < grand->size()) {
                NodeList base = std::move(editable.base);
                grand->erase(grand->begin() + st.child);
                int at = st.child;
                for (auto& n : base)
                    grand->insert(grand->begin() + at++, std::move(n));
                path_ = up;
                index_ = at;
            }
        }
        clamp();
        refresh_auto_function_word();
        return true;
    }
    if (l[size_t(index_)] && l[size_t(index_)]->tag() == Node::kScript) {
        auto& script = static_cast<ScriptNode&>(*l[size_t(index_)]);
        const bool emptySub = script.hasSub && script.sub.empty();
        const bool emptySup = script.hasSup && script.sup.empty();
        if (emptySub || emptySup) {
            checkpoint("Delete");
            NodeList& cur = here();
            auto& editable = static_cast<ScriptNode&>(*cur[size_t(index_)]);
            if (editable.hasSub && editable.sub.empty()) editable.hasSub = false;
            if (editable.hasSup && editable.sup.empty()) editable.hasSup = false;
            if (!editable.hasSub && !editable.hasSup) {
                NodeList base = std::move(editable.base);
                cur.erase(cur.begin() + index_);
                size_t at = size_t(index_);
                for (auto& node : base)
                    cur.insert(cur.begin() + at++, std::move(node));
            }
            refresh_auto_function_word();
            return true;
        }
    }
    /* Same rule as Backspace, read the other way: one press removes one thing
     * the reader can see, from the first slot that still holds something --
     * and the caret follows it in, so the next character typed lands where
     * the deletion happened.  Left outside, Delete on x_1^2 emptied the base
     * and then typing put the character in front of the whole script, which
     * emitted a stray empty group. */
    if (l[size_t(index_)] && first_nonempty_slot(*l[size_t(index_)])) {
        checkpoint("Delete");
        NodeList& fresh = here();
        std::vector<CaretStep> steps;
        if (NodeList* slot = first_nonempty_path(fresh, index_, steps)) {
            slot->erase(slot->begin());
            for (const CaretStep& step : steps) path_.push_back(step);
            index_ = 0;
        }
        refresh_auto_function_word();
        return true;
    }

    int functionPrefixEnd = -1;
    if (l[size_t(index_)] && l[size_t(index_)]->tag() == Node::kChar &&
        !static_cast<const CharNode*>(l[size_t(index_)].get())
             ->automaticFunction)
        functionPrefixEnd = automatic_function_prefix_end(l, index_);
    checkpoint("Delete");
    NodeList& cur = here();
    cur.erase(cur.begin() + index_);
    refresh_auto_function_word(functionPrefixEnd);
    return true;
}

/* ------------------------------------------------------------ commands */

const std::vector<std::string>& Equation::templates() {
    /* Every entry has a home on one of the palettes (see palettes.cpp), so
     * nothing here is reachable only by typing its command. */
    static const std::vector<std::string> kAll = {
        "frac", "slashfrac", "sqrt", "nthroot",
        "sub", "sup", "subsup", "over", "under", "lim",
        "paren", "bracket", "brace", "abs", "angle",
        "norm", "floor", "ceil", "dirac",
        "int", "iint", "iiint", "oint",
        "sum", "prod", "coprod", "bigcup", "bigcap",
        "overline", "underline", "overbrace", "underbrace",
        "overrightarrow", "overleftarrow", "overleftrightarrow",
        "hat", "vec", "bar", "tilde", "dot", "ddot", "dddot",
        "prime", "dprime", "tprime", "strike", "frown", "smile",
        "matrix1x2", "matrix1x3", "matrix2x1", "matrix2x2", "matrix2x3",
        "matrix3x1", "matrix3x2", "matrix3x3", "matrix4x4", "matrix5x5",
        "matrix6x6", "cases",
    };
    return kAll;
}

const std::vector<Equation::SequenceBinding>& Equation::sequence_shortcuts() {
    static const std::vector<SequenceBinding> kSequences = [] {
        std::vector<SequenceBinding> result = {
        /* Eqnedit32's template family. */
        {'T', 'F', false, "template.frac",       "fraction"},
        {'T', '/', false, "template.slashfrac",  "slashed fraction"},
        {'T', 'H', false, "template.sup",        "superscript"},
        {'T', 'L', false, "template.sub",        "subscript"},
        {'T', 'J', false, "template.subsup",     "subscript and superscript"},
        {'T', 'I', false, "template.int",        "integral"},
        {'T', '|', false, "template.abs",        "absolute value"},
        {'T', 'R', false, "template.sqrt",       "square root"},
        {'T', 'N', false, "template.nthroot",    "n-th root"},
        {'T', 'S', false, "template.sum",        "summation"},
        {'T', 'P', false, "template.prod",       "product"},
        {'T', 'M', false, "template.matrix3x3",  "3 by 3 matrix"},
        {'T', 'U', false, "template.under",      "annotation underneath"},

        /* Eqnedit64 additions that do not replace an old mnemonic. */
        {'T', 'D', false, "template.iint",       "double integral"},
        {'T', 'O', false, "template.oint",       "contour integral"},
        {'T', '2', false, "template.matrix2x2",  "2 by 2 matrix"},
        {'T', '3', false, "template.matrix3x3",  "3 by 3 matrix"},
        {'T', 'C', false, "template.cases",      "cases"},
        {'T', 'A', false, "template.abs",        "absolute value"},
        {'T', 'B', false, "template.bracket",    "brackets"},
        {'T', 'G', false, "template.brace",      "braces"},
        {'T', 'E', false, "template.angle",      "angle brackets"},
        {'T', 'V', false, "template.vec",        "vector accent"},
        {'T', 'Q', false, "template.hat",        "hat accent"},
        {'T', '-', false, "template.overline",   "overline"},
        {'T', '_', false, "template.underline",  "underline"},
        {'T', 'W', false, "template.tilde",      "tilde accent"},

        /* Eqnedit32's Ctrl+K symbol family. */
        {'K', 'I', false, "symbol.\\infty",      "infinity"},
        {'K', 'A', false, "symbol.\\rightarrow", "right arrow"},
        {'K', 'D', false, "symbol.\\partial",    "partial derivative"},
        {'K', '<', false, "symbol.\\leq",        "less than or equal"},
        {'K', '>', false, "symbol.\\geq",        "greater than or equal"},
        {'K', 'T', false, "symbol.\\times",      "multiplication"},
        {'K', 'E', false, "symbol.\\in",         "element of"},
        {'K', 'E', true,  "symbol.\\notin",      "not an element of"},
        {'K', 'C', false, "symbol.\\subset",     "subset"},
        {'K', 'C', true,  "symbol.\\not\\subset", "not a subset"},

        {'G', 'A', false, "symbol.\\alpha",   "alpha"},
        {'G', 'B', false, "symbol.\\beta",    "beta"},
        {'G', 'G', false, "symbol.\\gamma",   "gamma"},
        {'G', 'D', false, "symbol.\\delta",   "delta"},
        {'G', 'E', false, "symbol.\\epsilon", "epsilon"},
        {'G', 'Z', false, "symbol.\\zeta",    "zeta"},
        {'G', 'H', false, "symbol.\\eta",     "eta"},
        {'G', 'T', false, "symbol.\\theta",   "theta"},
        {'G', 'I', false, "symbol.\\iota",    "iota"},
        {'G', 'K', false, "symbol.\\kappa",   "kappa"},
        {'G', 'L', false, "symbol.\\lambda",  "lambda"},
        {'G', 'M', false, "symbol.\\mu",      "mu"},
        {'G', 'N', false, "symbol.\\nu",      "nu"},
        {'G', 'X', false, "symbol.\\xi",      "xi"},
        {'G', 'P', false, "symbol.\\pi",      "pi"},
        {'G', 'R', false, "symbol.\\rho",     "rho"},
        {'G', 'S', false, "symbol.\\sigma",   "sigma"},
        {'G', 'U', false, "symbol.\\upsilon", "upsilon"},
        {'G', 'F', false, "symbol.\\phi",     "phi"},
        {'G', 'C', false, "symbol.\\chi",     "chi"},
        {'G', 'Y', false, "symbol.\\psi",     "psi"},
        {'G', 'O', false, "symbol.\\omega",   "omega"},

        {'G', 'G', true, "symbol.\\Gamma",   "capital Gamma"},
        {'G', 'D', true, "symbol.\\Delta",   "capital Delta"},
        {'G', 'T', true, "symbol.\\Theta",   "capital Theta"},
        {'G', 'L', true, "symbol.\\Lambda",  "capital Lambda"},
        {'G', 'X', true, "symbol.\\Xi",      "capital Xi"},
        {'G', 'P', true, "symbol.\\Pi",      "capital Pi"},
        {'G', 'S', true, "symbol.\\Sigma",   "capital Sigma"},
        {'G', 'U', true, "symbol.\\Upsilon", "capital Upsilon"},
        {'G', 'F', true, "symbol.\\Phi",     "capital Phi"},
        {'G', 'Y', true, "symbol.\\Psi",     "capital Psi"},
        {'G', 'O', true, "symbol.\\Omega",   "capital Omega"},
        };

        /* Ctrl+B makes the next character a matrix/vector character, as in
         * Eqnedit32.  Keeping it in the same table makes it testable and lets
         * the shortcut coach teach it. */
        for (char key = 'A'; key <= 'Z'; ++key) {
            const char lower = static_cast<char>(key - 'A' + 'a');
            result.push_back({'B', key, false,
                std::string("latex.\\mathbf{") + lower + "}",
                std::string("vector ") + lower});
            result.push_back({'B', key, true,
                std::string("latex.\\mathbf{") + key + "}",
                std::string("vector ") + key});
        }
        return result;
    }();
    return kSequences;
}

const Equation::SequenceBinding* Equation::resolve_sequence(
        char prefix, char key, bool shift) {
    for (const auto& binding : sequence_shortcuts())
        if (binding.prefix == prefix && binding.key == key &&
            binding.shift == shift)
            return &binding;
    return nullptr;
}

/* Direct chords plus discoverable two-stroke families.  Prefix sequences are
 * generated from the same table used by the GUI so Help, tests, and actual
 * key handling cannot drift apart. */
const std::vector<Equation::Binding>& Equation::shortcuts() {
    static const std::vector<Binding> kBindings = [] {
        std::vector<Binding> result = {
            {"Ctrl+F",       "template.frac",      "fraction"},
            {"Ctrl+/",       "template.slashfrac", "slashed fraction"},
            {"Ctrl+R",       "template.sqrt",      "square root"},
            {"Ctrl+Alt+R",   "template.nthroot",   "n-th root"},
            {"Ctrl+H",       "template.sup",       "superscript (high)"},
            {"Ctrl+L",       "template.sub",       "subscript (low)"},
            {"Ctrl+J",       "template.subsup",    "subscript and superscript"},
            {"Ctrl+9",       "template.paren",     "parentheses"},
            {"Ctrl+0",       "template.paren",     "parentheses"},
            {"Ctrl+[",       "template.bracket",   "brackets"},
            {"Ctrl+]",       "template.bracket",   "brackets"},
            {"Ctrl+{",       "template.brace",     "braces"},
            {"Ctrl+}",       "template.brace",     "braces"},
            {"Ctrl+I",       "template.int",       "integral"},
            {"Ctrl+Alt+I",   "template.iint",      "double integral"},
            {"Ctrl+Shift+O", "template.oint",      "contour integral"},
            {"Ctrl+M",       "template.matrix2x2", "2 by 2 matrix"},
            {"Ctrl+Shift+M", "template.matrix3x3", "3 by 3 matrix"},
            {"Ctrl+Alt+Down", "matrix.add_row",     "add matrix row"},
            {"Ctrl+Alt+Up",   "matrix.remove_row",  "remove matrix row"},
            {"Ctrl+Alt+Right","matrix.add_column",  "add matrix column"},
            {"Ctrl+Alt+Left", "matrix.remove_column","remove matrix column"},
            {"Ctrl+-",       "template.bar",       "overbar"},
            {"Ctrl+Shift+~", "template.tilde",     "tilde accent"},
            {"Ctrl+Alt+-",   "template.vec",       "vector accent"},
            {"Ctrl+Alt+.",   "template.dot",       "dot accent"},
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
        for (const auto& sequence : sequence_shortcuts()) {
            std::string chord = "Ctrl+";
            chord.push_back(sequence.prefix);
            chord += ", ";
            if (sequence.shift) chord += "Shift+";
            chord.push_back(sequence.key);
            result.push_back({chord, sequence.command, sequence.label});
        }
        return result;
    }();
    return kBindings;
}

bool Equation::command(const std::string& name) {
    if (name.compare(0, 9, "template.") == 0)
        return insert_template(name.substr(9));
    if (name.compare(0, 7, "symbol.") == 0)
        return insert_symbol(name.substr(7));
    if (name.compare(0, 6, "latex.") == 0)
        return insert_latex(name.substr(6));
    if (name == "caret.left")       return move_left();
    if (name == "caret.right")      return move_right();
    if (name == "caret.next_slot")  return next_slot();
    if (name == "caret.prev_slot")  return prev_slot();
    if (name == "caret.out")        return move_out();
    if (name == "caret.home")       { move_home(); return true; }
    if (name == "caret.end")        { move_end();  return true; }
    if (name == "caret.up")         return move_up();
    if (name == "caret.down")       return move_down();
    if (name == "edit.new_line")    return new_line();
    if (name == "edit.alignment")   return alignment_tab();
    if (name == "edit.backspace")   return backspace();
    if (name == "edit.delete")      return erase();
    if (name == "edit.undo")        return undo();
    if (name == "edit.redo")        return redo();
    if (name == "edit.select_all")  return select_all();
    if (name == "matrix.add_row")       return matrix_add_row();
    if (name == "matrix.remove_row")    return matrix_remove_row();
    if (name == "matrix.add_column")    return matrix_add_column();
    if (name == "matrix.remove_column") return matrix_remove_column();
    return false;
}

}  // namespace eqnedit
