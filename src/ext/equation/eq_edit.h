/*
 * eq_edit.h -- the equation editing model
 *
 * What makes Equation Editor 3.0 pleasant is not its renderer, it is the way
 * typing works: an insertion point that lives *inside* the structure, templates
 * whose empty slots you tab through, and backspace that unwraps a template
 * instead of deleting a character you cannot see.  That behaviour is a data
 * structure plus a set of operations, entirely separate from any window
 * toolkit, so it lives here and is testable headlessly.
 *
 * The stored form of an equation is its LaTeX.  The tree is the working form:
 * loaded from LaTeX, mutated by these operations, written back as LaTeX.
 *
 * Undo snapshots the LaTeX rather than inverting each command.  Equations are
 * tens of characters, so a full snapshot costs nothing and cannot drift out of
 * step with the tree the way command inversion can.  It relies on
 * LaTeX -> tree -> LaTeX reaching a fixed point, which is checked by a test.
 */
#ifndef EQ_EDIT_H
#define EQ_EDIT_H

#include "mtef_node.h"
#include "mtef_omml.h"
#include "mtef_svg.h"

#include <memory>
#include <string>
#include <vector>

namespace mtef {

/* One step down the tree: into child `child` of the current slot, then into
 * that child's slot `slot`. */
struct CaretStep {
    int child;
    int slot;
};

class Equation {
public:
    Equation();

    /* ---- document ---------------------------------------------------- */
    void load_latex(const std::string& latex);   /* clears history */
    std::string latex() const;
    std::string omml(const OmmlOptions& opt = OmmlOptions()) const;
    std::string svg(const SvgStyle& style = SvgStyle()) const;

    /* ---- caret ------------------------------------------------------- */
    /* Path/index as text, e.g. "1.0/0:2" -- child 1 slot 0, then index 2.
     * Used by tests and by a status line; the GUI keeps the caret opaque. */
    std::string caret() const;
    bool move_left();
    bool move_right();
    bool next_slot();      /* Tab      */
    bool prev_slot();      /* Shift+Tab */
    bool move_out();       /* leave the innermost slot */
    void move_home();
    void move_end();

    /* ---- editing ----------------------------------------------------- */
    void insert_text(const std::string& utf8);      /* literal characters */
    bool insert_symbol(const std::string& cmd);     /* "\\alpha", "\\nabla" */
    bool insert_template(const std::string& kind);  /* "frac", "sqrt", ... */
    bool backspace();
    bool erase();                                    /* forward delete */

    /* ---- history ----------------------------------------------------- */
    bool undo();
    bool redo();

    /* ---- command dispatch -------------------------------------------- */
    /* One name per editing action, so a GUI binds keys to names and the model
     * never knows about keyboards.  Unknown names return false. */
    bool command(const std::string& name);

    /* The Equation Editor 3.0 key chords, as name/chord pairs.  Kept here so
     * the familiar chords survive into the new editor rather than being
     * reinvented per front end. */
    struct Binding { const char* chord; const char* command; const char* label; };
    static const std::vector<Binding>& shortcuts();

    /* Template names this build understands, for a palette or a test. */
    static const std::vector<std::string>& templates();

private:
    std::unique_ptr<LineNode> root_;
    std::vector<CaretStep> path_;
    int index_ = 0;

    struct Snapshot { std::string latex; std::string caret; };
    std::vector<Snapshot> undo_, redo_;

    NodeList* slot_at(const std::vector<CaretStep>& path) const;
    NodeList& here();
    Node* parent_node(const std::vector<CaretStep>& path) const;

    void checkpoint();                  /* push the current state for undo */
    void restore(const Snapshot& s);
    std::string caret_text() const;
    void set_caret_text(const std::string& s);
    void enter_first_empty_slot(Node& n);
    void clamp();
};

/* Slots of a node in visual order.  Empty for a leaf. */
std::vector<NodeList*> node_slots(Node& n);

}  // namespace mtef

#endif /* EQ_EDIT_H */
