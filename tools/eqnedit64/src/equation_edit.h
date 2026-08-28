/*
 * equation_edit.h -- the equation editing model
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
#ifndef EQUATION_EDIT_H
#define EQUATION_EDIT_H

#include "equation_node.h"
#include "equation_render.h"

#include <memory>
#include <string>
#include <vector>

namespace eqnedit {

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
    bool load_latex(const std::string& latex);   /* clears history */
    /* Replace the equation from TeX that was edited in place, keeping the
     * history.  load_latex() is a document load and drops Undo, so routing
     * raw-TeX-pane keystrokes through it threw away everything the canvas
     * had done.  checkpointFirst is true only for the first keystroke of an
     * editing burst, so one Undo steps back to the equation as it stood
     * before the pane was touched rather than one character at a time. */
    bool replace_latex(const std::string& latex, bool checkpointFirst);
    std::string latex() const;
    std::string svg(const SvgStyle& style = SvgStyle()) const;
    /* Empty after a successful operation.  A rejected edit leaves the model,
     * caret, selection, and Undo/Redo stacks unchanged. */
    const std::string& last_error() const { return lastError_; }

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

    /* Mouse/selection support for the native 64-bit editor.  Selection is a
     * contiguous range in one structural slot (the common case for visual
     * drag selection); Ctrl+A selects the root equation. */
    bool hit_test(double x_points, double y_points,
                  const SvgStyle& style = SvgStyle(), bool extend = false);
    bool new_line();          /* Enter: create/append an aligned row */
    bool alignment_tab();     /* &: split the current aligned row */
    bool move_up();
    bool move_down();
    bool is_multiline() const;
    /* Ordinary matrix editing.  Insertion accepts practical arbitrary sizes
     * (1..99 in each direction); add/remove operates on the row or column
     * containing the caret and is a single named Undo step. */
    bool matrix_dimensions(int* rows, int* cols) const;
    bool matrix_add_row();
    bool matrix_remove_row();
    bool matrix_add_column();
    bool matrix_remove_column();
    void begin_selection();
    void clear_selection();
    /* Extend the selection by one whole atom within the current slot.  Unlike
     * move_left/move_right, these never descend into a structure: Shift+arrow
     * selects the script or fraction next to the caret as a unit, and once the
     * caret is inside a slot (a superscript, say) it selects the content
     * there.  Descending mid-drag was why Shift+arrow silently selected
     * nothing when it crossed a structural boundary. */
    bool select_step_left();
    bool select_step_right();
    /* Select the contents of the structural slot containing the caret.
     * On the outermost slot this is equivalent to select_all().  This is the
     * legacy editor's double-click rule and is more useful than selecting the
     * entire equation from every nested slot. */
    bool select_current_slot();
    /* Select the innermost structural node containing the caret as one atom
     * in its parent slot.  This is the safe TeX-tree meaning of the legacy
     * Ctrl+click on a template symbol: a fraction, radical, fence, script, or
     * matrix is cut/replaced together with its contents, never torn apart. */
    bool select_containing_structure();
    bool select_all();
    bool has_selection() const;
    std::string selection_latex() const;
    bool delete_selection();

    /* ---- editing ----------------------------------------------------- */
    void insert_text(const std::string& utf8);      /* literal characters */
    bool insert_styled_text(const std::string& utf8,
                            const std::string& style); /* text/function/vector/variable */
    /* Restyle every character in the current selection -- Ctrl+B over a
     * selected symbol makes it a vector, rather than only setting the style
     * for the next character typed.  No selection: returns false. */
    bool restyle_selection(const std::string& style);
    bool insert_latex(const std::string& latex);    /* structural paste */
    bool insert_symbol(const std::string& cmd);     /* "\\alpha", "\\nabla" */
    bool insert_template(const std::string& kind);  /* "frac", "sqrt", ... */
    bool backspace();
    bool erase();                                    /* forward delete */

    /* ---- history ----------------------------------------------------- */
    bool can_undo() const;
    bool can_redo() const;
    std::string undo_name() const;
    std::string redo_name() const;
    bool undo();
    bool redo();

    /* ---- command dispatch -------------------------------------------- */
    /* One name per editing action, so a GUI binds keys to names and the model
     * never knows about keyboards.  Unknown names return false. */
    bool command(const std::string& name);

    RenderMetrics metrics(const SvgStyle& style = SvgStyle()) const;
    bool caret_geometry(CaretGeometry* geometry,
                        const SvgStyle& style = SvgStyle()) const;
#ifdef _WIN32
    void draw_gdi(HDC hdc, double left, double top, double scale,
                  const SvgStyle& style = SvgStyle(),
                  bool show_placeholders = true,
                  bool show_caret = true,
                  bool text_as_outlines = false) const;
#endif

    /* The Equation Editor 3.0 key chords, as name/chord pairs.  Kept here so
     * the familiar chords survive into the new editor rather than being
     * reinvented per front end. */
    struct Binding {
        std::string chord;
        std::string command;
        std::string label;
    };
    struct SequenceBinding {
        char prefix;
        char key;
        bool shift;
        std::string command;
        std::string label;
    };
    static const std::vector<Binding>& shortcuts();
    static const std::vector<SequenceBinding>& sequence_shortcuts();
    static const SequenceBinding* resolve_sequence(char prefix, char key,
                                                   bool shift = false);

    /* Template names this build understands, for a palette or a test. */
    /* Palette template names.  insert_template also accepts matrixRxC for
     * every 1 <= R,C <= 99, so arbitrary sizes need no giant catalogue. */
    static const std::vector<std::string>& templates();

private:
    std::unique_ptr<LineNode> root_;
    std::vector<CaretStep> path_;
    int index_ = 0;
    bool selecting_ = false;
    std::vector<CaretStep> anchorPath_;
    int anchorIndex_ = 0;

    struct Snapshot {
        std::string latex;
        std::string caret;
        std::string action;
    };
    std::vector<Snapshot> undo_, redo_;
    std::string lastError_;

    NodeList* slot_at(const std::vector<CaretStep>& path) const;
    NodeList& here();
    Node* parent_node(const std::vector<CaretStep>& path) const;

    void checkpoint(const char* action); /* push current state for named undo */
    void restore(const Snapshot& s);
    std::string caret_text() const;
    void set_caret_text(const std::string& s);
    void enter_first_empty_slot(Node& n);
    void clamp();
    void insert_text_typeface(const std::string& utf8, int forcedTypeface);
    void refresh_auto_function_word(int leftBoundary = -1);
    bool find_slot_path(const NodeList* target, NodeList& list,
                        std::vector<CaretStep>& path) const;
    bool selection_range(NodeList** slot, int* first, int* last) const;
    bool matrix_context(size_t* depth, MatrixNode** matrix, int* cell,
                        int layoutKind = -1) const;
};

/* Slots of a node in visual order.  Empty for a leaf. */
std::vector<NodeList*> node_slots(Node& n);

}  // namespace eqnedit

#endif /* EQUATION_EDIT_H */
