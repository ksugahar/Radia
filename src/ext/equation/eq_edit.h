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
#include "math_layout.h"
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
    bool operator==(const CaretStep& o) const {
        return child == o.child && slot == o.slot;
    }
    bool operator!=(const CaretStep& o) const { return !(*this == o); }
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
    /* One whole item, without stepping inside it: Ctrl+Left and Ctrl+Right.
     * `extend` makes it a selection instead. */
    bool move_item(int dir, bool extend);
    /* Enter: break the line here, making a stack of two if there is not one
     * already. */
    bool newline();
    /* One slot of the enclosing node: Ctrl+Up / Ctrl+Down. */
    bool move_slot(int dir);
    void move_home();
    void move_end();

    /* ---- style --------------------------------------------------------- */
    /* Equation Editor's Style menu, applied to the selection.
     *
     * The tree has carried these typefaces all along and nothing could set
     * them, so a bold vector -- most of what this lab writes -- was
     * unreachable from the editor.
     *
     * "greek" is not a typeface but a keyboard: Equation Editor's Greek style
     * maps the Latin letters onto Greek ones, so `a` becomes alpha.  Treating
     * it as a font would leave the letter Latin and merely restyle it.
     *
     * Names: math, text, function, variable, vector, greek. */
    /* Equation Editor's Style menu is a MODE, not an operation on what is
     * highlighted: its status bar reads "Style: Math" and stays that way, so
     * picking Greek and carrying on typing gives Greek.  With a selection this
     * also restyles it, which is the other half of what the menu does. */
    bool set_style(const std::string& name);
    const std::string& style() const { return style_; }
    static const std::vector<std::string>& styles();

    /* ---- selection ----------------------------------------------------- */
    /* A range within ONE slot.
     *
     * A template is a single item in its parent's slot, so selecting a whole
     * fraction is the same operation as selecting a run of characters -- which
     * means the useful cases are covered without a model that spans levels of
     * the tree.  Shift at a slot boundary stops there rather than jumping out,
     * because a selection that silently changes what it is anchored to is
     * worse than one that will not grow.
     *
     * The anchor is where the selection started; the caret is its other end.
     * They may be in either order. */
    bool has_selection() const;
    void clear_selection();
    void select_all();                  /* everything in the current slot */

    bool extend_left();                 /* Shift+Left      */
    bool extend_right();                /* Shift+Right     */
    void extend_home();                 /* Shift+Home      */
    void extend_end();                  /* Shift+End       */
    bool extend_to_point(double x, double y, const SvgStyle& style = SvgStyle());

    /* The selected nodes as LaTeX.  Not const: it borrows the nodes out of the
     * tree to write them and puts them back, which is cheaper and less
     * error-prone than teaching every node type to clone itself. */
    std::string selected_latex();
    bool delete_selection();

    /* What to highlight.  `found` is false when nothing is selected. */
    struct SelectionBox {
        bool found = false;
        double x0 = 0, x1 = 0, top = 0, bottom = 0;
    };
    SelectionBox selection_geometry(const SvgStyle& style = SvgStyle()) const;

    /* ---- editing ----------------------------------------------------- */
    /* Every one of these replaces the selection when there is one, which is
     * what typing over a highlighted range has to do. */
    void insert_text(const std::string& utf8);      /* literal characters */
    bool insert_symbol(const std::string& cmd);     /* "\\alpha", "\\nabla" */

    /* Paste: parse LaTeX and put it in at the caret.
     *
     * Goes through the model rather than replacing the document, so it lands
     * inside whatever slot the caret is in, replaces a selection, and undoes
     * in one step.  Surrounding $...$ are stripped by the parser, so an
     * equation copied out of a Markdown file arrives as an equation. */
    bool insert_latex(const std::string& latex);
    bool insert_template(const std::string& kind);  /* "frac", "sqrt", ... */
    bool backspace();
    bool erase();                                    /* forward delete */

    /* ---- geometry ---------------------------------------------------- */
    /* Where the caret is, in points relative to the equation's own origin
     * (on the baseline at its left edge).  `found` is false when the layout
     * has no position for the current caret, which happens inside a construct
     * the layout does not draw yet -- an editor should leave the caret where it
     * was rather than pretend. */
    struct CaretGeometry {
        bool found = false;
        double x = 0, top = 0, bottom = 0;
    };
    CaretGeometry caret_geometry(const SvgStyle& style = SvgStyle()) const;

    /* Put the caret where a click landed.  Coordinates are relative to the same
     * origin.  False when the equation is empty. */
    bool move_to_point(double x, double y, const SvgStyle& style = SvgStyle());

    /* Up (-1) and down (+1).  Structural motion cannot answer this -- a
     * numerator is not above a denominator in the tree, only on the page --
     * so it is decided from the laid-out geometry. */
    bool move_vertical(int dir, const SvgStyle& style = SvgStyle());

    /* The equation's own box, so a caller can place it and convert a click. */
    void extents(double& w, double& asc, double& desc,
                 const SvgStyle& style = SvgStyle()) const;

    /* The display list, for a window that draws the equation being edited.
     * The same layout the caret geometry and the pictures come from, so the
     * editor cannot show something the paste does not contain.
     *
     * Recomputed on every call rather than cached and patched.  Equation
     * Editor 3.0 needed a Redraw command because its display could go stale;
     * a layout derived fresh from the tree has nothing to go stale. */
    Layout layout(const SvgStyle& style = SvgStyle()) const;

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

    /* The palette bar, following Equation Editor 3.0's arrangement: a row of
     * symbol palettes above a row of template palettes.
     *
     * The contents are BUILT from the shared command table and the template
     * list, never written out again here.  A palette listed by hand drifts:
     * it starts offering symbols the editor cannot insert, and stops showing
     * ones that were added.  This is the same reason the key chords come from
     * shortcuts() rather than from the window.
     *
     * A group with nothing in it is not returned at all -- an empty dropdown
     * is worse than an absent one. */
    struct PaletteItem {
        std::string command;      /* "\\alpha", or a template name    */
        uint32_t code = 0;        /* the glyph to show on the button  */
        bool is_template = false;
    };
    struct PaletteGroup {
        std::string name;
        std::vector<PaletteItem> items;
    };
    static const std::vector<PaletteGroup>& symbol_palettes();
    static const std::vector<PaletteGroup>& template_palettes();

private:
    std::unique_ptr<LineNode> root_;
    std::vector<CaretStep> path_;
    int index_ = 0;
    int anchor_ = -1;              /* where a selection started; -1 for none */
    /* WHICH SLOT the anchor is an index into.
     *
     * Without this the anchor was a bare number, and every operation that
     * moved to another slot or rebuilt the tree -- Tab, Ctrl+Up, undo -- left
     * it pointing into a slot it no longer belonged to.  A shorter slot then
     * made take_selection() read past the end of the vector and erase a range
     * that was not there: a heap corruption that killed the editor some
     * seconds later, somewhere else entirely.  It is what the window
     * self-test found first.
     *
     * A selection belongs to the slot it was made in.  Leave that slot, or
     * change the tree under it, and there is NO selection -- not a selection
     * quietly shrunk to fit, which would delete a different range than the
     * one the user highlighted. */
    std::vector<CaretStep> anchor_path_;
    std::string style_ = "math";   /* the style typing goes in */

public:
    /* A chord half-entered.  It lives here rather than beside the window
     * because a two-step chord is part of what the EQUATION is in the middle
     * of, and because a table keyed on the object's address hands one
     * equation's pending prefix to whatever is allocated there next -- which
     * is exactly what happened. */
    struct PendingChord { std::vector<int> vk, mods; };
    PendingChord& pending_chord() { return pending_; }
private:
    PendingChord pending_;

    struct Snapshot { std::string latex; std::string caret; };
    std::vector<Snapshot> undo_, redo_;

    NodeList* slot_at(const std::vector<CaretStep>& path) const;
    NodeList& here();
    Node* parent_node(const std::vector<CaretStep>& path) const;

    /* Start a selection here, remembering which slot "here" is. */
    void set_anchor(int at);
    /* Begin one at the caret if there is not one already -- what every
     * Shift+move does first. */
    void start_selection_here();
    /* The selection as a range in the CURRENT slot, validated against it.
     * False when there is none -- including when the anchor was made in
     * another slot, or no longer fits the one it was made in.  This is the
     * ONE place lo/hi are computed, so no caller can index past the end. */
    bool selection_range(int& lo, int& hi) const;

    /* Remove the selected nodes without pushing an undo step, so replacing a
     * selection is ONE undo rather than two.  Returns them, because a template
     * wraps what was selected instead of discarding it. */
    NodeList take_selection();

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
