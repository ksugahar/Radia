/* Toolbar palettes: the catalogue of everything that can be inserted.
 *
 * Eqnedt32 puts every symbol and every template on one of nineteen toolbar
 * palettes -- ten of symbols, nine of templates -- and has no Insert menu at
 * all.  Two clicks reach any of its ~300 items.  Eqnedit64 had flat menus
 * instead, which left 117 of its 160 symbols with no mouse route and no
 * shortcut; the only way in was to know and type the command.  This table is
 * the same catalogue idea, and tests/test_palettes.py is what keeps it
 * honest: every symbol and every template must have exactly one home here.
 */
#ifndef EQNEDIT_PALETTES_H
#define EQNEDIT_PALETTES_H

#include <cstddef>
#include <string>
#include <vector>

namespace eqnedit {

struct PaletteItem {
    /* "template.<kind>", "symbol.<command>", "latex.<source>", or a
     * contextual "matrix.<operation>" -- the same command vocabulary the
     * keyboard shortcuts use. */
    std::string command;
    std::string face;    /* what the cell shows */
    std::string label;   /* what the status bar explains */
};

struct Palette {
    std::string title;   /* palette name, shown as the popup heading */
    std::string face;    /* toolbar button face */
    int columns = 3;     /* grid width, as in Eqnedt32's drop-downs */
    std::vector<PaletteItem> items;
};

/* Five user-facing groups shared by the menu bar and the compact palette
 * tabs.  Every palette index belongs to exactly one group. */
struct PaletteCategory {
    std::string title;
    std::vector<std::size_t> paletteIndices;
};

/* Nineteen palettes: the first eleven are symbols, the last eight templates. */
const std::vector<Palette>& palettes();

/* Number of leading entries of palettes() that are symbol palettes. */
int symbol_palette_count();

/* Basic / analysis / sets-symbols / geometry / Greek, in display order. */
const std::vector<PaletteCategory>& palette_categories();

}  // namespace eqnedit

#endif
