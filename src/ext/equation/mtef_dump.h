/*
 * mtef_dump.h -- indented text dump of the parsed node tree (diagnostic).
 *
 * Used to settle which field of the tree actually carries a piece of content
 * when the LaTeX, SVG and OMML emitters disagree.  `run_passes` selects the
 * raw parse or the post-pass tree the emitters actually see.
 */
#ifndef MTEF_DUMP_H
#define MTEF_DUMP_H

#include "mtef_node.h"

#include <cstddef>
#include <cstdint>
#include <string>

namespace mtef {

std::string dump_tree(const uint8_t* data, size_t len, bool run_passes = true);

/* Same dump for a tree that is already built (e.g. from the LaTeX parser).
 * No passes are run: those repair MTEF's layout and would corrupt this tree. */
std::string dump_latex_tree(const LineNode& root);

}  // namespace mtef

#endif /* MTEF_DUMP_H */
