/*
 * line_pass.h -- Pass pipeline for LINE node children processing
 *
 * Each LinePass operates on a LineNode's children with a skip bitset.
 * Passes may:
 *   - Mark children as skipped (skip[i] = true)
 *   - Mutate node fields (e.g., fill fence content from siblings)
 *   - Set display_lower/display_upper on integral/bigop nodes
 *
 * Passes are executed in order: 0a → 0 → 1 → 2 → 2b → 2c → 3
 */
#ifndef LINE_PASS_H
#define LINE_PASS_H

#include "mtef_node.h"
#include <vector>

namespace mtef {

/* Skip bitset for LINE children */
using SkipSet = std::vector<bool>;

/* ============================================================
 * Base class
 * ============================================================ */
class LinePass {
public:
    virtual ~LinePass() = default;

    /* Process LINE children. May modify nodes and skipSet.
     * depth: nesting depth (0 = root LINE). */
    virtual void run(NodeList& children, SkipSet& skip, int depth, int prodVer) = 0;
};

/* ============================================================
 * Pass 0a: Display-fraction reassembly
 *
 * EQNEDT32 writes a DISPLAY fraction across the template boundary: slot[1]
 * holds only the first line of the numerator, and the rest of it -- and the
 * denominator -- follow as siblings, separated by full-size records.  Read
 * straight that gives \dfrac{first}{} with the remainder beside it.
 *
 * Does nothing unless the shape is unmistakable: a fraction with a numerator
 * and NO denominator, a separator, and a terminator.
 * ============================================================ */
class DisplayFractionPass : public LinePass {
public:
    void run(NodeList& children, SkipSet& skip, int depth, int prodVer) override;
};

/* ============================================================
 * Pass 1: Fence/Decoration empty-slot merging
 *
 * EQNEDT32 pattern: TMPL(fence/deco, slot[0]=empty) followed by
 * sibling LINE (content) and display chars. This pass:
 *   - Fills FenceNode.content from the next sibling LINE
 *   - Fills DecorationNode.content from the next sibling LINE
 *   - Skips absorbed siblings and display chars
 * ============================================================ */
class FenceMergePass : public LinePass {
public:
    void run(NodeList& children, SkipSet& skip, int depth, int prodVer) override;
    static bool isFenceDisplayChar(const Node* n, int prodVer);
    static bool isBigOpDisplayChar(const Node* n, int prodVer);
};

/* ============================================================
 * Pass 2: BigOp remote display data detection (SIZE_SUB path)
 *
 * Pattern: SIZE_SUB → LINE(limit)* → SIZE_SYM → CHAR(display)*
 * Reverse-scans to find the owning BigOp/Integral template,
 * sets displayLower/displayUpper, skips the block.
 * ============================================================ */
class BigOpDisplayPass : public LinePass {
public:
    void run(NodeList& children, SkipSet& skip, int depth, int prodVer) override;
};

/* ============================================================
 * Pass 2b: BigOp display data without SIZE_SUB
 *
 * EQNEDT32 native pattern: TMPL(BigOp, slot0=empty) → LINE(content)
 * → [LINE(limit)]* → NULL_LINE → SIZE_SYM → CHAR(display)
 * ============================================================ */
class BigOpDisplayAltPass : public LinePass {
public:
    void run(NodeList& children, SkipSet& skip, int depth, int prodVer) override;
};

/* ============================================================
 * Pipeline: runs all passes in order
 * ============================================================ */
class PassPipeline {
public:
    PassPipeline();

    /* Run all passes on a LINE's children */
    void process(NodeList& children, int depth, int prodVer);

private:
    std::vector<std::unique_ptr<LinePass>> passes_;
};

} /* namespace mtef */

#endif /* LINE_PASS_H */
