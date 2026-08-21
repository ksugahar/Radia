/*
 * line_pass.cpp -- Pass pipeline implementation
 *
 * Ports the 7-pass convert_line system from mtef2tex.c into
 * independent, testable pass classes.
 */
#include "line_pass.h"
#include <iterator>
#include <vector>
#include <cstring>

namespace mtef {

/* ============================================================
 * Helpers
 * ============================================================ */
static bool hasOpenFence(const NodeList& list, int prodVer);

static bool isSlotEmpty(const NodeList& content) {
    if (content.empty()) return true;
    /* Check if all children are null LINEs */
    for (auto& c : content) {
        if (c->tag() != Node::kLine) return false;
        auto* ln = static_cast<const LineNode*>(c.get());
        if (!ln->isNull && !ln->children.empty()) return false;
    }
    return true;
}

/* BigOp display char codes */
static const uint16_t BIGOP_DISPLAY_CHARS[] = {
    0x222B, /* ∫ */  0x222C, /* ∬ */  0x222D, /* ∭ */
    0x222E, /* ∮ */  0x222F, /* ∯ */  0x2230, /* ∰ */
    0x2211, /* ∑ */  0x220F, /* ∏ */  0x2210, /* ∐ */
    0x22C2, /* ⋂ */  0x22C3, /* ⋃ */
    0xFE37, /* ︷ */  0xFE38, /* ︸ */
    0x00E5, /* Σ symbol font */
};
#define BIGOP_DISPLAY_N (sizeof(BIGOP_DISPLAY_CHARS)/sizeof(BIGOP_DISPLAY_CHARS[0]))

static bool inSet(uint16_t code, const uint16_t* set, int n) {
    for (int i = 0; i < n; i++)
        if (set[i] == code) return true;
    return false;
}

bool FenceMergePass::isFenceDisplayChar(const Node* n, int /*prodVer*/) {
    if (n->tag() != Node::kChar) return false;
    auto* ch = static_cast<const CharNode*>(n);
    return (ch->typeface == (TF_DISPLAY & 0x7F));  /* typeface 22 = display */
}

bool FenceMergePass::isBigOpDisplayChar(const Node* n, int /*prodVer*/) {
    if (n->tag() != Node::kChar) return false;
    auto* ch = static_cast<const CharNode*>(n);
    return inSet(ch->charCode, BIGOP_DISPLAY_CHARS, BIGOP_DISPLAY_N);
}

/* ============================================================
 * Pass 1: Fence/Decoration merging
 *
 * For each fence/decoration template with empty slot[0]:
 *   - Find the next non-skipped content LINE
 *   - Move it into the template's content
 *   - Skip fence display chars after the content
 * ============================================================ */
void FenceMergePass::run(NodeList& children, SkipSet& skip, int depth, int prodVer) {
    int n = (int)children.size();

    for (int idx = 0; idx < n; idx++) {
        if (skip[idx]) continue;
        Node* node = children[idx].get();

        /* Only process fence and decoration templates */
        bool isFence = false, isDeco = false;
        int sel = -1;
        NodeList* contentPtr = nullptr;

        if (node->tag() == Node::kFence) {
            auto* f = static_cast<FenceNode*>(node);
            if (isSlotEmpty(f->content)) {
                isFence = true;
                sel = f->selector;
                contentPtr = &f->content;
            }
        } else if (node->tag() == Node::kDecoration) {
            auto* d = static_cast<DecorationNode*>(node);
            if (isSlotEmpty(d->content)) {
                isDeco = true;
                sel = d->selector;
                contentPtr = &d->content;
            }
        } else if (node->tag() == Node::kSqrt) {
            /* Sqrt with empty content (EQNEDT32 pattern: slot[0]=empty,
             * content LINE + SIZE_SUB + LINE(index) as siblings) */
            auto* sq = static_cast<SqrtNode*>(node);
            if (isSlotEmpty(sq->content)) {
                /* Find content LINE */
                for (int j = idx + 1; j < n; j++) {
                    if (skip[j]) continue;
                    if (children[j]->tag() == Node::kLine) {
                        auto* ln = static_cast<LineNode*>(children[j].get());
                        if (!ln->isNull && !ln->children.empty()) {
                            for (auto& c : ln->children)
                                sq->content.push_back(std::move(c));
                            skip[j] = true;
                            break;
                        }
                    }
                    if (children[j]->tag() == Node::kSize) continue;
                    break;
                }
                /* Find index (after SIZE_SUB) */
                for (int j = idx + 1; j < n; j++) {
                    if (skip[j]) continue;
                    if (children[j]->tag() == Node::kSize) {
                        auto* sz = static_cast<SizeNode*>(children[j].get());
                        if (sz->sizeType == SIZETYPE_SUB) {
                            skip[j] = true;
                            /* Next LINE is the index */
                            for (int k = j + 1; k < n; k++) {
                                if (skip[k]) continue;
                                if (children[k]->tag() == Node::kLine) {
                                    auto* ln = static_cast<LineNode*>(children[k].get());
                                    if (!ln->isNull && !ln->children.empty()) {
                                        for (auto& c : ln->children)
                                            sq->index.push_back(std::move(c));
                                        sq->hasIndex = true;
                                        skip[k] = true;
                                    }
                                    break;
                                }
                                if (children[k]->tag() == Node::kSize) { skip[k] = true; continue; }
                                break;
                            }
                        }
                        break;
                    }
                    break;
                }
            }
            continue;
        } else if (node->tag() == Node::kBigOp) {
            /* BigOp with empty body — skip display data chars after it */
            auto* b = static_cast<BigOpNode*>(node);
            if (isSlotEmpty(b->body)) {
                /* Don't merge content, but skip display chars */
                for (int j = idx + 1; j < n; j++) {
                    if (skip[j]) continue;
                    if (isBigOpDisplayChar(children[j].get(), prodVer)) {
                        skip[j] = true;
                    } else if (children[j]->tag() == Node::kSize) {
                        skip[j] = true;
                    } else {
                        break;
                    }
                }
            }
            continue;
        }

        if (!contentPtr) continue;

        /* Find next content LINE */
        int contentIdx = -1;
        for (int j = idx + 1; j < n; j++) {
            if (skip[j]) continue;
            if (children[j]->tag() == Node::kLine) {
                auto* ln = static_cast<LineNode*>(children[j].get());
                if (!ln->isNull) {
                    contentIdx = j;
                    break;
                }
            }
            /* SIZE nodes can appear between template and content */
            if (children[j]->tag() == Node::kSize) continue;
            break;
        }

        if (contentIdx < 0) continue;

        /* Move content LINE into the template */
        auto contentLine = std::move(children[contentIdx]);
        auto* ln = static_cast<LineNode*>(contentLine.get());

        /* For fences: move the LINE's children into fence content */
        contentPtr->clear();
        for (auto& child : ln->children)
            contentPtr->push_back(std::move(child));
        skip[contentIdx] = true;

        /* A fence closes at its own display characters.  Everything between
         * the content LINE and them belongs inside it -- Equation Editor
         * chunks the content, so |r - r_1| arrives as LINE(r), then -, then
         * r_1, then the bars, and one LINE closed the fence round the r.
         *
         * If no display character follows, nothing is absorbed and the
         * one-LINE reading stands: that is what the corpus relies on, and a
         * fence swallowing the rest of the line would be far worse than one
         * that stops early. */
        if (isFence) {
            int close = -1;
            for (int j = contentIdx + 1; j < n; j++) {
                if (skip[j]) continue;
                if (isFenceDisplayChar(children[j].get(), prodVer)) { close = j; break; }
                if (children[j]->tag() == Node::kFence) break;
                if (children[j]->tag() == Node::kFrac) break;
            }
            for (int j = contentIdx + 1; close > 0 && j < close; j++) {
                if (skip[j]) continue;
                if (children[j]->tag() != Node::kSize)
                    contentPtr->push_back(std::move(children[j]));
                skip[j] = true;
            }
            for (int j = (close > 0 ? close : contentIdx + 1); j < n; j++) {
                if (skip[j]) continue;
                if (isFenceDisplayChar(children[j].get(), prodVer)) {
                    skip[j] = true;
                } else if (isBigOpDisplayChar(children[j].get(), prodVer)) {
                    skip[j] = true;
                } else if (children[j]->tag() == Node::kSize) {
                    skip[j] = true;
                } else {
                    break;
                }
            }
        }

        /* Decorations: merge only ONE sibling LINE */
        /* (already done above — contentIdx was the first LINE) */
    }
}

/* ============================================================
 * Pass 2: BigOp remote display data detection (SIZE_SUB path)
 *
 * Scans for SIZE_SUB → LINE(limit)* → SIZE_SYM → CHAR(display)*
 * Reverse-scans to find the owning template and sets
 * displayLower/displayUpper.
 * ============================================================ */
/* The deepest operator in here that has no limits yet.
 *
 * Equation Editor writes the display blocks innermost-last, so when a block
 * turns up with only a LINE before it, the operator it belongs to is the
 * deepest one inside that line that is still bare.
 *
 * "Bare" means no display block has been given to it yet, and none is coming
 * later in its own list.  It deliberately does NOT consult hasLower/hasUpper:
 * those record which SLOTS the template wrote, and an integral written with
 * variation 2 -- by far the common case -- puts its INTEGRAND in the slot this
 * reader calls `upper`.  Testing them made every such integral look like it
 * already had limits, so an operator sitting inside a LINE could never be
 * given the limit block that followed the line.  The sibling path a few lines
 * below takes the first integral it meets with no such test; the two are now
 * consistent, which is what fixed the Harrington plate integrals -- Equation
 * Editor draws them correctly, so the documents were never the problem. */
/* Both shapes of display block end with a switch to symbol size, so an
 * operator with one after it in its own list has its limits coming and must
 * not be given somebody else's. */
static bool symFollows(const NodeList& list, size_t from) {
    for (size_t j = from + 1; j < list.size(); ++j) {
        const Node* n = list[j].get();
        if (n && n->tag() == Node::kSize &&
            static_cast<const SizeNode*>(n)->sizeType == SIZETYPE_SYM)
            return true;
    }
    return false;
}

static Node* deepestBareBigOp(NodeList& list) {
    Node* found = nullptr;
    for (size_t i = 0; i < list.size(); ++i) {
        Node* n = list[i].get();
        if (!n) continue;
        Node* deeper = nullptr;
        if (n->tag() == Node::kLine) {
            auto* ln = static_cast<LineNode*>(n);
            if (!ln->isNull) deeper = deepestBareBigOp(ln->children);
        } else if (n->tag() == Node::kBigOp) {
            auto* b = static_cast<BigOpNode*>(n);
            deeper = deepestBareBigOp(b->body);
            if (!deeper && !b->displayLower && !b->displayUpper &&
                !symFollows(list, i))
                deeper = n;
        } else if (n->tag() == Node::kIntegral) {
            auto* ig = static_cast<IntegralNode*>(n);
            deeper = deepestBareBigOp(ig->body);
            if (!deeper && !ig->displayLower && !ig->displayUpper &&
                !symFollows(list, i))
                deeper = n;
        }
        if (deeper) found = deeper;      /* the last one is the innermost */
    }
    return found;
}

void BigOpDisplayPass::run(NodeList& children, SkipSet& skip, int depth, int prodVer) {
    int n = (int)children.size();

    for (int idx = 0; idx < n; idx++) {
        if (skip[idx]) continue;

        if (children[idx]->tag() != Node::kSize) continue;
        auto* sz = static_cast<SizeNode*>(children[idx].get());
        if (sz->sizeType != SIZETYPE_SUB) continue;

        /* Collect display data block */
        std::vector<int> block;
        std::vector<LineNode*> limitLines;
        bool symFound = false;

        block.push_back(idx);

        for (int j = idx + 1; j < n && block.size() < 128; j++) {
            if (skip[j]) continue;
            Node* c = children[j].get();

            if (c->tag() == Node::kLine) {
                block.push_back(j);
                auto* ln = static_cast<LineNode*>(c);
                if (!ln->isNull && !ln->children.empty())
                    limitLines.push_back(ln);
                continue;
            }
            if (c->tag() == Node::kSize) {
                block.push_back(j);
                auto* s2 = static_cast<SizeNode*>(c);
                if (s2->sizeType == SIZETYPE_SYM) symFound = true;
                else if (s2->sizeType == SIZETYPE_FULL && symFound) break;
                continue;
            }
            if (c->tag() == Node::kChar && symFound) {
                if (FenceMergePass::isBigOpDisplayChar(c, prodVer)) {
                    block.push_back(j);
                    continue;
                }
            }
            break;
        }

        if (!symFound) continue;

        /* Find nearest preceding display template (reverse scan) */
        Node* targetTmpl = nullptr;
        int contentIdx = -1;
        for (int j = idx - 1; j >= 0; j--) {
            if (skip[j]) continue;
            Node* c = children[j].get();

            if (c->tag() == Node::kIntegral || c->tag() == Node::kBigOp) {
                targetTmpl = c;
                break;
            }
            if (contentIdx < 0 && c->tag() == Node::kLine) {
                auto* ln = static_cast<LineNode*>(c);
                if (!ln->isNull) contentIdx = j;
            }
        }

        /* Nothing but a LINE before the block: the operator is inside it. */
        bool targetInsideLine = false;
        if (!targetTmpl && contentIdx >= 0 && !skip[contentIdx]) {
            auto* ln = static_cast<LineNode*>(children[contentIdx].get());
            if (!ln->isNull) {
                targetTmpl = deepestBareBigOp(ln->children);
                targetInsideLine = (targetTmpl != nullptr);
            }
        }

        /* Still nothing, and the nearest line held no operator.  Keep looking
         * back for one that does.
         *
         * The scan above takes the NEAREST line and hopes.  That is enough when
         * the block follows its own content, and wrong when something else sits
         * between: a document with TWO integrals in one line and their two
         * blocks at the end has the second block reverse-scan straight into the
         * FIRST block's limit lines, take `a` for the content, find no operator
         * in it, and give up -- leaving both blocks unclaimed and both
         * integrals bare.  Asking for a line that actually contains a bare
         * operator makes the search purposeful rather than positional.
         *
         * `symFollows` still protects an operator whose own block is coming
         * later in its list, so this cannot steal one that is already spoken
         * for.  The content-merge below is deliberately not attempted for these:
         * a line reached across other material is not this operator's body. */
        if (!targetTmpl) {
            for (int j = idx - 1; j >= 0; j--) {
                if (skip[j]) continue;
                Node* c = children[j].get();
                if (c->tag() != Node::kLine) continue;
                auto* ln = static_cast<LineNode*>(c);
                if (ln->isNull || ln->children.empty()) continue;
                if (Node* found = deepestBareBigOp(ln->children)) {
                    targetTmpl = found;
                    targetInsideLine = true;
                    break;
                }
            }
        }

        if (!targetTmpl) {
            /* No owner.  If the block carries no limits either -- Equation
             * Editor repeats it at the enclosing level with both lines empty
             * -- it can draw nothing, so it goes.  A block with real limits
             * and no owner stays: that is a fault, and it should show. */
            if (limitLines.empty())
                for (int bi : block) skip[bi] = true;
            continue;
        }

        /* Merge content LINE into template body -- but not when the target was
         * found INSIDE that line: the line is the enclosing operator's body,
         * and folding it into the operator it contains would swallow it. */
        if (!targetInsideLine && contentIdx >= 0 && !skip[contentIdx]) {
            auto contentLine = std::move(children[contentIdx]);
            auto* ln = static_cast<LineNode*>(contentLine.get());

            if (targetTmpl->tag() == Node::kIntegral) {
                auto* integ = static_cast<IntegralNode*>(targetTmpl);
                if (isSlotEmpty(integ->body)) {
                    integ->body.clear();
                    for (auto& child : ln->children)
                        integ->body.push_back(std::move(child));
                }
            } else if (targetTmpl->tag() == Node::kBigOp) {
                auto* bigop = static_cast<BigOpNode*>(targetTmpl);
                if (isSlotEmpty(bigop->body)) {
                    bigop->body.clear();
                    for (auto& child : ln->children)
                        bigop->body.push_back(std::move(child));
                }
            }
            skip[contentIdx] = true;
        }

        /* Set limit lines on the target template */
        if (!limitLines.empty()) {
            if (targetTmpl->tag() == Node::kIntegral) {
                auto* integ = static_cast<IntegralNode*>(targetTmpl);
                if (limitLines.size() >= 1) {
                    integ->displayLower = std::make_unique<LineNode>();
                    for (auto& c : limitLines[0]->children)
                        integ->displayLower->children.push_back(std::move(c));
                }
                if (limitLines.size() >= 2) {
                    integ->displayUpper = std::make_unique<LineNode>();
                    for (auto& c : limitLines[1]->children)
                        integ->displayUpper->children.push_back(std::move(c));
                }
            } else if (targetTmpl->tag() == Node::kBigOp) {
                auto* bigop = static_cast<BigOpNode*>(targetTmpl);
                if (limitLines.size() >= 1) {
                    bigop->displayLower = std::make_unique<LineNode>();
                    for (auto& c : limitLines[0]->children)
                        bigop->displayLower->children.push_back(std::move(c));
                }
                if (limitLines.size() >= 2) {
                    bigop->displayUpper = std::make_unique<LineNode>();
                    for (auto& c : limitLines[1]->children)
                        bigop->displayUpper->children.push_back(std::move(c));
                }
            }
        }

        /* Skip all display data */
        for (int bi : block) skip[bi] = true;
    }
}

/* ============================================================
 * Pass 2b: BigOp display data without SIZE_SUB
 * ============================================================ */
void BigOpDisplayAltPass::run(NodeList& children, SkipSet& skip, int depth, int prodVer) {
    int n = (int)children.size();

    for (int idx = 0; idx < n; idx++) {
        if (skip[idx]) continue;

        Node* node = children[idx].get();
        bool isDisplayTmpl = false;

        if (node->tag() == Node::kIntegral || node->tag() == Node::kBigOp)
            isDisplayTmpl = true;
        if (!isDisplayTmpl) continue;

        /* Check if slot[0]/body is empty */
        bool bodyEmpty = false;
        if (node->tag() == Node::kIntegral) {
            bodyEmpty = isSlotEmpty(static_cast<IntegralNode*>(node)->body);
        } else {
            bodyEmpty = isSlotEmpty(static_cast<BigOpNode*>(node)->body);
        }
        if (!bodyEmpty) continue;

        /* Find next content LINE (may be empty) */
        int contentIdx = -1;
        for (int j = idx + 1; j < n; j++) {
            if (skip[j]) continue;
            if (children[j]->tag() == Node::kLine) {
                auto* ln = static_cast<LineNode*>(children[j].get());
                if (!ln->isNull) { contentIdx = j; break; }
            }
            if (children[j]->tag() == Node::kSize) continue;
            break;
        }
        if (contentIdx < 0) continue;

        /* Collect limit LINEs and display data */
        std::vector<LineNode*> limitLines;
        std::vector<int> block;
        bool dispDataFound = false;

        for (int j = contentIdx + 1; j < n && block.size() < 128; j++) {
            if (skip[j]) continue;
            Node* c = children[j].get();

            if (c->tag() == Node::kLine && static_cast<LineNode*>(c)->isNull) {
                dispDataFound = true;
                block.push_back(j);
                continue;
            }
            if (c->tag() == Node::kLine && !static_cast<LineNode*>(c)->isNull) {
                if (!dispDataFound && limitLines.size() < 2) {
                    limitLines.push_back(static_cast<LineNode*>(c));
                    block.push_back(j);
                    continue;
                }
                break;
            }
            if (c->tag() == Node::kSize) {
                block.push_back(j);
                auto* s2 = static_cast<SizeNode*>(c);
                if (s2->sizeType == SIZETYPE_SYM) dispDataFound = true;
                continue;
            }
            if (c->tag() == Node::kChar && dispDataFound &&
                FenceMergePass::isBigOpDisplayChar(c, prodVer)) {
                block.push_back(j);
                continue;
            }
            break;
        }

        if (!dispDataFound) continue;

        /* Merge content LINE into body */
        auto contentLine = std::move(children[contentIdx]);
        auto* ln = static_cast<LineNode*>(contentLine.get());

        if (node->tag() == Node::kIntegral) {
            auto* integ = static_cast<IntegralNode*>(node);
            integ->body.clear();
            for (auto& child : ln->children)
                integ->body.push_back(std::move(child));
        } else {
            auto* bigop = static_cast<BigOpNode*>(node);
            bigop->body.clear();
            for (auto& child : ln->children)
                bigop->body.push_back(std::move(child));
        }
        skip[contentIdx] = true;

        /* Set limits */
        if (!limitLines.empty()) {
            if (node->tag() == Node::kIntegral) {
                auto* integ = static_cast<IntegralNode*>(node);
                if (limitLines.size() >= 1) {
                    integ->displayLower = std::make_unique<LineNode>();
                    for (auto& c : limitLines[0]->children)
                        integ->displayLower->children.push_back(std::move(c));
                }
                if (limitLines.size() >= 2) {
                    integ->displayUpper = std::make_unique<LineNode>();
                    for (auto& c : limitLines[1]->children)
                        integ->displayUpper->children.push_back(std::move(c));
                }
            } else {
                auto* bigop = static_cast<BigOpNode*>(node);
                if (limitLines.size() >= 1) {
                    bigop->displayLower = std::make_unique<LineNode>();
                    for (auto& c : limitLines[0]->children)
                        bigop->displayLower->children.push_back(std::move(c));
                }
                if (limitLines.size() >= 2) {
                    bigop->displayUpper = std::make_unique<LineNode>();
                    for (auto& c : limitLines[1]->children)
                        bigop->displayUpper->children.push_back(std::move(c));
                }
            }
        }

        /* Skip all display data */
        for (int bi : block) skip[bi] = true;
    }
}



/* Does this list hold a bracket that has not closed?  A fence arrives with an
 * empty slot and is filled from its siblings later, so an empty one that has
 * not met its display characters yet is a bracket still open. */

static bool hasOpenFenceNode(const Node* n, int prodVer) {
    if (!n) return false;
    if (n->tag() == Node::kFence)
        return static_cast<const FenceNode*>(n)->content.empty();
    if (n->tag() == Node::kLine)
        return hasOpenFence(static_cast<const LineNode*>(n)->children, prodVer);
    return false;
}

static bool hasOpenFence(const NodeList& list, int prodVer) {
    bool open = false;
    for (const auto& n : list) {
        if (hasOpenFenceNode(n.get(), prodVer)) open = true;
        if (n && FenceMergePass::isFenceDisplayChar(n.get(), prodVer)) open = false;
        if (n && n->tag() == Node::kLine &&
            !hasOpenFence(static_cast<const LineNode*>(n.get())->children, prodVer))
            continue;
    }
    return open;
}



/* ============================================================
 * Pass 0b: Matrix cell redistribution
 * ============================================================ */
void MatrixCellPass::run(NodeList& children, SkipSet& skip,
                         int /*depth*/, int /*prodVer*/) {
    for (size_t idx = 0; idx < children.size(); ++idx) {
        if (skip[idx] || !children[idx]) continue;
        if (children[idx]->tag() != Node::kMatrix) continue;
        auto* m = static_cast<MatrixNode*>(children[idx].get());

        const size_t cells = m->elements.size();
        if (cells < 2) continue;

        /* Exactly one cell with anything in it. */
        size_t full = cells;
        size_t nonEmpty = 0;
        for (size_t i = 0; i < cells; ++i) {
            const Node* e = m->elements[i].get();
            if (!e || e->tag() != Node::kLine) return;   /* not the shape */
            if (!static_cast<const LineNode*>(e)->children.empty()) {
                ++nonEmpty;
                full = i;
            }
        }
        if (nonEmpty != 1) continue;

        /* Holding a line per cell, and no more than there are cells. */
        auto* src = static_cast<LineNode*>(m->elements[full].get());
        size_t lines = 0;
        for (const auto& c : src->children) {
            if (!c || c->tag() != Node::kLine) { lines = 0; break; }
            ++lines;
        }
        if (lines < 2 || lines > cells) continue;

        NodeList rows = std::move(src->children);
        src->children.clear();
        size_t at = 0;
        for (auto& row : rows) {
            auto* dst = static_cast<LineNode*>(m->elements[at].get());
            auto* rl = static_cast<LineNode*>(row.get());
            dst->children = std::move(rl->children);
            ++at;
        }
    }
}

/* ============================================================
 * Pass 0: Integral/BigOp slot splitting
 * ============================================================ */

/* The line's children, or the node itself if it is not a line. */
static void appendFlat(NodeList& into, NodePtr node) {
    if (node && node->tag() == Node::kLine) {
        auto* ln = static_cast<LineNode*>(node.get());
        if (!ln->isNull) {
            for (auto& c : ln->children) into.push_back(std::move(c));
            return;
        }
    }
    if (node) into.push_back(std::move(node));
}

/* Split one slot that holds "integrand, SUB, lower, upper, SYM, signs".
 * Returns false and touches nothing when the slot is not that shape. */
static bool splitDisplaySlot(NodeList& slot, NodeList& body,
                             NodeList& lower, bool& hasLower,
                             NodeList& upper, bool& hasUpper) {
    int subIdx = -1;
    for (size_t i = 0; i < slot.size(); ++i) {
        Node* n = slot[i].get();
        if (n && n->tag() == Node::kSize &&
            static_cast<SizeNode*>(n)->sizeType == SIZETYPE_SUB) {
            subIdx = int(i);
            break;
        }
    }
    if (subIdx <= 0) return false;          /* no marker, or nothing before it */

    /* The slot is one of the fields about to be written -- the variation sent
     * the content to `upper`, so `upper` and `slot` are the same list.  Take
     * it out first; writing the limits and THEN clearing the slot cleared the
     * limit with it. */
    NodeList src = std::move(slot);
    slot.clear();

    NodeList newBody;
    for (int i = 0; i < subIdx; ++i) {
        if (src[i] && src[i]->tag() == Node::kSize) continue;
        appendFlat(newBody, std::move(src[i]));
    }
    if (newBody.empty()) {                  /* nothing to be a body: undo */
        slot = std::move(src);
        return false;
    }

    /* After the marker: the limit lines, lower first, until symbol size. */
    NodeList lines[2];
    int seen = 0;
    for (size_t i = size_t(subIdx) + 1; i < src.size(); ++i) {
        Node* n = src[i].get();
        if (!n) continue;
        if (n->tag() == Node::kSize) {
            if (static_cast<SizeNode*>(n)->sizeType == SIZETYPE_SYM) break;
            continue;
        }
        if (n->tag() != Node::kLine) break;
        auto* ln = static_cast<LineNode*>(n);
        if (seen < 2) {
            if (!ln->isNull)
                for (auto& c : ln->children) lines[seen].push_back(std::move(c));
            ++seen;
        }
    }

    /* The body goes back as a LINE.  A slot is a plain list and the pass
     * pipeline runs on a line's children, so left loose nothing would run on
     * it -- and a second integral nested in the first needs the same repair
     * this one just had. */
    auto wrap = std::make_unique<LineNode>();
    wrap->children = std::move(newBody);
    body.clear();
    body.push_back(std::move(wrap));

    lower = std::move(lines[0]);
    hasLower = !lower.empty();
    upper = std::move(lines[1]);
    hasUpper = !upper.empty();
    return true;
}

/* SUB, then limit lines, then SYM: a display block waiting in this list. */
static bool displayBlockFollows(const NodeList& list, size_t from) {
    for (size_t j = from + 1; j < list.size(); ++j) {
        const Node* n = list[j].get();
        if (!n) continue;
        if (n->tag() == Node::kLine) continue;
        if (n->tag() != Node::kSize) return false;
        const int t = static_cast<const SizeNode*>(n)->sizeType;
        if (t == SIZETYPE_SYM) return true;
        if (t != SIZETYPE_SUB) return false;
    }
    return false;
}

void IntegralSlotPass::run(NodeList& children, SkipSet& skip,
                           int /*depth*/, int /*prodVer*/) {
    for (size_t idx = 0; idx < children.size(); ++idx) {
        if (skip[idx] || !children[idx]) continue;
        Node* n = children[idx].get();

        if (n->tag() == Node::kIntegral) {
            auto* ig = static_cast<IntegralNode*>(n);
            /* Whichever slot the variation sent it to. */
            NodeList* cand[3] = {&ig->body, &ig->upper, &ig->lower};
            bool split = false;
            for (NodeList* sl : cand) {
                if (sl->empty()) continue;
                if (splitDisplaySlot(*sl, ig->body, ig->lower, ig->hasLower,
                                     ig->upper, ig->hasUpper)) {
                    split = true;
                    break;
                }
            }
            /* Not that shape -- but if the body is empty, exactly one limit
             * holds something, and a display block is coming in this list to
             * supply the real limits, then what is in that slot is the
             * integrand and belongs in the body. */
            if (!split && ig->body.empty() &&
                ig->upper.empty() != ig->lower.empty() &&
                displayBlockFollows(children, idx)) {
                NodeList& used = ig->upper.empty() ? ig->lower : ig->upper;
                ig->body = std::move(used);
                used.clear();
                ig->hasLower = ig->hasUpper = false;
            }
        } else if (n->tag() == Node::kBigOp) {
            auto* b = static_cast<BigOpNode*>(n);
            NodeList* cand[3] = {&b->body, &b->upper, &b->lower};
            for (NodeList* sl : cand) {
                if (sl->empty()) continue;
                if (splitDisplaySlot(*sl, b->body, b->lower, b->hasLower,
                                     b->upper, b->hasUpper)) {
                    b->hasLimits = true;
                    break;
                }
            }
        }
    }
}

/* ============================================================
 * Pass 0a: Display-fraction reassembly
 * ============================================================ */
void DisplayFractionPass::run(NodeList& children, SkipSet& skip,
                              int /*depth*/, int prodVer) {
    const int n = int(children.size());

    auto nextLive = [&](int from) {
        for (int j = from; j < n; ++j)
            if (!skip[j]) return j;
        return -1;
    };

    for (int idx = 0; idx < n; ++idx) {
        if (skip[idx]) continue;
        if (children[idx]->tag() != Node::kFrac) continue;
        auto* f = static_cast<FracNode*>(children[idx].get());
        if (f->numer.empty()) continue;

        if (!f->denom.empty()) {
            /* An inline fraction usually carries both parts and needs
             * nothing -- unless a bracket inside its denominator has not
             * closed, in which case the denominator carries on past the
             * template and the closing characters are the end of it. */
            if (!hasOpenFence(f->denom, prodVer)) continue;
            int close = -1;
            for (int j = idx + 1; j < n; ++j) {
                if (skip[j]) continue;
                if (FenceMergePass::isFenceDisplayChar(children[j].get(), prodVer)) {
                    close = j;
                    break;
                }
                if (children[j]->tag() == Node::kFrac) break;
                if (children[j]->tag() == Node::kFence) break;
            }
            if (close < 0) continue;      /* never closes: leave it alone */
            /* An exponent on the bracket follows its closing characters. */
            int last = close;
            for (int j = close + 1; j < n; ++j) {
                if (skip[j]) continue;
                if (children[j]->tag() == Node::kScript) last = j;
                break;
            }
            NodeList flat;
            for (auto& c : f->denom) flat.push_back(std::move(c));
            for (int j = idx + 1; j <= last; ++j) {
                if (skip[j]) continue;
                skip[j] = true;
                if (children[j]->tag() == Node::kSize) continue;
                flat.push_back(std::move(children[j]));
            }
            auto wrap = std::make_unique<LineNode>();
            wrap->children = std::move(flat);
            f->denom.clear();
            f->denom.push_back(std::move(wrap));
            continue;
        }

        /* The separator between the two parts. */
        int sep = -1;
        for (int j = idx + 1; j < n; ++j) {
            if (skip[j]) continue;
            if (children[j]->tag() == Node::kSize) { sep = j; break; }
            if (children[j]->tag() == Node::kFrac) {
                auto* g = static_cast<FracNode*>(children[j].get());
                if (g->denom.empty()) { sep = -2; break; }  /* another one */
            }
        }
        /* No separator at all, and the very next thing is one LINE: that
         * line IS the denominator.
         *
         * This is the third shape Equation Editor writes a display fraction
         * in -- numerator inside the template, denominator immediately after
         * it, no size markers between.  The two shapes below need a separator
         * and a terminator to know where the denominator ENDS, and refuse
         * without them, which is right when the boundary is a guess.  Here it
         * is not a guess: a LINE is one self-contained chunk, and the
         * denominator is exactly that chunk.
         *
         * Narrow on purpose: the denominator must be empty (so nothing is
         * displaced), and the very next live sibling must be a real LINE.
         * Anything else still refuses. */
        if (sep < 0) {
            const int m = nextLive(idx + 1);
            if (m < 0 || children[m]->tag() != Node::kLine) continue;
            auto* ln = static_cast<LineNode*>(children[m].get());
            if (ln->isNull || ln->children.empty()) continue;
            f->denom.clear();
            f->denom.push_back(std::move(children[m]));
            skip[m] = true;
            f->display = true;
            f->styleOverride = 1;
            continue;
        }

        /* The terminator.  A size record ends the denominator UNLESS the
         * display characters of a fence opened inside it follow -- the fence
         * has not closed, so neither has the denominator. */
        int endIdx = -1;
        for (int k = sep + 1; k < n; ++k) {
            if (skip[k]) continue;
            if (children[k]->tag() != Node::kSize) continue;
            const int m = nextLive(k + 1);
            if (m >= 0 &&
                FenceMergePass::isFenceDisplayChar(children[m].get(), prodVer))
                continue;
            endIdx = k;
            break;
        }
        if (endIdx < 0) continue;

        /* Nothing between the separators means there is no denominator to
         * recover, and moving an empty one would only hide the fault. */
        bool anyDenom = false;
        for (int k = sep + 1; k < endIdx; ++k)
            if (!skip[k] && children[k]->tag() != Node::kSize) anyDenom = true;
        if (!anyDenom) continue;

        /* Commit.  Size records are markers, not content: dropping them here
         * is also what stops a \displaystyle appearing in the middle of the
         * LaTeX for every one of these. */
        /* Each part is gathered into ONE line.
         *
         * A bare LINE among the pieces is another of Equation Editor's
         * chunks; its contents have to become siblings of what follows or the
         * fence in the denominator cannot close round the "- r_1" after it.
         * But a slot is a plain list and the pass pipeline runs on a LINE's
         * children, so splicing them in loose meant no pass ran on them at
         * all and the fence closed round nothing.  One line gives both. */
        auto gather = [&](NodeList& into, int from, int to) {
            /* The pieces are moved AS THEY ARE.  process() has already
             * spliced the chunk LINEs, and it is the one that knows which
             * LINE a template is waiting for; opening them again here took
             * the content away from the bracket that was waiting for it. */
            NodeList flat;
            for (auto& n : into) {
                if (n && n->tag() == Node::kLine) {
                    auto* ln = static_cast<LineNode*>(n.get());
                    if (!ln->isNull && !into.empty() && into.size() == 1) {
                        for (auto& c : ln->children) flat.push_back(std::move(c));
                        continue;
                    }
                }
                flat.push_back(std::move(n));
            }
            for (int j = from; j < to; ++j) {
                if (skip[j]) continue;
                skip[j] = true;
                if (children[j]->tag() == Node::kSize) continue;
                flat.push_back(std::move(children[j]));
            }
            auto wrap = std::make_unique<LineNode>();
            wrap->children = std::move(flat);
            into.clear();
            into.push_back(std::move(wrap));
        };
        gather(f->numer, idx + 1, sep);
        gather(f->denom, sep + 1, endIdx);
        skip[sep] = true;
        skip[endIdx] = true;
        f->display = true;
        f->styleOverride = 1;
    }
}

/* ============================================================
 * Pipeline construction
 * ============================================================ */
PassPipeline::PassPipeline() {
    /* Order matters: Pass 0a → Pass 1 → Pass 2 → Pass 2b.  The display
     * fraction goes first: it puts the fence back inside the denominator,
     * where Pass 1 can then fill it from its own siblings. */
    passes_.push_back(std::make_unique<MatrixCellPass>());
    passes_.push_back(std::make_unique<IntegralSlotPass>());
    passes_.push_back(std::make_unique<DisplayFractionPass>());
    passes_.push_back(std::make_unique<FenceMergePass>());
    passes_.push_back(std::make_unique<BigOpDisplayPass>());
    passes_.push_back(std::make_unique<BigOpDisplayAltPass>());
}


/* Does anything after this LINE continue it?  A script with no base has to
 * attach to an atom inside it, and a bracket's display characters close a
 * bracket opened in it: neither can begin a line, so the LINE is a chunk of
 * something larger.  Anything else after it and the LINE is a line. */
static bool continuedBy(const NodeList& list, size_t i, int prodVer) {
    for (size_t j = i + 1; j < list.size(); ++j) {
        const Node* n = list[j].get();
        if (!n) continue;
        /* Markers say nothing either way.  A block of remote limits was
         * briefly treated as a continuation here, which spliced the line open
         * and took the outer operator apart to reach the inner one; the block
         * is found by reading INTO the line instead (Pass 2). */
        if (n->tag() == Node::kSize) continue;
        if (FenceMergePass::isFenceDisplayChar(n, prodVer)) return true;
        if (n->tag() == Node::kScript)
            return static_cast<const ScriptNode*>(n)->base.empty();
        return false;
    }
    return false;
}

/* Is this a template still waiting for its content to arrive as the next
 * sibling?  Those LINEs are not chunks and must stay where they are. */
static bool awaitsContent(const Node* n) {
    if (!n) return false;
    switch (n->tag()) {
        case Node::kFence:
            return static_cast<const FenceNode*>(n)->content.empty();
        case Node::kDecoration:
            return static_cast<const DecorationNode*>(n)->content.empty();
        case Node::kSqrt:
            return static_cast<const SqrtNode*>(n)->content.empty();
        case Node::kBigOp:
            return static_cast<const BigOpNode*>(n)->body.empty();
        case Node::kIntegral:
            return static_cast<const IntegralNode*>(n)->body.empty();
        case Node::kBraceDeco:
            return static_cast<const BraceDecoNode*>(n)->content.empty();
        case Node::kScript:
            return static_cast<const ScriptNode*>(n)->base.empty();
        case Node::kFrac: {
            const auto* f = static_cast<const FracNode*>(n);
            return f->numer.empty() || f->denom.empty();
        }
        case Node::kSize:
            return false;
        default:
            return false;
    }
}

/* The slots of a node, for walking into.  Not every kind has one; the ones
 * that do are the templates whose contents a pass may need to see. */
static void collectSlots(Node* n, std::vector<NodeList*>& out) {
    if (!n) return;
    switch (n->tag()) {
        case Node::kFence: {
            auto* f = static_cast<FenceNode*>(n);
            out.push_back(&f->content);
            if (f->hasMiddle) out.push_back(&f->content2);
            break;
        }
        case Node::kFrac: {
            auto* f = static_cast<FracNode*>(n);
            out.push_back(&f->numer);
            out.push_back(&f->denom);
            break;
        }
        case Node::kSqrt: {
            auto* q = static_cast<SqrtNode*>(n);
            out.push_back(&q->content);
            if (q->hasIndex) out.push_back(&q->index);
            break;
        }
        case Node::kScript: {
            auto* sc = static_cast<ScriptNode*>(n);
            out.push_back(&sc->base);
            if (sc->hasSub) out.push_back(&sc->sub);
            if (sc->hasSup) out.push_back(&sc->sup);
            break;
        }
        case Node::kIntegral: {
            auto* ig = static_cast<IntegralNode*>(n);
            out.push_back(&ig->body);
            if (ig->hasLower) out.push_back(&ig->lower);
            if (ig->hasUpper) out.push_back(&ig->upper);
            break;
        }
        case Node::kBigOp: {
            auto* b = static_cast<BigOpNode*>(n);
            out.push_back(&b->body);
            if (b->hasLower) out.push_back(&b->lower);
            if (b->hasUpper) out.push_back(&b->upper);
            break;
        }
        case Node::kEmbell:
            out.push_back(&static_cast<EmbellNode*>(n)->content);
            break;
        case Node::kDecoration:
            out.push_back(&static_cast<DecorationNode*>(n)->content);
            break;
        case Node::kBraceDeco: {
            auto* b = static_cast<BraceDecoNode*>(n);
            out.push_back(&b->content);
            out.push_back(&b->label);
            break;
        }
        case Node::kOverset: {
            auto* o = static_cast<OversetNode*>(n);
            out.push_back(&o->over);
            out.push_back(&o->base);
            break;
        }
        case Node::kPhantom:
            out.push_back(&static_cast<PhantomNode*>(n)->content);
            break;
        default:
            break;
    }
}

/* A slot that is one LINE is walked by whoever walks that line; anything else
 * has to be walked here or no pass ever sees it. */
static bool slotNeedsWalking(const NodeList& slot) {
    if (slot.empty()) return false;
    if (slot.size() == 1 && slot[0] && slot[0]->tag() == Node::kLine)
        return false;
    return true;
}

void PassPipeline::process(NodeList& children, int depth, int prodVer) {
    if (children.empty()) return;

    /* A pile of ONE line is Equation Editor chunking its stream, not a stack:
     * it renders as its line either way.  Splicing it open first is what puts
     * a display fraction and its continuation in the same list -- in the
     * lecture file the fraction sits inside the pile and its denominator
     * after it, so no pass working on one list could see both. */
    for (size_t i = 0; i < children.size(); ++i) {
        if (!children[i] || children[i]->tag() != Node::kPile) continue;
        auto* p = static_cast<PileNode*>(children[i].get());
        if (p->lines.size() != 1 || !p->lines[0]) continue;
        if (p->lines[0]->tag() != Node::kLine) continue;
        auto* ln = static_cast<LineNode*>(p->lines[0].get());
        if (ln->isNull) continue;
        NodeList inner = std::move(ln->children);
        children.erase(children.begin() + long(i));
        children.insert(children.begin() + long(i),
                        std::make_move_iterator(inner.begin()),
                        std::make_move_iterator(inner.end()));
        --i;   /* look at what was spliced in, in case it is another one */
    }

    /* Chunk LINEs go the same way -- except the one a template is waiting
     * for.  Pass 1 finds a fence's content by it being the next sibling, so
     * opening that LINE would take the content away from its own bracket. */
    for (size_t i = 0; i < children.size(); ++i) {
        if (!children[i] || children[i]->tag() != Node::kLine) continue;
        auto* ln = static_cast<LineNode*>(children[i].get());
        if (ln->isNull || ln->children.empty()) continue;
        if (i > 0 && awaitsContent(children[i - 1].get())) continue;
        if (!continuedBy(children, i, prodVer) &&
            !hasOpenFence(ln->children, prodVer)) continue;
        NodeList inner = std::move(ln->children);
        children.erase(children.begin() + long(i));
        children.insert(children.begin() + long(i),
                        std::make_move_iterator(inner.begin()),
                        std::make_move_iterator(inner.end()));
        --i;
    }

    SkipSet skip(children.size(), false);

    for (auto& pass : passes_)
        pass->run(children, skip, depth, prodVer);

    /* A drop to script size that opens nothing but an empty line, and never
     * reaches a symbol-size switch, is not a display block -- it is the
     * remains of one, and left alone it sets the rest of the equation small. */
    for (size_t i = 0; i + 1 < children.size(); ++i) {
        if (skip[i] || !children[i]) continue;
        if (children[i]->tag() != Node::kSize) continue;
        if (static_cast<SizeNode*>(children[i].get())->sizeType != SIZETYPE_SUB)
            continue;
        std::vector<size_t> empties;
        bool block = false, ok = true;
        for (size_t j = i + 1; j < children.size(); ++j) {
            if (skip[j] || !children[j]) continue;
            Node* n = children[j].get();
            if (n->tag() == Node::kSize) {
                if (static_cast<SizeNode*>(n)->sizeType == SIZETYPE_SYM)
                    block = true;
                break;
            }
            if (n->tag() != Node::kLine) break;
            auto* ln = static_cast<LineNode*>(n);
            if (ln->isNull || ln->children.empty()) { empties.push_back(j); continue; }
            ok = false;                     /* a real limit: leave it alone */
            break;
        }
        if (block || !ok || empties.empty()) continue;
        skip[i] = true;
        for (size_t j : empties) skip[j] = true;
    }

    /* A size marker at the very end switches a style nothing follows, and an
     * empty line at the end sets nothing: both are the remains of a display
     * block whose operator turned out to be elsewhere.  Drop them, while
     * something else is left to keep. */
    for (int i = int(children.size()) - 1; i > 0; --i) {
        if (skip[size_t(i)]) continue;
        Node* n = children[size_t(i)].get();
        if (!n) continue;
        const bool dead =
            n->tag() == Node::kSize ||
            (n->tag() == Node::kLine &&
             (static_cast<LineNode*>(n)->isNull ||
              static_cast<LineNode*>(n)->children.empty()));
        if (!dead) break;
        skip[size_t(i)] = true;
    }

    /* Remove skipped children */
    NodeList filtered;
    for (size_t i = 0; i < children.size(); i++) {
        if (!skip[i])
            filtered.push_back(std::move(children[i]));
    }
    children = std::move(filtered);

    /* And into the slots, which nothing else walks. */
    if (depth < 32) {
        for (auto& child : children) {
            std::vector<NodeList*> slots;
            collectSlots(child.get(), slots);
            for (NodeList* sl : slots)
                if (slotNeedsWalking(*sl))
                    process(*sl, depth + 1, prodVer);
        }
    }
}

} /* namespace mtef */
