/*
 * line_pass.cpp -- Pass pipeline implementation
 *
 * Ports the 7-pass convert_line system from mtef2tex.c into
 * independent, testable pass classes.
 */
#include "line_pass.h"
#include <cstring>

namespace mtef {

/* ============================================================
 * Helpers
 * ============================================================ */
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

        /* Skip fence display chars after content */
        if (isFence) {
            for (int j = contentIdx + 1; j < n; j++) {
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

        if (!targetTmpl) continue;

        /* Merge content LINE into template body */
        if (contentIdx >= 0 && !skip[contentIdx]) {
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

/* ============================================================
 * Pipeline construction
 * ============================================================ */
PassPipeline::PassPipeline() {
    /* Order matters: Pass 1 → Pass 2 → Pass 2b */
    passes_.push_back(std::make_unique<FenceMergePass>());
    passes_.push_back(std::make_unique<BigOpDisplayPass>());
    passes_.push_back(std::make_unique<BigOpDisplayAltPass>());
}

void PassPipeline::process(NodeList& children, int depth, int prodVer) {
    if (children.empty()) return;

    SkipSet skip(children.size(), false);

    for (auto& pass : passes_)
        pass->run(children, skip, depth, prodVer);

    /* Remove skipped children */
    NodeList filtered;
    for (size_t i = 0; i < children.size(); i++) {
        if (!skip[i])
            filtered.push_back(std::move(children[i]));
    }
    children = std::move(filtered);
}

} /* namespace mtef */
