/*
 * mtef_parser.cpp -- MTEF binary parser implementation
 *
 * Faithfully ports the parsing logic from mtef2tex.c (parse_record,
 * parse_char, parse_tmpl, parse_line, parse_pile, etc.) into the
 * unified C++ node hierarchy.
 */
#include "mtef_parser.h"
#include <cstring>

namespace mtef {

/* ============================================================
 * Construction
 * ============================================================ */
MtefParser::MtefParser(const uint8_t* data, size_t len)
    : data_(data), len_(len) {}

/* ============================================================
 * Low-level byte reading
 * ============================================================ */
int MtefParser::readByte() {
    if (pos_ >= len_) return -1;
    return data_[pos_++];
}

int MtefParser::readUint16() {
    if (pos_ + 1 >= len_) return -1;
    int lo = data_[pos_++];
    int hi = data_[pos_++];
    return lo | (hi << 8);
}

void MtefParser::skipNudge() {
    /* Nudge data: 4 bytes (dx_lo, dx_hi, dy_lo, dy_hi) */
    if (pos_ + 4 <= len_) pos_ += 4;
}

/* ============================================================
 * Object list parser (reads until REC_END)
 * ============================================================ */
NodeList MtefParser::parseObjectList() {
    NodeList list;
    while (pos_ < len_) {
        int tag = readByte();
        if (tag < 0) break;
        int recType = tag & 0x0F;
        int options = (tag >> 4) & 0x0F;
        if (recType == REC_END) break;

        auto node = parseRecord(recType, options);
        if (node) list.push_back(std::move(node));
    }
    return list;
}

/* ============================================================
 * Record dispatch
 * ============================================================ */
NodePtr MtefParser::parseRecord(int recType, int options) {
    switch (recType) {
    case REC_LINE:   return parseLine(options);
    case REC_CHAR:   return parseChar(options);
    case REC_TMPL:   return parseTmpl(options);
    case REC_PILE:   return parsePile(options);
    case REC_MATRIX: return parseMatrix(options);
    case REC_EMBELL: return parseEmbell(options);
    case REC_SIZE:   return parseSize(options);
    case REC_FULL:   { auto s = std::make_unique<SizeNode>(); s->sizeType = SIZETYPE_FULL; return s; }
    case REC_SUB:    { auto s = std::make_unique<SizeNode>(); s->sizeType = SIZETYPE_SUB; return s; }
    case REC_SUB2:   { auto s = std::make_unique<SizeNode>(); s->sizeType = SIZETYPE_SUB2; return s; }
    case REC_SYM:    { auto s = std::make_unique<SizeNode>(); s->sizeType = SIZETYPE_SYM; return s; }
    case REC_SUBSYM: { auto s = std::make_unique<SizeNode>(); s->sizeType = SIZETYPE_SUBSYM; return s; }
    case REC_FONT:   return parseFont(options);
    case REC_RULER:
        /* Skip ruler data: n_stops * 3 bytes */
        {
            int n = readByte();
            if (n > 0) {
                for (int i = 0; i < n; i++) {
                    readByte();   /* type */
                    readUint16(); /* offset */
                }
            }
        }
        return nullptr;
    default:
        return nullptr;
    }
}

/* ============================================================
 * LINE parser
 * ============================================================ */
std::unique_ptr<LineNode> MtefParser::parseLine(int options) {
    if (options & OPT_NUDGE) skipNudge();
    if (options & OPT_LINE_LSPACE) readUint16();

    auto n = std::make_unique<LineNode>();

    if (options & OPT_LINE_NULL) {
        n->isNull = true;
        return n;
    }

    n->children = parseObjectList();
    return n;
}

/* ============================================================
 * CHAR parser
 * ============================================================ */
std::unique_ptr<CharNode> MtefParser::parseChar(int options) {
    if (options & OPT_NUDGE) skipNudge();

    int tfByte = readByte();
    if (tfByte < 0) return nullptr;
    int typeface = tfByte - 128;
    int charCode = readUint16();
    if (charCode < 0) return nullptr;

    auto n = std::make_unique<CharNode>();
    n->typeface = typeface;
    n->charCode = static_cast<uint16_t>(charCode);

    if (options & OPT_CHAR_EMBELL)
        n->embells = parseEmbellList();

    return n;
}

/* ============================================================
 * EMBELL parser
 * ============================================================ */
std::unique_ptr<EmbellNode> MtefParser::parseEmbell(int options) {
    if (options & OPT_NUDGE) skipNudge();
    int etype = readByte();
    if (etype < 0) return nullptr;

    auto n = std::make_unique<EmbellNode>();
    n->embellType = etype;
    return n;
}

/* Embellishment chain: reads consecutive EMBELL records */
std::vector<int> MtefParser::parseEmbellList() {
    std::vector<int> embells;
    while (pos_ < len_) {
        int tag = readByte();
        if (tag < 0) break;
        int recType = tag & 0x0F;
        int options = (tag >> 4) & 0x0F;
        if (recType == REC_END) { pos_--; break; }
        if (recType == REC_EMBELL) {
            if (options & OPT_NUDGE) skipNudge();
            int etype = readByte();
            if (etype >= 0) embells.push_back(etype);
        } else {
            pos_--;
            break;
        }
    }
    return embells;
}

/* ============================================================
 * SIZE parser
 * ============================================================ */
std::unique_ptr<SizeNode> MtefParser::parseSize(int options) {
    auto n = std::make_unique<SizeNode>();
    /* typesize is encoded in options (high nibble of tag) */
    n->sizeType = options;

    /* lsize == 101 (0x65): explicit size in half-points follows */
    if (pos_ < len_) {
        int peek = data_[pos_];
        if (peek == 101) {
            pos_++;
            readUint16(); /* skip the explicit size value */
        }
    }
    return n;
}

/* ============================================================
 * FONT parser
 * ============================================================ */
std::unique_ptr<FontNode> MtefParser::parseFont(int options) {
    if (options & OPT_NUDGE) skipNudge();

    auto n = std::make_unique<FontNode>();
    int enc_def = readByte();
    if (enc_def < 0) return n;

    /* Font number (optional) */
    if (enc_def & 0x04) {
        int fontnum = readByte();
        if (fontnum >= 0) n->fontIndex = fontnum;
    }

    /* Font name */
    if (enc_def & 0x01) {
        int nameLen = 0;
        while (pos_ < len_ && nameLen < 63) {
            int ch = data_[pos_++];
            if (ch == 0) break;
            n->name[nameLen++] = static_cast<char>(ch);
        }
        n->name[nameLen] = 0;
    }

    /* Style */
    if (enc_def & 0x02) {
        n->style = readByte();
    }

    return n;
}

/* ============================================================
 * Template slot count (from mtef2tex.c get_template_slot_count)
 * ============================================================ */
int MtefParser::getSlotCount(int selector, int variation) {
    /* Fences */
    if (selector >= tmANGLE && selector <= tmLPRB) return 1;

    /* Root */
    if (selector == tmROOT) return (variation == 0) ? 2 : 1;

    /* Fractions */
    if (selector == tmFRACT || selector == tmSLFRACT) return 2;

    /* Scripts */
    if (selector == tmSCRIPT || selector == tmLSCRIPT) return 2;

    /* Decorations */
    if (selector == tmUBAR || selector == tmOBAR) return 1;
    if (selector == tmLARROW || selector == tmRARROW || selector == tmBARROW) return 1;

    /* Integrals */
    if (selector >= tmSINT && selector <= tmTSINT) {
        int hasLower = variation & 0x01;
        int hasUpper = variation & 0x02;
        return 1 + (hasLower ? 1 : 0) + (hasUpper ? 1 : 0);
    }

    /* Horizontal braces */
    if (selector == tmUHBRACE || selector == tmLHBRACE) return 2;

    /* Sums/Products */
    if (selector >= tmSUM && selector <= tmIINTER) return 1;

    /* Lim */
    if (selector == tmLIM) return (variation == 2) ? 3 : 2;

    /* Long division */
    if (selector == tmLDIV) return (variation == 0) ? 2 : 1;

    /* IntOp/SumOp */
    if (selector == tmINTOP || selector == tmSUMOP) {
        int hasLower = variation & 0x01;
        int hasUpper = variation & 0x02;
        return (hasLower ? 1 : 0) + (hasUpper ? 1 : 0) + 1;
    }

    /* Dirac */
    if (selector == tmDIRAC) return (variation == 0) ? 2 : 1;

    /* Arrows */
    if (selector == tmUARROW || selector == tmOARROW || selector == tmOARC) return 1;

    return 1; /* default */
}

/* ============================================================
 * TMPL parser — reads a generic template, stores as NODE_TMPL
 * (Phase 4+ will map to typed nodes; for now keep compatible)
 * ============================================================ */
NodePtr MtefParser::parseTmpl(int options) {
    if (options & OPT_NUDGE) skipNudge();

    int selector = readByte();
    if (selector < 0) selector = 0;
    int variation = readByte();
    if (variation < 0) variation = 0;

    /* Extended variation (>= 128) */
    if (variation >= 128) {
        int v2 = readByte();
        if (v2 >= 0) variation = ((variation - 128) << 8) | v2;
    }

    int numSlots = getSlotCount(selector, variation);

    /* For now, store as a generic template with LineNode slots.
     * The selector and variation are stored in a FenceNode/FracNode/etc.
     * wrapper that we'll create later. For Phase 3, we use a temporary
     * "GenericTmplNode" approach: store in a LineNode with metadata.
     * Actually, let's create the appropriate typed node right here. */

    /* --- Map selector to typed node --- */

    /* Fences */
    if (selector >= tmANGLE && selector <= tmLPRB) {
        auto n = std::make_unique<FenceNode>();
        n->selector = selector;
        n->variation = variation;
        /* Read slot[0] content */
        n->content = parseObjectList();
        return n;
    }

    /* Fractions */
    if (selector == tmFRACT) {
        auto n = std::make_unique<FracNode>();
        auto slot0 = parseObjectList();
        auto slot1 = parseObjectList();
        if (slot0.empty() && !slot1.empty()) {
            /* EQNEDT32 display fraction: slot[0]=empty, slot[1] has
             * both numer and denom as sub-LINEs. Split at the first LINE
             * boundary: first LINE = numer, second LINE = denom. */
            n->display = true;
            bool foundFirst = false;
            for (auto& item : slot1) {
                if (!foundFirst) {
                    n->numer.push_back(std::move(item));
                    foundFirst = true;
                } else {
                    n->denom.push_back(std::move(item));
                }
            }
        } else {
            n->numer = std::move(slot0);
            n->denom = std::move(slot1);
        }
        return n;
    }
    if (selector == tmSLFRACT) {
        auto n = std::make_unique<FracNode>();
        n->slashed = true;
        auto slot0 = parseObjectList();
        auto slot1 = parseObjectList();
        n->numer = std::move(slot0);
        n->denom = std::move(slot1);
        return n;
    }

    /* Root (sqrt) */
    if (selector == tmROOT) {
        auto n = std::make_unique<SqrtNode>();
        if (variation == 0) {
            /* var=0: 2 slots. slot[0]=empty(END), slot[1]=content LINE.
             * tex2mtef emits: TMPL(ROOT,0) END LINE(content) END [REC_END]
             * The content is in slot[1]. No index. */
            auto slot0 = parseObjectList();
            auto slot1 = parseObjectList();

            bool slot0Empty = slot0.empty();
            if (!slot0Empty) {
                slot0Empty = true;
                for (auto& c : slot0) {
                    if (c->tag() != Node::kLine) { slot0Empty = false; break; }
                    if (!static_cast<LineNode*>(c.get())->isNull) { slot0Empty = false; break; }
                }
            }

            if (slot0Empty) {
                /* slot[0]=empty: content is in slot[1] */
                /* Extract LINE children as content */
                for (auto& item : slot1) {
                    if (item->tag() == Node::kLine) {
                        auto* ln = static_cast<LineNode*>(item.get());
                        for (auto& c : ln->children) n->content.push_back(std::move(c));
                    } else {
                        n->content.push_back(std::move(item));
                    }
                }
                n->hasIndex = false;
            } else {
                n->content = std::move(slot0);
                n->hasIndex = !slot1.empty();
                n->index = std::move(slot1);
            }
        } else {
            /* var=1: 1 slot = content only. Index follows as sibling (SIZE_SUB + LINE).
             * tex2mtef: TMPL(ROOT,1) END LINE(content) END SIZE_SUB LINE(index) END */
            auto slot0 = parseObjectList();
            bool slot0Empty = slot0.empty();
            if (slot0Empty) {
                /* Empty slot: content is next LINE (handled by FenceMergePass) */
                n->hasIndex = true; /* index will come from sibling */
            } else {
                /* Content in slot[0] */
                for (auto& item : slot0) {
                    if (item->tag() == Node::kLine) {
                        auto* ln = static_cast<LineNode*>(item.get());
                        for (auto& c : ln->children) n->content.push_back(std::move(c));
                    } else {
                        n->content.push_back(std::move(item));
                    }
                }
                n->hasIndex = true;
            }
        }
        return n;
    }

    /* Scripts */
    if (selector == tmSCRIPT || selector == tmLSCRIPT) {
        auto n = std::make_unique<ScriptNode>();
        auto slot0 = parseObjectList();
        auto slot1 = parseObjectList();
        n->base = std::move(slot0);

        /* Split slot1 into sub and sup based on variation.
         * slot1 contains non-null LINEs (content) and null LINEs (placeholders).
         * variation: 0=sup-only, 1=sub-only, 2=both sub+sup */
        std::vector<NodePtr> nonNullLines;
        for (auto& item : slot1) {
            if (item->tag() == Node::kLine) {
                auto* ln = static_cast<LineNode*>(item.get());
                if (!ln->isNull && !ln->children.empty())
                    nonNullLines.push_back(std::move(item));
            }
            /* SIZE nodes in slot1 are ignored */
        }

        if (variation == 2 && nonNullLines.size() >= 2) {
            /* Both: first non-null = sub, second = sup */
            n->hasSub = true; n->hasSup = true;
            auto* subLn = static_cast<LineNode*>(nonNullLines[0].get());
            for (auto& c : subLn->children) n->sub.push_back(std::move(c));
            auto* supLn = static_cast<LineNode*>(nonNullLines[1].get());
            for (auto& c : supLn->children) n->sup.push_back(std::move(c));
        } else if (variation == 1 && !nonNullLines.empty()) {
            /* Sub only */
            n->hasSub = true;
            auto* subLn = static_cast<LineNode*>(nonNullLines[0].get());
            for (auto& c : subLn->children) n->sub.push_back(std::move(c));
        } else if (variation == 0 && !nonNullLines.empty()) {
            /* Sup only */
            n->hasSup = true;
            auto* supLn = static_cast<LineNode*>(nonNullLines[0].get());
            for (auto& c : supLn->children) n->sup.push_back(std::move(c));
        } else if (!nonNullLines.empty()) {
            /* Fallback: first line = sub */
            n->hasSub = true;
            auto* subLn = static_cast<LineNode*>(nonNullLines[0].get());
            for (auto& c : subLn->children) n->sub.push_back(std::move(c));
        }
        return n;
    }

    /* Decorations */
    if (selector == tmUBAR || selector == tmOBAR ||
        selector == tmLARROW || selector == tmRARROW || selector == tmBARROW) {
        auto n = std::make_unique<DecorationNode>();
        n->selector = selector;
        n->variation = variation;
        n->content = parseObjectList();
        return n;
    }

    /* Integrals */
    if (selector >= tmSINT && selector <= tmTSINT) {
        auto n = std::make_unique<IntegralNode>();
        n->selector = selector;
        n->variation = variation;
        /* Slots depend on variation */
        n->body = parseObjectList();  /* slot[0] always */
        if (variation & 0x01) {
            n->hasLower = true;
            n->lower = parseObjectList();
        }
        if (variation & 0x02) {
            n->hasUpper = true;
            n->upper = parseObjectList();
        }
        return n;
    }

    /* Horizontal braces */
    if (selector == tmUHBRACE || selector == tmLHBRACE) {
        auto n = std::make_unique<BraceDecoNode>();
        n->selector = selector;
        auto slot0 = parseObjectList();
        auto slot1 = parseObjectList();
        if (slot0.empty() && !slot1.empty()) {
            /* EQNEDT32 pattern: slot[0]=empty, slot[1] has content+label.
             * First non-null LINE = content, rest = label + display data */
            bool foundContent = false;
            for (auto& item : slot1) {
                if (!foundContent && item->tag() == Node::kLine) {
                    auto* ln = static_cast<LineNode*>(item.get());
                    if (!ln->isNull && !ln->children.empty()) {
                        for (auto& c : ln->children) n->content.push_back(std::move(c));
                        foundContent = true;
                        continue;
                    }
                }
                if (foundContent && item->tag() == Node::kLine) {
                    auto* ln = static_cast<LineNode*>(item.get());
                    if (!ln->isNull && !ln->children.empty()) {
                        for (auto& c : ln->children) n->label.push_back(std::move(c));
                        continue;
                    }
                }
                /* Skip SIZE and display chars */
            }
        } else {
            n->content = std::move(slot0);
            n->label = std::move(slot1);
        }
        return n;
    }

    /* BigOps (Sum, Product, etc.) */
    if (selector >= tmSUM && selector <= tmIINTER) {
        auto n = std::make_unique<BigOpNode>();
        n->selector = selector;
        n->body = parseObjectList();
        return n;
    }

    /* Lim */
    if (selector == tmLIM) {
        auto n = std::make_unique<LimNode>();
        n->content = parseObjectList();
        /* Additional slots for variation 2 */
        if (variation == 2) {
            parseObjectList(); /* slot[1] */
            parseObjectList(); /* slot[2] */
        } else {
            parseObjectList(); /* slot[1] */
        }
        return n;
    }

    /* Dirac */
    if (selector == tmDIRAC) {
        auto n = std::make_unique<DiracNode>();
        n->variation = variation;
        if (variation == 0) {
            n->bra = parseObjectList();
            n->ket = parseObjectList();
        } else {
            n->bra = parseObjectList();
        }
        return n;
    }

    /* Arrows (over/under) */
    if (selector == tmUARROW || selector == tmOARROW || selector == tmOARC) {
        auto n = std::make_unique<DecorationNode>();
        n->selector = selector;
        n->variation = variation;
        n->content = parseObjectList();
        return n;
    }

    /* Default: generic — read slots and discard */
    for (int i = 0; i < numSlots; i++)
        parseObjectList();
    return nullptr;
}

/* ============================================================
 * PILE parser
 * ============================================================ */
std::unique_ptr<PileNode> MtefParser::parsePile(int options) {
    if (options & OPT_NUDGE) skipNudge();

    bool hasRuler = (options & 0x02) != 0;

    int halign = readByte();
    if (halign < 0) halign = 0;

    if (hasRuler) {
        int nStops = readByte();
        if (nStops > 0) {
            for (int i = 0; i < nStops; i++) {
                readByte();   /* tab stop type */
                readUint16(); /* tab stop offset */
            }
        }
    }

    auto pile = std::make_unique<PileNode>();
    pile->halign = halign;

    /* Read LINE records until double-END or structural break */
    int sawEnd = 0;
    while (pos_ < len_) {
        int tag = readByte();
        if (tag < 0) break;
        int recType = tag & 0x0F;
        int opts = (tag >> 4) & 0x0F;

        if (recType == REC_END) {
            if (pos_ < len_) {
                int nextByte = data_[pos_];
                int nextRType = nextByte & 0x0F;
                if (nextRType == REC_END) break;
                if (nextRType == REC_LINE) { sawEnd = 1; continue; }
                if (!pile->lines.empty() &&
                    nextRType >= REC_CHAR && nextRType <= REC_FONT) {
                    break;
                }
            } else {
                break;
            }
            /* Lookahead heuristic for edge cases */
            int hasMore = 0;
            size_t scanEnd = pos_ + 256;
            if (scanEnd > len_) scanEnd = len_;
            for (size_t j = pos_; j < scanEnd; j++) {
                if ((data_[j] & 0x0F) == REC_LINE) { hasMore = 1; break; }
            }
            if (!hasMore) break;
            sawEnd = 1;
            continue;
        }

        if (recType == REC_LINE) {
            auto line = parseLine(opts);
            if (line) pile->lines.push_back(std::move(line));
        } else {
            /* Non-LINE record in PILE → parse as record and attach to last line */
            auto node = parseRecord(recType, opts);
            if (node && !pile->lines.empty()) {
                auto& lastLine = pile->lines.back();
                auto* ln = static_cast<LineNode*>(lastLine.get());
                ln->children.push_back(std::move(node));
            }
        }
    }

    return pile;
}

/* ============================================================
 * MATRIX parser
 * ============================================================ */
std::unique_ptr<MatrixNode> MtefParser::parseMatrix(int options) {
    if (options & OPT_NUDGE) skipNudge();

    auto n = std::make_unique<MatrixNode>();

    int valign = readByte();  /* vertical alignment */
    int hjust = readByte();   /* horizontal justification */
    int vjust = readByte();   /* vertical justification */
    n->rows = readByte();
    n->cols = readByte();

    if (n->rows < 0) n->rows = 0;
    if (n->cols < 0) n->cols = 0;

    /* Skip the row/column partition lines.  There are (rows + 1) of them and
     * (cols + 1) of them, at TWO BITS each, packed into bytes -- not one byte
     * apiece.  Read as bytes, a 2x1 matrix consumed five where the file has
     * two, and the three it took past them were its own first cell: the cells
     * came out empty with their contents left outside the matrix. */
    const int rowPartBytes = ((n->rows + 1) * 2 + 7) / 8;
    for (int i = 0; i < rowPartBytes; i++) readByte();
    const int colPartBytes = ((n->cols + 1) * 2 + 7) / 8;
    for (int i = 0; i < colPartBytes; i++) readByte();

    /* Read elements (row-major order) */
    int total = n->rows * n->cols;
    for (int i = 0; i < total; i++) {
        auto line = std::make_unique<LineNode>();
        line->children = parseObjectList();
        n->elements.push_back(std::move(line));
    }

    return n;
}

/* ============================================================
 * Public API: parse MTEF binary
 * ============================================================ */
MtefParser::Result MtefParser::parse(const uint8_t* data, size_t len) {
    Result result;
    if (!data || len < 5) return result;

    MtefParser parser(data, len);

    /* Read 5-byte header */
    int mtefVer = parser.readByte();    /* MTEF version (3) */
    int platform = parser.readByte();   /* platform (1=Windows) */
    int product = parser.readByte();    /* product (1=EqEd) */
    int prodMaj = parser.readByte();    /* product major */
    int prodMin = parser.readByte();    /* product minor */

    if (mtefVer < 0) return result;
    result.prodVer = prodMaj;

    /* Parse the body as a root LINE's children.
     *
     * Not ONE object list: Equation Editor closes a list with END and carries
     * on, so a document is a run of them.  Taking the first and stopping read
     * 86 bytes of an 1805-byte lecture file and returned a short equation
     * rather than a truncated one -- a wrong answer with no error in it. */
    auto root = std::make_unique<LineNode>();
    while (parser.pos_ < len) {
        const size_t before = parser.pos_;
        NodeList more = parser.parseObjectList();
        for (auto& n : more)
            if (n) root->children.push_back(std::move(n));
        if (parser.pos_ == before) break;   /* no progress: stop, do not spin */
    }
    result.root = std::move(root);

    return result;
}

} /* namespace mtef */
