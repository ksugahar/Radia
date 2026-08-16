/*
 * mtef_parser.h -- MTEF binary → unified node tree parser
 *
 * Phase 3: Parses MTEF v3 binary format into the unified node hierarchy
 * defined in mtef_node.h. Faithfully reproduces the behavior of
 * parse_record/parse_char/parse_tmpl/parse_line/parse_pile from mtef2tex.c.
 */
#ifndef MTEF_PARSER_H
#define MTEF_PARSER_H

#include "mtef_node.h"
#include <cstdint>
#include <cstddef>
#include <memory>

namespace mtef {

class MtefParser {
public:
    struct Result {
        std::unique_ptr<LineNode> root;
        int prodVer = 0;   /* 3 or 10 */
    };

    /* Parse MTEF binary (raw DS Equation data, with 5-byte header).
     * Returns root LineNode or nullptr on failure. */
    static Result parse(const uint8_t* data, size_t len);

private:
    MtefParser(const uint8_t* data, size_t len);

    /* Low-level byte reading */
    int readByte();
    int readUint16();       /* little-endian uint16 */
    void skipNudge();       /* skip 4-byte nudge data */

    /* Record dispatch */
    NodePtr parseRecord(int recType, int options);
    NodeList parseObjectList();

    /* Typed parsers */
    std::unique_ptr<LineNode>  parseLine(int options);
    std::unique_ptr<CharNode>  parseChar(int options);
    NodePtr                    parseTmpl(int options);
    std::unique_ptr<PileNode>  parsePile(int options);
    std::unique_ptr<MatrixNode> parseMatrix(int options);
    std::unique_ptr<EmbellNode> parseEmbell(int options);
    std::unique_ptr<SizeNode>   parseSize(int options);
    std::unique_ptr<FontNode>   parseFont(int options);

    /* Embellishment chain parser */
    std::vector<int> parseEmbellList();

    /* Template slot count calculator */
    static int getSlotCount(int selector, int variation);

    /* State */
    const uint8_t* data_;
    size_t len_;
    size_t pos_ = 0;
    int prodVer_ = 0;
};

} /* namespace mtef */

#endif /* MTEF_PARSER_H */
