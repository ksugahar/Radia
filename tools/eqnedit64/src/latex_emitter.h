/* TeX syntax tree -> normalized TeX text. */
#ifndef LATEX_EMITTER_H
#define LATEX_EMITTER_H

#include "equation_node.h"
#include <string>
#include <memory>

namespace eqnedit {

class LaTeXEmitter {
public:
    LaTeXEmitter() = default;

    /* Main entry: node tree → LaTeX string */
    std::string emit(const LineNode& root);

    /* Emit a contiguous range from one editor slot.  This is used for
     * structural copy/cut: no temporary clone of the unique_ptr tree is
     * required, and the same serializer is used for documents and
     * selections. */
    std::string emit_range(const NodeList& nodes, size_t first, size_t last);

private:
    /* Per-node emission (recursive) */
    void emitNode(const Node* node, std::string& out);
    void emitLine(const LineNode& line, std::string& out);
    void emitSequence(const NodeList& nodes, size_t first, size_t last,
                      std::string& out);
    void emitChar(const CharNode& ch, std::string& out);
    void emitFence(const FenceNode& fence, std::string& out);
    void emitFrac(const FracNode& frac, std::string& out);
    void emitSqrt(const SqrtNode& sq, std::string& out);
    void emitScript(const ScriptNode& script, std::string& out);
    void emitIntegral(const IntegralNode& integ, std::string& out);
    void emitBigOp(const BigOpNode& bigop, std::string& out);
    void emitDecoration(const DecorationNode& deco, std::string& out);
    void emitEmbell(const EmbellNode& embell, std::string& out);
    void emitBraceDeco(const BraceDecoNode& bd, std::string& out);
    void emitDirac(const DiracNode& dirac, std::string& out);
    void emitLim(const LimNode& lim, std::string& out);
    void emitOverset(const OversetNode& stacked, std::string& out);
    void emitPile(const PileNode& pile, std::string& out);
    void emitMatrix(const MatrixNode& mat, std::string& out);

    /* Helper: emit a NodeList to string */
    std::string emitNodes(const NodeList& nodes);

    /* Character lookup */
    const char* lookupSymbol(uint16_t code) const;
    const char* lookupGreekLower(uint16_t code) const;
    const char* lookupGreekUpper(uint16_t code) const;

    /* Post-processing */
    static std::string postProcess(const std::string& raw);
};

/* Serialise a tree the LaTeX parser built back to LaTeX.  This is what saves
 * an edited equation: the stored form of an equation is its LaTeX. */
std::string tree_to_latex(const LineNode& root);

} /* namespace eqnedit */

#endif /* LATEX_EMITTER_H */
