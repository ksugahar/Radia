/*
 * latex_emitter.h -- Node tree → LaTeX text emitter
 *
 * Phase 4: Multi-pass conversion pipeline that transforms the unified
 * node tree (from MtefParser) into LaTeX text. Ports the 7-pass
 * convert_line system from mtef2tex.c.
 */
#ifndef LATEX_EMITTER_H
#define LATEX_EMITTER_H

#include "mtef_node.h"
#include "line_pass.h"
#include <string>
#include <memory>

namespace mtef {

class LaTeXEmitter {
public:
    /* run_passes repairs EQNEDT32's sibling layout and belongs only to trees
     * that came from MTEF.  A tree from the LaTeX parser already has every
     * slot filled; running the passes on it would move content that is
     * already where it belongs. */
    explicit LaTeXEmitter(int prodVer = 3, bool run_passes = true);

    /* Main entry: node tree → LaTeX string */
    std::string emit(const LineNode& root);

private:
    int prodVer_;
    /* How deep in fractions the walk is, so a fraction goes out at the size
     * it was DRAWN: \dfrac outermost, rac inside.  LaTeX's own rule then
     * agrees at every level and the paste matches the picture. */
    int fracDepth_ = 0;
    /* The style in force, so a marker that repeats it writes nothing. */
    int sizeStyle_ = 0;                 /* SIZETYPE_FULL */
    static const char* arrow_command(unsigned cp);
    bool runPasses_ = true;
    int depth_ = 0;
    PassPipeline pipeline_;

    /* Per-node emission (recursive) */
    void emitNode(const Node* node, std::string& out);
    void emitLine(const LineNode& line, std::string& out);
    void emitChar(const CharNode& ch, std::string& out);
    void emitFence(const FenceNode& fence, std::string& out);
    void emitFrac(const FracNode& frac, std::string& out);
    void emitSqrt(const SqrtNode& sq, std::string& out);
    void emitScript(const ScriptNode& script, std::string& out);
    void emitIntegral(const IntegralNode& integ, std::string& out);
    void emitBigOp(const BigOpNode& bigop, std::string& out);
    void emitDecoration(const DecorationNode& deco, std::string& out);
    void emitEmbell(const EmbellNode& em, std::string& out);
    void emitBraceDeco(const BraceDecoNode& bd, std::string& out);
    void emitDirac(const DiracNode& dirac, std::string& out);
    void emitLim(const LimNode& lim, std::string& out);
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

} /* namespace mtef */

#endif /* LATEX_EMITTER_H */
