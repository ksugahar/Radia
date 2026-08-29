/* Structural syntax tree used by the TeX parser, editor, and renderer. */
#ifndef EQUATION_NODE_H
#define EQUATION_NODE_H

#include "equation_types.h"
#include <stdint.h>
#include <memory>
#include <vector>
#include <string>

namespace eqnedit {

/* Parser, editor, layout, and emit all recurse over the structural tree.
 * Keep one shared limit so no editing route can create a tree that the parser
 * would truncate when Undo restores its LaTeX snapshot. */
inline constexpr int kMaxNestingDepth = 200;

/* ============================================================
 * Forward declarations
 * ============================================================ */
class Node;
class NodeVisitor;

using NodeList = std::vector<std::unique_ptr<Node>>;
using NodePtr  = std::unique_ptr<Node>;

/* Slots of a node in visual order.  Empty for a leaf.  Parser and editor both
 * use this traversal so a math alphabet reaches nested fractions, scripts,
 * matrices, and decorations consistently. */
std::vector<NodeList*> node_slots(Node& n);

/* ============================================================
 * Base Node class
 * ============================================================ */
class Node {
public:
    enum Tag : uint8_t {
        /* Structural */
        kLine, kPile, kMatrix,
        /* Character */
        kChar, kEmbell, kSize, kFont,
        /* Templates — fences */
        kFence,
        /* Templates — structural */
        kFrac, kSqrt, kScript,
        /* Templates — integrals & BigOps */
        kIntegral, kBigOp,
        /* Templates — decorations */
        kDecoration, kBraceDeco,
        /* Templates — special */
        kDirac, kLim,
        /* LaTeX-specific (only produced by LaTeX parser) */
        kFunction, kText, kMathbf, kGroup, kPrime, kDegree, kOverset,
        /* Sentinel for RM char */
        kRM,
    };

    virtual ~Node() = default;
    Tag tag() const { return tag_; }
    virtual void accept(NodeVisitor& v) = 0;

protected:
    explicit Node(Tag t) : tag_(t) {}

private:
    Tag tag_;
};

/* ============================================================
 * Leaf / simple nodes
 * ============================================================ */

class LineNode : public Node {
public:
    NodeList children;
    bool isNull = false;

    LineNode() : Node(kLine) {}
    void accept(NodeVisitor& v) override;
};

class CharNode : public Node {
public:
    int typeface = 0;       /* logical Typeface value */
    /* True only for Function style inferred by Math input.  Explicit or
     * parsed Function style stays false, so extending an automatically
     * recognised word can safely revert just the inferred letters. */
    bool automaticFunction = false;
    /* A Unicode scalar value.  Keeping this at 16 bits corrupted UTF-8 input
     * byte-by-byte and made supplementary-plane characters impossible. */
    uint32_t charCode = 0;
    char ch = 0;            /* for ND_CHAR: ASCII char */
    std::string latex;      /* for ND_SYMBOL: LaTeX command (e.g. "\\nabla") */
    std::vector<int> embells; /* chained embellishment types */

    CharNode() : Node(kChar) {}
    void accept(NodeVisitor& v) override;
};

class EmbellNode : public Node {
public:
    int embellType = 0;
    NodeList content;       /* for ND_EMBELL (LaTeX: \hat{content}) */

    EmbellNode() : Node(kEmbell) {}
    void accept(NodeVisitor& v) override;
};

class SizeNode : public Node {
public:
    int sizeType = 0;       /* SIZETYPE_FULL..SIZETYPE_SUBSYM */

    SizeNode() : Node(kSize) {}
    void accept(NodeVisitor& v) override;
};

class FontNode : public Node {
public:
    int fontIndex = 0;
    int style = 0;
    char name[64] = {};

    FontNode() : Node(kFont) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Container nodes
 * ============================================================ */

class PileNode : public Node {
public:
    int halign = 0;         /* 0=gathered, 1=aligned, 20-24=matrix types */
    int kind = 0;           /* ND_ENVIRONMENT: 0=gathered,1=aligned,2..7=matrix variants */
    int ncols = 0;          /* for aligned/matrix */
    NodeList lines;         /* each is a LineNode */

    PileNode() : Node(kPile) {}
    void accept(NodeVisitor& v) override;
};

class MatrixNode : public Node {
public:
    /* kCasesLayout is a matrix whose cells are left-aligned, which is what
     * `cases` does.  Centring them made the canvas disagree with the PDF,
     * and it is also what tells a real `cases` apart from someone's literal
     * \left\{ \begin{matrix} ... \right. . */
    enum LayoutKind { kMatrixLayout = 0, kAlignedLayout = 1, kCasesLayout = 2 };
    int rows = 0, cols = 0;
    NodeList elements;      /* row-major LineNodes */
    LayoutKind layoutKind = kMatrixLayout;

    MatrixNode() : Node(kMatrix) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — fences
 * ============================================================ */

class FenceNode : public Node {
public:
    int selector = 0;       /* tmANGLE..tmLPRB: the opening delimiter */
    int variation = 0;      /* 0=both, 1=left-only, 2=right-only */
    /* Closing delimiter when it differs from the opening one, as in the
     * half-open interval \left( x \right].  -1 means "same as selector",
     * which is the common matched-pair case. */
    int rightSelector = -1;
    NodeList content;

    int right_selector() const {
        return rightSelector >= 0 ? rightSelector : selector;
    }

    FenceNode() : Node(kFence) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — structural
 * ============================================================ */

class FracNode : public Node {
public:
    bool slashed = false;   /* tmSLFRACT */
    bool display = false;   /* display fraction (Pass 0a) */
    NodeList numer;
    NodeList denom;

    FracNode() : Node(kFrac) {}
    void accept(NodeVisitor& v) override;
};

class SqrtNode : public Node {
public:
    NodeList content;
    NodeList index;         /* optional nth-root index */
    bool hasIndex = false;

    SqrtNode() : Node(kSqrt) {}
    void accept(NodeVisitor& v) override;
};

class ScriptNode : public Node {
public:
    NodeList base;
    NodeList sub;
    NodeList sup;
    bool hasSub = false;
    bool hasSup = false;

    ScriptNode() : Node(kScript) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — integrals & BigOps
 * ============================================================ */

class IntegralNode : public Node {
public:
    int selector = tmSINT;  /* tmSINT..tmTSINT, tmSSINT..tmTSINT */
    int variation = 0;
    NodeList body;
    NodeList lower;
    NodeList upper;
    bool hasLower = false;
    bool hasUpper = false;
    bool hasLimits = false;

    IntegralNode() : Node(kIntegral) {}
    void accept(NodeVisitor& v) override;
};

class BigOpNode : public Node {
public:
    int selector = tmSUM;   /* tmSUM..tmIINTER */
    NodeList body;
    NodeList lower;
    NodeList upper;
    bool hasLower = false;
    bool hasUpper = false;
    bool hasLimits = false;

    BigOpNode() : Node(kBigOp) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — decorations
 * ============================================================ */

class DecorationNode : public Node {
public:
    int selector = 0;       /* tmUBAR, tmOBAR, tmLARROW, tmRARROW, tmBARROW */
    int variation = 0;
    NodeList content;

    DecorationNode() : Node(kDecoration) {}
    void accept(NodeVisitor& v) override;
};

class BraceDecoNode : public Node {
public:
    int selector = 0;       /* tmUHBRACE or tmLHBRACE */
    NodeList content;
    NodeList label;

    BraceDecoNode() : Node(kBraceDeco) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — special
 * ============================================================ */

class DiracNode : public Node {
public:
    int variation = 0;
    NodeList bra;
    NodeList ket;

    DiracNode() : Node(kDirac) {}
    void accept(NodeVisitor& v) override;
};

class LimNode : public Node {
public:
    NodeList content;

    LimNode() : Node(kLim) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * TeX command and grouping nodes
 * ============================================================ */

class FunctionNode : public Node {
public:
    char name[64] = {};
    bool isOperator = false;

    FunctionNode() : Node(kFunction) {}
    void accept(NodeVisitor& v) override;
};

class TextNode : public Node {
public:
    char text[256] = {};

    TextNode() : Node(kText) {}
    void accept(NodeVisitor& v) override;
};

class MathbfNode : public Node {
public:
    NodeList content;

    MathbfNode() : Node(kMathbf) {}
    void accept(NodeVisitor& v) override;
};

class GroupNode : public Node {
public:
    NodeList children;

    GroupNode() : Node(kGroup) {}
    void accept(NodeVisitor& v) override;
};

class PrimeNode : public Node {
public:
    int count = 1;

    PrimeNode() : Node(kPrime) {}
    void accept(NodeVisitor& v) override;
};

class DegreeNode : public Node {
public:
    DegreeNode() : Node(kDegree) {}
    void accept(NodeVisitor& v) override;
};

class OversetNode : public Node {
public:
    NodeList over;
    NodeList base;
    bool under = false;       /* false: \\overset, true: \\underset */

    OversetNode() : Node(kOverset) {}
    void accept(NodeVisitor& v) override;
};

class RMNode : public Node {
public:
    char ch = 0;

    RMNode() : Node(kRM) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Visitor interface
 * ============================================================ */

class NodeVisitor {
public:
    virtual ~NodeVisitor() = default;
    virtual void visit(LineNode&) = 0;
    virtual void visit(CharNode&) = 0;
    virtual void visit(EmbellNode&) = 0;
    virtual void visit(SizeNode&) = 0;
    virtual void visit(FontNode&) = 0;
    virtual void visit(PileNode&) = 0;
    virtual void visit(MatrixNode&) = 0;
    virtual void visit(FenceNode&) = 0;
    virtual void visit(FracNode&) = 0;
    virtual void visit(SqrtNode&) = 0;
    virtual void visit(ScriptNode&) = 0;
    virtual void visit(IntegralNode&) = 0;
    virtual void visit(BigOpNode&) = 0;
    virtual void visit(DecorationNode&) = 0;
    virtual void visit(BraceDecoNode&) = 0;
    virtual void visit(DiracNode&) = 0;
    virtual void visit(LimNode&) = 0;
    virtual void visit(FunctionNode&) = 0;
    virtual void visit(TextNode&) = 0;
    virtual void visit(MathbfNode&) = 0;
    virtual void visit(GroupNode&) = 0;
    virtual void visit(PrimeNode&) = 0;
    virtual void visit(DegreeNode&) = 0;
    virtual void visit(OversetNode&) = 0;
    virtual void visit(RMNode&) = 0;
};

/* ============================================================
 * Utility: NodeList helpers
 * ============================================================ */

/* Move a node into a NodeList */
inline void push(NodeList& list, NodePtr node) {
    list.push_back(std::move(node));
}

/* Create a typed node (convenience template) */
template<typename T, typename... Args>
std::unique_ptr<T> make_node(Args&&... args) {
    return std::make_unique<T>(std::forward<Args>(args)...);
}

} /* namespace eqnedit */

#endif /* EQUATION_NODE_H */
