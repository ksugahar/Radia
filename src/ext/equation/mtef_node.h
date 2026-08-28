/*
 * mtef_node.h -- structural equation node hierarchy
 *
 * The filename and namespace are retained for source compatibility with the
 * mature layout engine. TeX is now the sole producer and supported source.
 */
#ifndef MTEF_NODE_H
#define MTEF_NODE_H

#include "mtef_common.h"
#include <stdint.h>
#include <memory>
#include <vector>
#include <string>

namespace mtef {

/* ============================================================
 * Forward declarations
 * ============================================================ */
class Node;
class NodeVisitor;

using NodeList = std::vector<std::unique_ptr<Node>>;
using NodePtr  = std::unique_ptr<Node>;

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
        kDirac, kLim, kPhantom,
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
    int typeface = 0;       /* logical typeface (1-11) */
    uint16_t charCode = 0;
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
    /* 0=centred (gathered), 1=aligned at &, 2=flush left, 3=flush right,
     * 20-24=matrix types.  The Format menu's Align Left / Center / Right. */
    int halign = 0;
    int kind = 0;           /* ND_ENVIRONMENT: 0=gathered,1=aligned,2..7=matrix variants */
    int ncols = 0;          /* for aligned/matrix */
    NodeList lines;         /* each is a LineNode */

    PileNode() : Node(kPile) {}
    void accept(NodeVisitor& v) override;
};

class MatrixNode : public Node {
public:
    int rows = 0, cols = 0;
    NodeList elements;      /* row-major LineNodes */

    MatrixNode() : Node(kMatrix) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — fences
 * ============================================================ */

class FenceNode : public Node {
public:
    int selector = 0;       /* tmANGLE..tmLPRB */
    int variation = 0;      /* 0=both, 1=left-only, 2=right-only */
    NodeList content;
    /* A bra-ket has a bar down the middle that grows with the fence.  TeX
     * spells it \middle|, and it is a third DELIMITER, not a character in the
     * content -- a character would not grow. */
    bool hasMiddle = false;
    int middle = 0;         /* the delimiter selector for it */
    NodeList content2;

    FenceNode() : Node(kFence) {}
    void accept(NodeVisitor& v) override;
};

/* ============================================================
 * Template nodes — structural
 * ============================================================ */

class FracNode : public Node {
public:
    bool slashed = false;   /* tmSLFRACT */
    /* False for inom: the parts stack with no rule between them. */
    bool ruled = true;
    bool display = false;   /* display fraction (Pass 0a) */
    /* Which size the author ASKED for, as opposed to the one the nesting
     * depth implies: 0 = by depth (\frac), +1 = \dfrac, -1 = \tfrac.
     * Without it \tfrac drew display-size at the top level, which is a
     * quiet wrong answer rather than a missing feature. */
    int styleOverride = 0;
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
    /* Scripts on the LEFT: the 14 and the 6 of a carbon isotope.  TeX writes
     * them as an empty-based script, {}^{14}_{6}C, and they hang off the
     * atom that FOLLOWS rather than the one before. */
    bool pre = false;

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
    bool isEcho = false;
    /* Display data (populated by mtef2tex passes) */
    std::unique_ptr<LineNode> displayLower;
    std::unique_ptr<LineNode> displayUpper;

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
    /* Display data (populated by mtef2tex passes) */
    std::unique_ptr<LineNode> displayLower;
    std::unique_ptr<LineNode> displayUpper;

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

/* \phantom: the room something would take, with nothing in it.
 *
 * Used to line one line of an equation up with another, and to reserve space
 * beside a tall thing so a script clears it.  \hphantom keeps only the width
 * and \vphantom only the height, which is what each is for. */
class PhantomNode : public Node {
public:
    NodeList content;
    bool keepWidth = true;
    bool keepHeight = true;

    PhantomNode() : Node(kPhantom) {}
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
 * LaTeX-specific nodes (only from LaTeX parser, not MTEF parser)
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
    bool under = false;     /* \underset puts it below instead */

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
    virtual void visit(PhantomNode&) = 0;
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

/* The function name a list spells, or empty when it is not one.
 *
 * A name arrives as a run of TF_FUNCTION characters -- the parser builds
 * "sin" as three of them -- so reading it back is how both the layout and the
 * writers recognise \lim without a node type of its own, and why they cannot
 * disagree about it. */
std::string function_name_of(const NodeList& list);

} /* namespace mtef */

#endif /* MTEF_NODE_H */
