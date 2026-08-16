/*
 * mtef_node.cpp -- NodeVisitor accept() implementations
 */
#include "mtef_node.h"

namespace mtef {

void LineNode::accept(NodeVisitor& v) { v.visit(*this); }
void CharNode::accept(NodeVisitor& v) { v.visit(*this); }
void EmbellNode::accept(NodeVisitor& v) { v.visit(*this); }
void SizeNode::accept(NodeVisitor& v) { v.visit(*this); }
void FontNode::accept(NodeVisitor& v) { v.visit(*this); }
void PileNode::accept(NodeVisitor& v) { v.visit(*this); }
void MatrixNode::accept(NodeVisitor& v) { v.visit(*this); }
void FenceNode::accept(NodeVisitor& v) { v.visit(*this); }
void FracNode::accept(NodeVisitor& v) { v.visit(*this); }
void SqrtNode::accept(NodeVisitor& v) { v.visit(*this); }
void ScriptNode::accept(NodeVisitor& v) { v.visit(*this); }
void IntegralNode::accept(NodeVisitor& v) { v.visit(*this); }
void BigOpNode::accept(NodeVisitor& v) { v.visit(*this); }
void DecorationNode::accept(NodeVisitor& v) { v.visit(*this); }
void BraceDecoNode::accept(NodeVisitor& v) { v.visit(*this); }
void DiracNode::accept(NodeVisitor& v) { v.visit(*this); }
void LimNode::accept(NodeVisitor& v) { v.visit(*this); }
void FunctionNode::accept(NodeVisitor& v) { v.visit(*this); }
void TextNode::accept(NodeVisitor& v) { v.visit(*this); }
void MathbfNode::accept(NodeVisitor& v) { v.visit(*this); }
void GroupNode::accept(NodeVisitor& v) { v.visit(*this); }
void PrimeNode::accept(NodeVisitor& v) { v.visit(*this); }
void DegreeNode::accept(NodeVisitor& v) { v.visit(*this); }
void OversetNode::accept(NodeVisitor& v) { v.visit(*this); }
void RMNode::accept(NodeVisitor& v) { v.visit(*this); }

} /* namespace mtef */
