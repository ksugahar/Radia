/*
 * mtef_dump.cpp -- indented text dump of the parsed node tree.
 *
 * The LaTeX emitter, the SVG renderer and the OMML emitter all traverse the
 * same tree, and every disagreement between them so far has come from reading
 * a different field of it rather than from the parse.  Printing the tree is
 * the fastest way to settle which field actually carries the content, so this
 * is a permanent diagnostic rather than scaffolding.
 *
 * Optionally runs the LINE pass pipeline first, which is what the emitters
 * see -- the raw parse and the post-pass tree differ substantially (fences and
 * big-operator limits only move into their slots during the passes).
 */
#include "mtef_dump.h"
#include "mtef_parser.h"
#include "line_pass.h"
#include "mtef_common.h"

#include <sstream>

namespace mtef {
namespace {

class Dumper {
public:
    Dumper(bool run_passes, int prodVer)
        : run_passes_(run_passes), prodVer_(prodVer) {}

    std::string run(const LineNode& root) {
        line(root, 0, "root");
        return o_.str();
    }

private:
    std::ostringstream o_;
    bool run_passes_;
    int prodVer_;
    PassPipeline pipeline_;

    void pad(int d) { for (int i = 0; i < d; ++i) o_ << "  "; }

    void slot(const NodeList& l, int d, const char* name) {
        if (l.empty()) return;
        pad(d); o_ << "." << name << ":\n";
        for (const auto& n : l) node(n.get(), d + 1);
    }

    void line(const LineNode& l, int d, const char* name) {
        pad(d); o_ << "LINE " << name << (l.isNull ? " [null]" : "")
                   << " n=" << l.children.size() << "\n";
        if (l.isNull) return;
        NodeList& kids = const_cast<NodeList&>(l.children);
        if (run_passes_) pipeline_.process(kids, d, prodVer_);
        for (const auto& n : kids) node(n.get(), d + 1);
    }

    void node(const Node* n, int d) {
        if (!n) { pad(d); o_ << "(null)\n"; return; }
        switch (n->tag()) {
            case Node::kLine:
                line(static_cast<const LineNode&>(*n), d, "");
                break;
            case Node::kChar: {
                const auto& c = static_cast<const CharNode&>(*n);
                pad(d); o_ << "CHAR tf=" << c.typeface
                           << " code=0x" << std::hex << c.charCode << std::dec;
                if (c.ch) o_ << " ch='" << c.ch << "'";
                if (!c.latex.empty()) o_ << " latex=" << c.latex;
                o_ << "\n";
                break;
            }
            case Node::kScript: {
                const auto& s = static_cast<const ScriptNode&>(*n);
                pad(d); o_ << "SCRIPT sub=" << s.hasSub << " sup=" << s.hasSup
                           << " base_n=" << s.base.size() << "\n";
                slot(s.base, d + 1, "base");
                slot(s.sub,  d + 1, "sub");
                slot(s.sup,  d + 1, "sup");
                break;
            }
            case Node::kFrac: {
                const auto& f = static_cast<const FracNode&>(*n);
                pad(d); o_ << "FRAC slashed=" << f.slashed << "\n";
                slot(f.numer, d + 1, "numer");
                slot(f.denom, d + 1, "denom");
                break;
            }
            case Node::kSqrt: {
                const auto& s = static_cast<const SqrtNode&>(*n);
                pad(d); o_ << "SQRT index=" << s.hasIndex << "\n";
                slot(s.index, d + 1, "index");
                slot(s.content, d + 1, "content");
                break;
            }
            case Node::kFence: {
                const auto& f = static_cast<const FenceNode&>(*n);
                pad(d); o_ << "FENCE sel=" << f.selector
                           << " var=" << f.variation << "\n";
                slot(f.content, d + 1, "content");
                break;
            }
            case Node::kIntegral: {
                const auto& i = static_cast<const IntegralNode&>(*n);
                pad(d); o_ << "INTEGRAL sel=" << i.selector
                           << " var=" << i.variation
                           << " hasLower=" << i.hasLower
                           << " hasUpper=" << i.hasUpper
                           << " limits=" << i.hasLimits
                           << " dispLo=" << (i.displayLower ? 1 : 0)
                           << " dispHi=" << (i.displayUpper ? 1 : 0) << "\n";
                slot(i.lower, d + 1, "lower");
                slot(i.upper, d + 1, "upper");
                slot(i.body,  d + 1, "body");
                if (i.displayLower) line(*i.displayLower, d + 1, "displayLower");
                if (i.displayUpper) line(*i.displayUpper, d + 1, "displayUpper");
                break;
            }
            case Node::kBigOp: {
                const auto& b = static_cast<const BigOpNode&>(*n);
                pad(d); o_ << "BIGOP sel=" << b.selector
                           << " hasLower=" << b.hasLower
                           << " hasUpper=" << b.hasUpper
                           << " limits=" << b.hasLimits
                           << " dispLo=" << (b.displayLower ? 1 : 0)
                           << " dispHi=" << (b.displayUpper ? 1 : 0) << "\n";
                slot(b.lower, d + 1, "lower");
                slot(b.upper, d + 1, "upper");
                slot(b.body,  d + 1, "body");
                if (b.displayLower) line(*b.displayLower, d + 1, "displayLower");
                if (b.displayUpper) line(*b.displayUpper, d + 1, "displayUpper");
                break;
            }
            case Node::kPile: {
                const auto& p = static_cast<const PileNode&>(*n);
                pad(d); o_ << "PILE n=" << p.lines.size()
                           << " ncols=" << p.ncols
                           << " halign=" << p.halign << "\n";
                for (const auto& l : p.lines) node(l.get(), d + 1);
                break;
            }
            case Node::kMatrix: {
                const auto& m = static_cast<const MatrixNode&>(*n);
                pad(d); o_ << "MATRIX " << m.rows << "x" << m.cols << "\n";
                for (const auto& e : m.elements) node(e.get(), d + 1);
                break;
            }
            case Node::kEmbell: {
                const auto& e = static_cast<const EmbellNode&>(*n);
                pad(d); o_ << "EMBELL type=" << e.embellType << "\n";
                slot(e.content, d + 1, "content");
                break;
            }
            case Node::kDecoration: {
                const auto& dn = static_cast<const DecorationNode&>(*n);
                pad(d); o_ << "DECO sel=" << dn.selector << "\n";
                slot(dn.content, d + 1, "content");
                break;
            }
            case Node::kBraceDeco: {
                const auto& b = static_cast<const BraceDecoNode&>(*n);
                pad(d); o_ << "BRACEDECO sel=" << b.selector << "\n";
                slot(b.content, d + 1, "content");
                slot(b.label, d + 1, "label");
                break;
            }
            case Node::kOverset: {
                const auto& o = static_cast<const OversetNode&>(*n);
                pad(d); o_ << "OVERSET under=" << (o.under ? 1 : 0) << "\n";
                slot(o.over, d + 1, "over");
                slot(o.base, d + 1, "base");
                break;
            }
            case Node::kPhantom: {
                const auto& p = static_cast<const PhantomNode&>(*n);
                pad(d); o_ << "PHANTOM w=" << (p.keepWidth ? 1 : 0)
                           << " h=" << (p.keepHeight ? 1 : 0) << "\n";
                slot(p.content, d + 1, "content");
                break;
            }
            case Node::kSize:
                pad(d); o_ << "SIZE\n"; break;
            case Node::kFont:
                pad(d); o_ << "FONT\n"; break;
            default:
                pad(d); o_ << "NODE tag=" << int(n->tag()) << "\n"; break;
        }
    }
};

}  // namespace

std::string dump_tree(const uint8_t* data, size_t len, bool run_passes) {
    MtefParser::Result res = MtefParser::parse(data, len);
    if (!res.root) return "(parse failed)\n";
    Dumper d(run_passes, res.prodVer);
    std::ostringstream head;
    head << "prodVer=" << res.prodVer << " passes=" << (run_passes ? "on" : "off") << "\n";
    return head.str() + d.run(*res.root);
}

std::string dump_latex_tree(const LineNode& root) {
    Dumper d(/*run_passes=*/false, 3);
    return d.run(root);
}

}  // namespace mtef
