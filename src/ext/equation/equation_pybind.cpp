/*
 * mtef_pybind.cpp -- Python binding for the MTEF codec + SVG renderer.
 *
 * This is the ONLY supported programmatic entry point.  The former stdin/stdout
 * CLIs were test scaffolding and are retired: one C++ source of truth, one
 * binding, and the Python tests import this module instead of spawning
 * processes.
 *
 * Architecture note: the codec and renderer are architecture-neutral and build
 * as x64 here to match CPython.  Only mtef2tex_hook.dll stays x86, because it
 * is loaded into EQNEDT32 -- which from now on is the GUI a human uses, not
 * something a program starts.
 */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <cstdlib>
#include <cstring>
#include <fstream>
#ifdef _WIN32
#include <windows.h>
#endif
#include <stdexcept>
#include <map>
#include <string>
#include <vector>

#include "mtef2tex.h"
#include "tex2mtef.h"
#include "mtef_svg.h"
#include "mtef_parser.h"
#include "mtef_omml.h"
#include "mtef_rtf.h"
#include "mtef_mathml.h"
#include "gvml_clip.h"
#include "mtef_gdi.h"
#include "mtef_dump.h"
#include "tex_parser.h"
#include "latex_emitter.h"
#include "eq_chords.h"
#include "math_font.h"
#include "eq_edit.h"
#include "md_doc.h"
#include "md_blocks.h"
#include "md_layout.h"

namespace py = pybind11;

namespace {

py::bytes tex_to_mtef_py(const std::string& latex) {
    int len = 0;
    uint8_t* out = tex_to_mtef(latex.c_str(), &len);
    if (!out || len <= 0) {
        if (out) free(out);
        throw std::runtime_error("tex_to_mtef: conversion produced no MTEF for: "
                                 + latex);
    }
    py::bytes result(reinterpret_cast<const char*>(out), size_t(len));
    free(out);
    return result;
}

std::string mtef_to_tex_py(const py::bytes& data) {
    std::string buf = data;
    char* out = mtef_to_latex_c(reinterpret_cast<const uint8_t*>(buf.data()),
                                buf.size());
    if (!out) throw std::runtime_error("mtef_to_tex: parse failed");
    std::string result(out);
    free(out);
    return result;
}

/* MTEF -> the editor's tree -> LaTeX.  mtef_to_tex answers with the legacy
 * converter's reading; this answers with the editor's, which is the one that
 * will be drawn and edited. */
std::string mtef_to_latex_py(const py::bytes& data) {
    std::string buf = data;
    mtef::MtefParser::Result res = mtef::MtefParser::parse(
        reinterpret_cast<const uint8_t*>(buf.data()), buf.size());
    if (!res.root) throw std::runtime_error("mtef_to_latex: parse failed");
    mtef::LaTeXEmitter em(res.prodVer, true);
    return em.emit(*res.root);
}

std::string mtef_to_svg_py(const py::bytes& data, const mtef::SvgStyle& style) {
    std::string buf = data;
    std::string svg = mtef::mtef_to_svg(
        reinterpret_cast<const uint8_t*>(buf.data()), buf.size(), style);
    if (svg.empty()) throw std::runtime_error("mtef_to_svg: parse failed");
    return svg;
}

std::string tex_to_svg_py(const std::string& latex, const mtef::SvgStyle& style) {
    std::string svg = mtef::tex_to_svg(latex, style);
    if (svg.empty()) throw std::runtime_error("tex_to_svg: render failed for: " + latex);
    return svg;
}

/* LaTeX -> tree -> LaTeX.  The tree is what an editor mutates, so this is the
 * shape an equation has after a round of editing; it must be stable. */
std::string tex_normalize_py(const std::string& latex) {
    std::unique_ptr<mtef::LineNode> root = mtef::parse_latex(latex);
    if (!root) throw std::runtime_error("parse_latex failed for: " + latex);
    return mtef::tree_to_latex(*root);
}

std::string mtef_to_omml_py(const py::bytes& data, const mtef::OmmlOptions& opt) {
    std::string buf = data;
    std::string omml = mtef::mtef_to_omml(
        reinterpret_cast<const uint8_t*>(buf.data()), buf.size(), opt);
    if (omml.empty()) throw std::runtime_error("mtef_to_omml: parse failed");
    return omml;
}

std::string tex_to_omml_py(const std::string& latex, const mtef::OmmlOptions& opt) {
    std::string omml = mtef::tex_to_omml(latex, opt);
    if (omml.empty()) throw std::runtime_error("tex_to_omml: emit failed for: " + latex);
    return omml;
}

std::string tex_to_rtf_py(const std::string& latex, const mtef::RtfOptions& opt) {
    std::string rtf = mtef::tex_to_rtf(latex, opt);
    if (rtf.empty()) throw std::runtime_error("tex_to_rtf: emit failed for: " + latex);
    return rtf;
}

std::string mtef_to_rtf_py(const py::bytes& data, const mtef::RtfOptions& opt) {
    std::string buf = data;
    std::string rtf = mtef::mtef_to_rtf(
        reinterpret_cast<const uint8_t*>(buf.data()), buf.size(), opt);
    if (rtf.empty()) throw std::runtime_error("mtef_to_rtf: parse failed");
    return rtf;
}

std::string tex_to_mathml_py(const std::string& latex,
                             const mtef::MathMLOptions& opt) {
    std::string m = mtef::tex_to_mathml(latex, opt);
    if (m.empty()) throw std::runtime_error("tex_to_mathml: emit failed for: " + latex);
    return m;
}

py::bytes tex_to_gvml_py(const std::string& latex, double size_pt,
                         bool display) {
    std::string pkg = mtef::tex_to_gvml(latex, size_pt, display);
    if (pkg.empty()) throw std::runtime_error("tex_to_gvml: emit failed for: " + latex);
    return py::bytes(pkg);
}

std::string mtef_to_mathml_py(const py::bytes& data, const mtef::MathMLOptions& opt) {
    std::string buf = data;
    std::string m = mtef::mtef_to_mathml(
        reinterpret_cast<const uint8_t*>(buf.data()), buf.size(), opt);
    if (m.empty()) throw std::runtime_error("mtef_to_mathml: parse failed");
    return m;
}

py::bytes tex_to_emf_py(const std::string& latex, const mtef::SvgStyle& style) {
    std::string emf = mtef::tex_to_emf(latex, style);
    if (emf.empty()) throw std::runtime_error("tex_to_emf: render failed for: " + latex);
    return py::bytes(emf);
}

py::bytes tex_to_png_py(const std::string& latex, const mtef::SvgStyle& style,
                        double scale) {
    std::string png = mtef::tex_to_png(latex, style, scale);
    if (png.empty()) throw std::runtime_error("tex_to_png: render failed for: " + latex);
    return py::bytes(png);
}

std::string tex_dump_tree_py(const std::string& latex) {
    std::unique_ptr<mtef::LineNode> root = mtef::parse_latex(latex);
    if (!root) throw std::runtime_error("parse_latex failed for: " + latex);
    return mtef::dump_latex_tree(*root);
}

std::string dump_tree_py(const py::bytes& data, bool run_passes) {
    std::string buf = data;
    return mtef::dump_tree(reinterpret_cast<const uint8_t*>(buf.data()),
                           buf.size(), run_passes);
}

/* A path arrives from Python as UTF-8.  Handing that to a narrow ifstream on
 * Windows reads it in the ANSI code page instead, so every file whose name is
 * not ASCII was unreachable -- which in this lab is most of them.  Widen it
 * and use the wide overload MSVC provides. */
#ifdef _WIN32
std::wstring widen(const std::string& utf8) {
    if (utf8.empty()) return std::wstring();
    const int n = MultiByteToWideChar(CP_UTF8, 0, utf8.data(),
                                      int(utf8.size()), nullptr, 0);
    std::wstring out(size_t(n > 0 ? n : 0), L'\0');
    if (n > 0)
        MultiByteToWideChar(CP_UTF8, 0, utf8.data(), int(utf8.size()),
                            &out[0], n);
    return out;
}
#define EQ_PATH(p) widen(p)
#else
#define EQ_PATH(p) (p)
#endif

py::bytes read_eqn_py(const std::string& path) {
    std::ifstream f(EQ_PATH(path), std::ios::binary);
    if (!f) throw std::runtime_error("read_eqn: cannot open " + path);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    return py::bytes(data);
}

void write_eqn_py(const std::string& path, const py::bytes& data) {
    std::string buf = data;
    std::ofstream f(EQ_PATH(path), std::ios::binary);
    if (!f) throw std::runtime_error("write_eqn: cannot write " + path);
    f.write(buf.data(), std::streamsize(buf.size()));
}

}  // namespace

PYBIND11_MODULE(_equation, m) {
    m.doc() = "MTEF (Equation Editor 3.x / MathType) codec and SVG renderer.\n"
              "Equation Editor itself is only a GUI for humans; every "
              "programmatic path runs here.";

    py::class_<mtef::SvgStyle>(m, "SvgStyle",
        "Type sizes in points, mirroring the editor's Sizes dialog.")
        .def(py::init<>())
        .def_readwrite("full", &mtef::SvgStyle::full)
        .def_readwrite("sub", &mtef::SvgStyle::sub)
        .def_readwrite("sub2", &mtef::SvgStyle::sub2)
        .def_readwrite("sym", &mtef::SvgStyle::sym)
        .def_readwrite("subsym", &mtef::SvgStyle::subsym)
        .def_readwrite("serif", &mtef::SvgStyle::serif)
        .def_readwrite("symbol", &mtef::SvgStyle::symbol)
        .def_readwrite("padding", &mtef::SvgStyle::padding)
        .def_readwrite("empty_slot_em", &mtef::SvgStyle::empty_slot_em,
                       "Room an EMPTY slot takes, as a fraction of the type "
                       "size.  0 -- the default -- is what TeX does, and is "
                       "what a picture wants.  The editor sets 0.55 so a "
                       "template nobody has typed into is visible.")
        .def("__repr__", [](const mtef::SvgStyle& s) {
            return "<SvgStyle full=" + std::to_string(s.full) +
                   " sub=" + std::to_string(s.sub) +
                   " sym=" + std::to_string(s.sym) + ">";
        });

    m.def("tex_to_mtef", &tex_to_mtef_py, py::arg("latex"),
          "LaTeX -> MTEF binary.");
    m.def("mtef_to_tex", &mtef_to_tex_py, py::arg("data"),
          "MTEF binary -> LaTeX.");
    m.def("mtef_to_latex", &mtef_to_latex_py, py::arg("data"),
          "MTEF binary -> LaTeX, through the EDITOR's parser and emitter -- "
          "the reading that will actually be drawn and edited.  mtef_to_tex "
          "answers with the legacy standalone converter instead.");
    m.def("mtef_to_svg", &mtef_to_svg_py, py::arg("data"),
          py::arg("style") = mtef::SvgStyle(), "MTEF binary -> SVG.");
    m.def("tex_to_svg", &tex_to_svg_py, py::arg("latex"),
          py::arg("style") = mtef::SvgStyle(), "LaTeX -> SVG.");
    py::class_<mtef::OmmlOptions>(m, "OmmlOptions",
        "Office Math (OMML) emission options.")
        .def(py::init<>())
        .def_readwrite("display", &mtef::OmmlOptions::display)
        .def_readwrite("declare_namespace", &mtef::OmmlOptions::declare_namespace)
        .def_readwrite("italic_variables", &mtef::OmmlOptions::italic_variables);

    m.def("mtef_to_omml", &mtef_to_omml_py, py::arg("data"),
          py::arg("options") = mtef::OmmlOptions(),
          "MTEF binary -> OMML (Office-native, editable equation).");
    m.def("tex_to_omml", &tex_to_omml_py, py::arg("latex"),
          py::arg("options") = mtef::OmmlOptions(),
          "LaTeX -> OMML (Office-native, editable equation).");

    py::class_<mtef::RtfOptions>(m, "RtfOptions",
        "How an equation is written into RTF, the clipboard form Office takes.")
        .def(py::init<>())
        .def_readwrite("display", &mtef::RtfOptions::display)
        .def_readwrite("font_size_pt", &mtef::RtfOptions::font_size_pt);

    m.def("tex_to_rtf", &tex_to_rtf_py, py::arg("latex"),
          py::arg("options") = mtef::RtfOptions(),
          "LaTeX -> a complete RTF document carrying one native equation.");
    m.def("mtef_to_rtf", &mtef_to_rtf_py, py::arg("data"),
          py::arg("options") = mtef::RtfOptions(),
          "MTEF binary -> a complete RTF document carrying one equation.");

    py::class_<mtef::MathMLOptions>(m, "MathMLOptions",
        "How an equation is written as MathML -- the one clipboard format the "
        "whole of Office reads as maths.")
        .def(py::init<>())
        .def_readwrite("display", &mtef::MathMLOptions::display)
        .def_readwrite("declare_namespace", &mtef::MathMLOptions::declare_namespace);

    m.def("tex_to_mathml", &tex_to_mathml_py, py::arg("latex"),
          py::arg("options") = mtef::MathMLOptions(),
          "LaTeX -> MathML (Word, PowerPoint and Excel all read it as maths).");

    m.attr("PASTE_SIZE_PT") = mtef::kPasteSizePt;
    m.def("tex_to_gvml", &tex_to_gvml_py, py::arg("latex"),
          py::arg("size_pt") = mtef::kPasteSizePt,
          py::arg("display") = false,
          "LaTeX -> an Art::GVML ClipFormat package: PowerPoint's own shape "
          "format, and the only clipboard format that can state the SIZE the "
          "equation should arrive at.  MathML cannot -- mathsize is ignored "
          "and the destination box wins, which is why a pasted equation used "
          "to come out at 18 pt.  Returns the OPC package bytes.");
    m.def("mtef_to_mathml", &mtef_to_mathml_py, py::arg("data"),
          py::arg("options") = mtef::MathMLOptions(),
          "MTEF binary -> MathML.");

    m.def("tex_to_emf", &tex_to_emf_py, py::arg("latex"),
          py::arg("style") = mtef::SvgStyle(),
          "LaTeX -> an enhanced metafile, the vector picture Office and "
          "Google Drawings both take.");
    m.def("tex_to_png", &tex_to_png_py, py::arg("latex"),
          py::arg("style") = mtef::SvgStyle(), py::arg("scale") = 4.0,
          "LaTeX -> a PNG, for the targets that accept nothing else.");

    m.def("tex_normalize", &tex_normalize_py, py::arg("latex"),
          "LaTeX -> tree -> LaTeX (the shape an edited equation is saved in).");

    m.def("tex_dump_tree", &tex_dump_tree_py, py::arg("latex"),
          "Indented text dump of the tree the LaTeX parser builds.");

    m.def("dump_tree", &dump_tree_py, py::arg("data"),
          py::arg("run_passes") = true,
          "Indented text dump of the parsed node tree (diagnostic).");

    py::class_<mtef::Equation::SelectionBox>(m, "SelectionBox",
        "The rectangle an editor highlights.  `found` is false when nothing "
        "is selected.")
        .def_readonly("found", &mtef::Equation::SelectionBox::found)
        .def_readonly("x0", &mtef::Equation::SelectionBox::x0)
        .def_readonly("x1", &mtef::Equation::SelectionBox::x1)
        .def_readonly("top", &mtef::Equation::SelectionBox::top)
        .def_readonly("bottom", &mtef::Equation::SelectionBox::bottom);

    py::class_<mtef::Equation>(m, "Equation",
        "The editing model: an insertion point that lives inside the structure, "
        "templates you tab through, backspace that unwraps.  LaTeX in, LaTeX out.")
        .def(py::init<>())
        .def("load_latex", &mtef::Equation::load_latex, py::arg("latex"))
        .def("latex", &mtef::Equation::latex)
        .def("omml", &mtef::Equation::omml,
             py::arg("options") = mtef::OmmlOptions())
        .def("svg", &mtef::Equation::svg, py::arg("style") = mtef::SvgStyle())
        .def("caret", &mtef::Equation::caret)
        .def("caret_geometry", [](const mtef::Equation& e, const mtef::SvgStyle& st) {
                 auto g = e.caret_geometry(st);
                 return py::make_tuple(g.found, g.x, g.top, g.bottom);
             }, py::arg("style") = mtef::SvgStyle(),
             "Where to draw the caret: (found, x, top, bottom), in points "
             "relative to the equation's origin.")
        .def("move_to_point", &mtef::Equation::move_to_point,
             py::arg("x"), py::arg("y"), py::arg("style") = mtef::SvgStyle(),
             "Put the caret where a click landed.")
        .def("extents", [](const mtef::Equation& e, const mtef::SvgStyle& st) {
                 double w = 0, a = 0, d = 0;
                 e.extents(w, a, d, st);
                 return py::make_tuple(w, a, d);
             }, py::arg("style") = mtef::SvgStyle(),
             "The equation's box: (width, ascent, descent) in points.")
        .def("move_left", &mtef::Equation::move_left)
        .def("move_right", &mtef::Equation::move_right)
        .def("next_slot", &mtef::Equation::next_slot)
        .def("prev_slot", &mtef::Equation::prev_slot)
        .def("move_out", &mtef::Equation::move_out)
        .def("move_home", &mtef::Equation::move_home)
        .def("move_end", &mtef::Equation::move_end)
        .def("insert_text", &mtef::Equation::insert_text, py::arg("text"))
        .def("insert_symbol", &mtef::Equation::insert_symbol, py::arg("command"))
        .def("insert_template", &mtef::Equation::insert_template, py::arg("kind"))
        .def("insert_latex", &mtef::Equation::insert_latex, py::arg("latex"),
             "Paste: parse LaTeX and put it in at the caret, replacing any "
             "selection, in one undo step.")
        .def("backspace", &mtef::Equation::backspace)
        .def("erase", &mtef::Equation::erase)
        /* Selection: a range within one slot.  Equation Editor's Select All,
         * Cut, Copy and Clear all rest on these. */
        .def("set_style", &mtef::Equation::set_style, py::arg("name"),
             "Equation Editor's Style menu, applied to the selection: math, "
             "text, function, variable, vector (bold), greek.")
        .def_static("styles", &mtef::Equation::styles)
        .def("has_selection", &mtef::Equation::has_selection)
        .def("clear_selection", &mtef::Equation::clear_selection)
        .def("select_all", &mtef::Equation::select_all)
        .def("extend_left", &mtef::Equation::extend_left)
        .def("extend_right", &mtef::Equation::extend_right)
        .def("extend_home", &mtef::Equation::extend_home)
        .def("extend_end", &mtef::Equation::extend_end)
        .def("extend_to_point", &mtef::Equation::extend_to_point,
             py::arg("x"), py::arg("y"), py::arg("style") = mtef::SvgStyle(),
             "Drag the selection to a point, keeping the anchor.")
        .def("selected_latex", &mtef::Equation::selected_latex,
             "The selected nodes as LaTeX, or an empty string.")
        .def("delete_selection", &mtef::Equation::delete_selection)
        .def("selection_geometry", &mtef::Equation::selection_geometry,
             py::arg("style") = mtef::SvgStyle(),
             "The box an editor highlights.")
        .def("undo", &mtef::Equation::undo)
        .def("redo", &mtef::Equation::redo)
        .def("command", &mtef::Equation::command, py::arg("name"))
        .def_static("templates", &mtef::Equation::templates)
        .def_static("shortcuts", []() {
            py::list out;
            for (const auto& b : mtef::Equation::shortcuts())
                out.append(py::make_tuple(b.chord, b.command, b.label));
            return out;
        }, "Equation Editor 3.0 key chords as (chord, command, label).")
        .def("style", &mtef::Equation::style,
             "The style typing goes in.  Equation Editor's Style menu is a "
             "mode, not an operation on what is highlighted.")
        .def("press", [](mtef::Equation& e, unsigned vk, bool ctrl,
                         bool shift, bool alt) {
            /* The half-entered chord is carried by the equation itself.  A
             * static table keyed on its address gave one equation's pending
             * prefix to whatever was allocated at that address next, which
             * showed up as a test that passed alone and failed in a run. */
            mtef::KeyState st;
            mtef::Equation::PendingChord& p = e.pending_chord();
            for (size_t i = 0; i < p.vk.size(); ++i) {
                mtef::Step s2;
                s2.vk = unsigned(p.vk[i]);
                s2.ctrl  = (p.mods[i] & 1) != 0;
                s2.shift = (p.mods[i] & 2) != 0;
                s2.alt   = (p.mods[i] & 4) != 0;
                st.pending.push_back(s2);
            }
            const mtef::KeyResult r =
                mtef::press_key(e, st, vk, ctrl, shift, alt);
            p.vk.clear();
            p.mods.clear();
            for (const mtef::Step& s2 : st.pending) {
                p.vk.push_back(int(s2.vk));
                p.mods.push_back((s2.ctrl ? 1 : 0) | (s2.shift ? 2 : 0) |
                                 (s2.alt ? 4 : 0));
            }
            return r == mtef::KeyResult::Ignored ? "ignored"
                 : r == mtef::KeyResult::Pending ? "pending" : "consumed";
        }, py::arg("vk"), py::arg("ctrl") = false, py::arg("shift") = false,
           py::arg("alt") = false,
           "Press a key, the way the window does -- but with the modifiers "
           "passed rather than read from the thread's input state, so this "
           "needs no window and no keyboard.")
        .def_static("chord_steps", [](const std::string& chord) {
            mtef::Chord c;
            py::list out;
            if (!mtef::parse_chord(chord, c)) return out;
            for (const mtef::Step& s : c.steps)
                out.append(py::make_tuple(s.vk, s.ctrl, s.shift, s.alt));
            return out;
        }, py::arg("chord"),
           "What the window will actually watch for, as (vk, ctrl, shift, alt) "
           "per key press.  Empty if the chord names no key this knows.")
        .def("__repr__", [](const mtef::Equation& e) {
            return "<Equation " + e.latex() + " caret=" + e.caret() + ">";
        });

    /* The typesetting parameters the math font itself carries.  The layout
     * used to guess these; exposing them is what lets a test say the guess is
     * gone and the font is being read. */
    m.def("math_constants", []() {
        const mtef::MathFont& f = mtef::MathFont::math();
        if (!f.ok()) throw std::runtime_error(f.why_not());
        const mtef::MathConstants& c = f.constants();
        py::dict d;
        d["script_percent"] = c.scriptPercentScaleDown;
        d["script_script_percent"] = c.scriptScriptPercentScaleDown;
        d["axis_height"] = c.axisHeight;
        d["fraction_rule_thickness"] = c.fractionRuleThickness;
        d["fraction_numerator_gap_min"] = c.fractionNumeratorGapMin;
        d["fraction_num_display_gap_min"] = c.fractionNumDisplayStyleGapMin;
        d["fraction_num_display_shift_up"] = c.fractionNumeratorDisplayStyleShiftUp;
        d["fraction_denom_display_shift_down"] = c.fractionDenominatorDisplayStyleShiftDown;
        d["radical_vertical_gap"] = c.radicalVerticalGap;
        d["radical_display_vertical_gap"] = c.radicalDisplayStyleVerticalGap;
        d["radical_rule_thickness"] = c.radicalRuleThickness;
        d["radical_extra_ascender"] = c.radicalExtraAscender;
        d["radical_kern_before_degree"] = c.radicalKernBeforeDegree;
        d["radical_kern_after_degree"] = c.radicalKernAfterDegree;
        d["radical_degree_bottom_raise"] = c.radicalDegreeBottomRaisePercent;
        d["overbar_rule_thickness"] = c.overbarRuleThickness;
        d["overbar_vertical_gap"] = c.overbarVerticalGap;
        d["accent_base_height"] = c.accentBaseHeight;
        d["superscript_shift_up"] = c.superscriptShiftUp;
        d["subscript_shift_down"] = c.subscriptShiftDown;
        d["display_operator_min_height"] = c.displayOperatorMinHeight;
        d["math_leading"] = c.mathLeading;
        d["upem"] = 1.0;   /* everything above is already in em */
        return d;
    }, "Cambria Math's MATH constants, in em.");

    m.def("math_glyph", [](uint32_t codepoint) {
        const mtef::MathFont& f = mtef::MathFont::math();
        if (!f.ok()) throw std::runtime_error(f.why_not());
        return f.glyph_for(codepoint);
    }, py::arg("codepoint"), "Glyph id, or 0 when the font has no such glyph.");

    m.def("math_stretch", [](uint32_t codepoint) {
        const mtef::MathFont& f = mtef::MathFont::math();
        if (!f.ok()) throw std::runtime_error(f.why_not());
        const uint16_t g = f.glyph_for(codepoint);
        py::dict d;
        d["glyph"] = g;
        py::list variants, assembly;
        if (const mtef::Stretch* s = f.vertical(g)) {
            for (const auto& v : s->variants)
                variants.append(py::make_tuple(v.first, v.second));
            for (const mtef::GlyphPart& p : s->assembly)
                assembly.append(py::make_tuple(p.glyph, p.fullAdvance,
                                               p.extender));
        }
        d["variants"] = variants;
        d["assembly"] = assembly;
        d["min_overlap"] = f.min_connector_overlap();
        return d;
    }, py::arg("codepoint"),
       "How a character grows: the font's ready-made taller drawings, and the "
       "parts to assemble one taller than any of them.");

    m.def("math_variant_for_height", [](uint32_t codepoint, double em) {
        const mtef::MathFont& f = mtef::MathFont::math();
        if (!f.ok()) throw std::runtime_error(f.why_not());
        double got = 0;
        const uint16_t g = f.vertical_variant(f.glyph_for(codepoint), em, &got);
        return py::make_tuple(g, got);
    }, py::arg("codepoint"), py::arg("em"),
       "The smallest drawing at least that tall, and how tall it is.");

    m.def("tex_metrics", [](const std::string& latex,
                            const mtef::SvgStyle& style) {
        double w = 0, asc = 0, desc = 0;
        if (!mtef::tex_box(latex, style, w, asc, desc))
            throw std::runtime_error("tex_metrics: render failed for: " + latex);
        return py::make_tuple(w, asc, desc);
    }, py::arg("latex"), py::arg("style"),
       "Width, height above the baseline and depth below it, in points -- the "
       "same three numbers TeX reports for a box, so the two can be compared "
       "without measuring a picture.");

    py::enum_<mtef::AtomKind>(m, "AtomKind", "TeX's atom classes.")
        .value("Ord", mtef::AtomOrd).value("Op", mtef::AtomOp)
        .value("Bin", mtef::AtomBin).value("Rel", mtef::AtomRel)
        .value("Open", mtef::AtomOpen).value("Close", mtef::AtomClose)
        .value("Punct", mtef::AtomPunct).value("Inner", mtef::AtomInner);

    m.def("atom_kind", [](uint32_t cp) { return mtef::atom_kind(cp); },
          py::arg("codepoint"),
          "Which class a character belongs to.  What decides whether \"a*b\" "
          "reads as a product or as three letters.");
    m.def("atom_space_mu", &mtef::atom_space_mu, py::arg("left"),
          py::arg("right"),
          "Space between two atoms, in eighteenths of an em.");

    py::class_<mtef::MdSegment> seg(m, "MdSegment",
        "One span of a Markdown file: prose, code, or an equation.");
    py::enum_<mtef::MdSegment::Kind>(seg, "Kind")
        .value("Text", mtef::MdSegment::kText)
        .value("InlineMath", mtef::MdSegment::kInlineMath)
        .value("DisplayMath", mtef::MdSegment::kDisplayMath)
        .value("CodeSpan", mtef::MdSegment::kCodeSpan)
        .value("CodeBlock", mtef::MdSegment::kCodeBlock)
        .export_values();
    seg.def_readonly("kind", &mtef::MdSegment::kind)
       .def_readonly("open", &mtef::MdSegment::open)
       .def_readonly("body", &mtef::MdSegment::body)
       .def_readonly("close", &mtef::MdSegment::close)
       .def_property_readonly("is_math", &mtef::MdSegment::is_math)
       .def_property_readonly("source", &mtef::MdSegment::source)
       .def("__repr__", [](const mtef::MdSegment& s) {
           return "<MdSegment " + std::to_string(int(s.kind)) + " " +
                  s.source().substr(0, 40) + ">";
       });

    py::class_<mtef::MarkdownDoc>(m, "MarkdownDoc",
        "A Markdown file seen as text with equations in it.  Loading and "
        "saving an untouched file reproduces it byte for byte.")
        .def(py::init<>())
        .def("load", &mtef::MarkdownDoc::load, py::arg("markdown"))
        .def("text", &mtef::MarkdownDoc::text)
        .def("segments", &mtef::MarkdownDoc::segments)
        .def("math_count", &mtef::MarkdownDoc::math_count)
        .def("math_latex", &mtef::MarkdownDoc::math_latex, py::arg("index"))
        .def("math_is_display", &mtef::MarkdownDoc::math_is_display, py::arg("index"))
        .def("set_math_latex", &mtef::MarkdownDoc::set_math_latex,
             py::arg("index"), py::arg("latex"))
        .def("math_segment_index", &mtef::MarkdownDoc::math_segment_index,
             py::arg("index"));

    py::class_<mtef::MdBlock> blk(m, "MdBlock",
        "One block of a Markdown file: a heading, a paragraph, a list item, "
        "fenced code, or a run of blank lines.");
    py::enum_<mtef::MdBlock::Kind>(blk, "Kind")
        .value("Paragraph", mtef::MdBlock::kParagraph)
        .value("Heading", mtef::MdBlock::kHeading)
        .value("Bullet", mtef::MdBlock::kBullet)
        .value("Numbered", mtef::MdBlock::kNumbered)
        .value("Code", mtef::MdBlock::kCode)
        .value("Blank", mtef::MdBlock::kBlank)
        .export_values();
    blk.def_readonly("kind", &mtef::MdBlock::kind)
       .def_readonly("level", &mtef::MdBlock::level)
       .def_readonly("text", &mtef::MdBlock::text)
       .def_readonly("info", &mtef::MdBlock::info)
       .def_readonly("source", &mtef::MdBlock::source)
       .def("__repr__", [](const mtef::MdBlock& b) {
           return "<MdBlock " + std::to_string(int(b.kind)) + " " +
                  b.text.substr(0, 40) + ">";
       });

    py::class_<mtef::Equation::PaletteItem>(m, "PaletteItem",
        "One cell of a palette: a symbol command, or a template name.")
        .def_readonly("command", &mtef::Equation::PaletteItem::command)
        .def_readonly("code", &mtef::Equation::PaletteItem::code)
        .def_readonly("is_template", &mtef::Equation::PaletteItem::is_template)
        .def("__repr__", [](const mtef::Equation::PaletteItem& i) {
            return "<PaletteItem " + i.command + ">";
        });

    py::class_<mtef::Equation::PaletteGroup>(m, "PaletteGroup",
        "One button on the palette bar and what its popup holds.")
        .def_readonly("name", &mtef::Equation::PaletteGroup::name)
        .def_readonly("items", &mtef::Equation::PaletteGroup::items)
        /* The member the bar button wears.  Exposed so a test can see what
         * the WINDOW sees: membership alone proves nothing, because a
         * mis-spelled name silently falls back to the first item and the
         * button quietly stays wrong -- which is exactly what happened when
         * the table reached the compiler with one backslash too few. */
        .def_property_readonly("icon", &mtef::Equation::PaletteGroup::icon,
             py::return_value_policy::reference_internal,
             "The member this group's palette button draws.")
        .def("__repr__", [](const mtef::Equation::PaletteGroup& g) {
            return "<PaletteGroup " + g.name + " n=" +
                   std::to_string(g.items.size()) + ">";
        });

    m.def("symbol_palettes", &mtef::Equation::symbol_palettes,
          py::return_value_policy::reference,
          "The symbol row of the palette bar, built from the shared command "
          "table so it cannot offer what the editor will not insert.");
    m.def("template_palettes", &mtef::Equation::template_palettes,
          py::return_value_policy::reference,
          "The template row of the palette bar.");

    m.def("tex_empty_slots", [](const std::string& latex,
                                const mtef::SvgStyle& st) {
              auto root = mtef::parse_latex(latex);
              mtef::Layout L = mtef::layout_math(*root, st);
              py::list out;
              for (const mtef::SlotBox& b : L.empty_slots)
                  out.append(py::make_tuple(b.x, b.y, b.w, b.h));
              return out;
          }, py::arg("latex"), py::arg("style") = mtef::SvgStyle(),
             "Where the empty slots of a half-typed equation are, as "
             "(x, y, w, h) in points.  The editor draws these as dotted boxes "
             "so the structure is visible; a picture must not.");

    m.def("tex_to_dib", [](const std::string& latex, const mtef::SvgStyle& st,
                           double scale) {
              return py::bytes(mtef::tex_to_dib(latex, st, scale));
          }, py::arg("latex"), py::arg("style") = mtef::SvgStyle(),
             py::arg("scale") = 4.0,
             "The equation as a packed DIB, for CF_DIB -- the image format an "
             "application reads when it pastes a picture.");

    m.def("md_blocks", &mtef::md_blocks, py::arg("markdown"),
          "Split Markdown into blocks.  Concatenating every source rebuilds "
          "the file exactly.");

    /* ---- the document, laid out ----------------------------------------- */

    py::class_<mtef::DocStyle>(m, "DocStyle",
        "Point sizes, fonts and spacing a Markdown document is set in.")
        .def(py::init<>())
        .def_readwrite("body", &mtef::DocStyle::body)
        .def_readwrite("mono", &mtef::DocStyle::mono)
        .def_readwrite("line_spacing", &mtef::DocStyle::line_spacing)
        .def_readwrite("para_gap", &mtef::DocStyle::para_gap)
        .def_readwrite("list_indent", &mtef::DocStyle::list_indent)
        .def_readwrite("margin", &mtef::DocStyle::margin)
        .def_readwrite("text_font", &mtef::DocStyle::text_font)
        .def_readwrite("mono_font", &mtef::DocStyle::mono_font)
        .def_readwrite("math_scale", &mtef::DocStyle::math_scale)
        .def_property("heading",
            [](const mtef::DocStyle& s) {
                return std::vector<double>(s.heading, s.heading + 6);
            },
            [](mtef::DocStyle& s, const std::vector<double>& v) {
                for (size_t i = 0; i < 6 && i < v.size(); ++i) s.heading[i] = v[i];
            },
            "Point size of each heading level, 1 to 6.");

    py::class_<mtef::DocRun>(m, "DocRun", "A run of text at a position.")
        .def_readonly("text", &mtef::DocRun::text)
        .def_readonly("x", &mtef::DocRun::x)
        .def_readonly("baseline", &mtef::DocRun::baseline)
        .def_readonly("size", &mtef::DocRun::size)
        .def_readonly("bold", &mtef::DocRun::bold)
        .def_readonly("italic", &mtef::DocRun::italic)
        .def_readonly("mono", &mtef::DocRun::mono)
        .def("__repr__", [](const mtef::DocRun& r) {
            return "<DocRun " + r.text.substr(0, 30) + ">";
        });

    py::class_<mtef::DocMath>(m, "DocMath",
        "An equation placed in the document, carrying its own math layout so "
        "the same layout serves the screen and the Office paste.")
        .def_readonly("x", &mtef::DocMath::x)
        .def_readonly("baseline", &mtef::DocMath::baseline)
        .def_readonly("latex", &mtef::DocMath::latex)
        .def_readonly("block", &mtef::DocMath::block)
        .def_readonly("index", &mtef::DocMath::index)
        .def_readonly("display", &mtef::DocMath::display)
        .def_property_readonly("width",
            [](const mtef::DocMath& d) { return d.layout.w; })
        .def_property_readonly("ascent",
            [](const mtef::DocMath& d) { return d.layout.asc; })
        .def_property_readonly("descent",
            [](const mtef::DocMath& d) { return d.layout.desc; })
        .def("__repr__", [](const mtef::DocMath& d) {
            return "<DocMath " + d.latex.substr(0, 30) + ">";
        });

    py::class_<mtef::DocBlockBox>(m, "DocBlockBox",
        "Where a block ended up, so a click can select it.")
        .def_readonly("block", &mtef::DocBlockBox::block)
        .def_readonly("kind", &mtef::DocBlockBox::kind)
        .def_readonly("top", &mtef::DocBlockBox::top)
        .def_readonly("bottom", &mtef::DocBlockBox::bottom);

    py::class_<mtef::DocLayout>(m, "DocLayout",
        "A Markdown document laid out to a width, in points.")
        .def_readonly("width", &mtef::DocLayout::width)
        .def_readonly("height", &mtef::DocLayout::height)
        .def_readonly("runs", &mtef::DocLayout::runs)
        .def_readonly("maths", &mtef::DocLayout::maths)
        .def_readonly("blocks", &mtef::DocLayout::blocks)
        .def("block_at", [](const mtef::DocLayout& d, double x, double y) {
                 return mtef::block_at(d, x, y);
             }, py::arg("x"), py::arg("y"),
             "Which block a point falls in, or -1.")
        .def("math_at", [](const mtef::DocLayout& d, double x, double y) {
                 return mtef::math_at(d, x, y);
             }, py::arg("x"), py::arg("y"),
             "Which equation a point falls in, or -1 -- the test for opening "
             "the equation widget rather than the text editor.");

    m.def("layout_markdown", &mtef::layout_markdown,
          py::arg("markdown"), py::arg("width"),
          py::arg("style") = mtef::DocStyle(),
          "Lay a Markdown document out to a width in points.");

    m.def("read_eqn", &read_eqn_py, py::arg("path"),
          "Read a .eqn file (raw MTEF, no OLE header).");
    m.def("write_eqn", &write_eqn_py, py::arg("path"), py::arg("data"),
          "Write a .eqn file the editor GUI can open.");
}
