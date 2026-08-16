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
#include <stdexcept>
#include <string>
#include <vector>

#include "mtef2tex.h"
#include "tex2mtef.h"
#include "mtef_svg.h"
#include "mtef_omml.h"
#include "mtef_dump.h"
#include "tex_parser.h"
#include "latex_emitter.h"
#include "eq_edit.h"
#include "md_doc.h"

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

py::bytes read_eqn_py(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("read_eqn: cannot open " + path);
    std::string data((std::istreambuf_iterator<char>(f)),
                     std::istreambuf_iterator<char>());
    return py::bytes(data);
}

void write_eqn_py(const std::string& path, const py::bytes& data) {
    std::string buf = data;
    std::ofstream f(path, std::ios::binary);
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
        .def("__repr__", [](const mtef::SvgStyle& s) {
            return "<SvgStyle full=" + std::to_string(s.full) +
                   " sub=" + std::to_string(s.sub) +
                   " sym=" + std::to_string(s.sym) + ">";
        });

    m.def("tex_to_mtef", &tex_to_mtef_py, py::arg("latex"),
          "LaTeX -> MTEF binary.");
    m.def("mtef_to_tex", &mtef_to_tex_py, py::arg("data"),
          "MTEF binary -> LaTeX.");
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

    m.def("tex_normalize", &tex_normalize_py, py::arg("latex"),
          "LaTeX -> tree -> LaTeX (the shape an edited equation is saved in).");

    m.def("tex_dump_tree", &tex_dump_tree_py, py::arg("latex"),
          "Indented text dump of the tree the LaTeX parser builds.");

    m.def("dump_tree", &dump_tree_py, py::arg("data"),
          py::arg("run_passes") = true,
          "Indented text dump of the parsed node tree (diagnostic).");

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
        .def("backspace", &mtef::Equation::backspace)
        .def("erase", &mtef::Equation::erase)
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
        .def("__repr__", [](const mtef::Equation& e) {
            return "<Equation " + e.latex() + " caret=" + e.caret() + ">";
        });

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

    m.def("read_eqn", &read_eqn_py, py::arg("path"),
          "Read a .eqn file (raw MTEF, no OLE header).");
    m.def("write_eqn", &write_eqn_py, py::arg("path"), py::arg("data"),
          "Write a .eqn file the editor GUI can open.");
}
