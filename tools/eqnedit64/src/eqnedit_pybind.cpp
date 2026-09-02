/* Headless test/programming surface for Eqnedit64's TeX editing core. */
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <string>

#include "equation_edit.h"
#include "equation_render.h"
#include "latex_emitter.h"
#include "math_symbols.h"
#include "mathml_emitter.h"
#include "palettes.h"
#include "tex_document.h"
#include "tex_parser.h"

namespace py = pybind11;

namespace {

std::string tex_to_svg_checked(const std::string& tex,
                               const eqnedit::SvgStyle& style) {
    std::string svg = eqnedit::tex_to_svg(tex, style);
    if (svg.empty()) throw std::runtime_error("TeX rendering failed");
    return svg;
}

std::string normalize_tree(const std::string& tex) {
    std::unique_ptr<eqnedit::LineNode> root = eqnedit::parse_latex(tex);
    if (!root) throw std::runtime_error("TeX parsing failed");
    return eqnedit::tree_to_latex(*root);
}

}  // namespace

PYBIND11_MODULE(eqnedit_core, module) {
    module.doc() = "Eqnedit64 TeX parser, structural editor, and renderer";
    module.attr("MAX_NESTING_DEPTH") = eqnedit::kMaxNestingDepth;

    py::class_<eqnedit::SvgStyle>(module, "SvgStyle")
        .def(py::init<>())
        .def_readwrite("full", &eqnedit::SvgStyle::full)
        .def_readwrite("sub", &eqnedit::SvgStyle::sub)
        .def_readwrite("sub2", &eqnedit::SvgStyle::sub2)
        .def_readwrite("sym", &eqnedit::SvgStyle::sym)
        .def_readwrite("subsym", &eqnedit::SvgStyle::subsym)
        .def_readwrite("serif", &eqnedit::SvgStyle::serif)
        .def_readwrite("symbol", &eqnedit::SvgStyle::symbol)
        .def_readwrite("cjk", &eqnedit::SvgStyle::cjk)
        .def_readwrite("padding", &eqnedit::SvgStyle::padding);

    module.def("tex_to_svg", &tex_to_svg_checked, py::arg("tex"),
               py::arg("style") = eqnedit::SvgStyle());
    module.def("tex_normalize", &normalize_tree, py::arg("tex"));
    module.def("tex_to_mathml", &eqnedit::latex_to_mathml, py::arg("tex"),
               py::arg("point_size") = 24.0);
    module.def("tex_to_office_mathml_fragment",
               &eqnedit::latex_to_office_mathml_fragment, py::arg("tex"),
               py::arg("point_size") = 18.0);
    module.def("normalize_paste", &eqnedit::normalize_tex_paste, py::arg("text"));
    /* The whole symbol table, so a test can sweep every command instead of
     * the handful a hand-written corpus happens to mention. */
    module.def("math_font_loaded", &eqnedit::math_font_loaded);
    module.def("symbol_commands", &eqnedit::latex_symbol_commands);
    /* The toolbar catalogue, so a test can prove every symbol and every
     * template has a place a mouse can reach. */
    module.def("palettes", []() {
        py::list result;
        for (const auto& palette : eqnedit::palettes()) {
            py::list items;
            for (const auto& item : palette.items)
                items.append(py::make_tuple(item.command, item.face,
                                            item.label));
            result.append(py::make_tuple(palette.title, palette.face,
                                         palette.columns, items));
        }
        return result;
    });
    module.def("symbol_palette_count", &eqnedit::symbol_palette_count);
    module.def("palette_categories", []() {
        py::list result;
        for (const auto& category : eqnedit::palette_categories())
            result.append(py::make_tuple(category.title,
                                         category.paletteIndices));
        return result;
    });
    module.def("compose_document",
               [](const std::string& body, bool numbered) {
                   return eqnedit::compose_tex_document(body, numbered);
               },
               py::arg("body"), py::arg("numbered") = true);

    py::class_<eqnedit::Equation>(module, "Equation")
        .def(py::init<>())
        .def("load_latex", &eqnedit::Equation::load_latex, py::arg("tex"))
        .def("replace_latex", &eqnedit::Equation::replace_latex,
             py::arg("tex"), py::arg("checkpoint_first") = true)
        .def("last_error", &eqnedit::Equation::last_error,
             py::return_value_policy::reference_internal)
        .def("latex", &eqnedit::Equation::latex)
        .def("svg", &eqnedit::Equation::svg,
             py::arg("style") = eqnedit::SvgStyle())
        .def("caret", &eqnedit::Equation::caret)
        .def("move_left", &eqnedit::Equation::move_left)
        .def("move_right", &eqnedit::Equation::move_right)
        .def("next_slot", &eqnedit::Equation::next_slot)
        .def("prev_slot", &eqnedit::Equation::prev_slot)
        .def("move_out", &eqnedit::Equation::move_out)
        .def("move_home", &eqnedit::Equation::move_home)
        .def("move_end", &eqnedit::Equation::move_end)
        .def("begin_pointer_selection",
             &eqnedit::Equation::begin_pointer_selection,
             py::arg("x_points"), py::arg("y_points"),
             py::arg("style") = eqnedit::SvgStyle())
        .def("extend_pointer_selection",
             &eqnedit::Equation::extend_pointer_selection,
             py::arg("x_points"), py::arg("y_points"),
             py::arg("style") = eqnedit::SvgStyle())
        .def("end_pointer_selection",
             &eqnedit::Equation::end_pointer_selection)
        .def("hit_test", &eqnedit::Equation::hit_test,
             py::arg("x_points"), py::arg("y_points"),
             py::arg("style") = eqnedit::SvgStyle(), py::arg("extend") = false)
        .def("new_line", &eqnedit::Equation::new_line)
        .def("alignment_tab", &eqnedit::Equation::alignment_tab)
        .def("move_up", &eqnedit::Equation::move_up)
        .def("move_down", &eqnedit::Equation::move_down)
        .def("is_multiline", &eqnedit::Equation::is_multiline)
        .def("matrix_dimensions", [](const eqnedit::Equation& equation)
                -> py::object {
            int rows = 0, cols = 0;
            if (!equation.matrix_dimensions(&rows, &cols)) return py::none();
            return py::object(py::make_tuple(rows, cols));
        })
        .def("matrix_add_row", &eqnedit::Equation::matrix_add_row)
        .def("matrix_remove_row", &eqnedit::Equation::matrix_remove_row)
        .def("matrix_add_column", &eqnedit::Equation::matrix_add_column)
        .def("matrix_remove_column", &eqnedit::Equation::matrix_remove_column)
        .def("begin_selection", &eqnedit::Equation::begin_selection)
        .def("select_step_left", &eqnedit::Equation::select_step_left)
        .def("select_step_right", &eqnedit::Equation::select_step_right)
        .def("clear_selection", &eqnedit::Equation::clear_selection)
        .def("select_current_slot", &eqnedit::Equation::select_current_slot)
        .def("select_containing_structure",
             &eqnedit::Equation::select_containing_structure)
        .def("select_all", &eqnedit::Equation::select_all)
        .def("has_selection", &eqnedit::Equation::has_selection)
        .def("selection_latex", &eqnedit::Equation::selection_latex)
        .def("delete_selection", &eqnedit::Equation::delete_selection)
        .def("insert_text", &eqnedit::Equation::insert_text, py::arg("text"))
        .def("insert_styled_text", &eqnedit::Equation::insert_styled_text,
             py::arg("text"), py::arg("style"))
        .def("restyle_selection", &eqnedit::Equation::restyle_selection,
             py::arg("style"))
        .def("insert_latex", &eqnedit::Equation::insert_latex, py::arg("tex"))
        .def("insert_symbol", &eqnedit::Equation::insert_symbol, py::arg("command"))
        .def("insert_template", &eqnedit::Equation::insert_template, py::arg("kind"))
        .def("backspace", &eqnedit::Equation::backspace)
        .def("erase", &eqnedit::Equation::erase)
        .def("can_undo", &eqnedit::Equation::can_undo)
        .def("can_redo", &eqnedit::Equation::can_redo)
        .def("undo_name", &eqnedit::Equation::undo_name)
        .def("redo_name", &eqnedit::Equation::redo_name)
        .def("undo", &eqnedit::Equation::undo)
        .def("redo", &eqnedit::Equation::redo)
        .def("command", &eqnedit::Equation::command, py::arg("name"))
        .def("metrics", [](const eqnedit::Equation& equation,
                           const eqnedit::SvgStyle& style) {
            auto metrics = equation.metrics(style);
            return py::make_tuple(metrics.width, metrics.height, metrics.baseline);
        }, py::arg("style") = eqnedit::SvgStyle())
        .def("caret_geometry", [](const eqnedit::Equation& equation,
                                  const eqnedit::SvgStyle& style) -> py::object {
            eqnedit::CaretGeometry caret;
            if (!equation.caret_geometry(&caret, style)) return py::none();
            return py::object(py::make_tuple(
                caret.x, caret.top, caret.bottom));
        }, py::arg("style") = eqnedit::SvgStyle())
        .def_static("templates", &eqnedit::Equation::templates)
        .def_static("shortcuts", []() {
            py::list result;
            for (const auto& binding : eqnedit::Equation::shortcuts())
                result.append(py::make_tuple(binding.chord, binding.command,
                                             binding.label));
            return result;
        })
        .def_static("sequence_shortcuts", []() {
            py::list result;
            for (const auto& binding : eqnedit::Equation::sequence_shortcuts())
                result.append(py::make_tuple(std::string(1, binding.prefix),
                                             std::string(1, binding.key),
                                             binding.shift, binding.command,
                                             binding.label));
            return result;
        })
        .def_static("resolve_sequence", [](const std::string& prefix,
                                            const std::string& key,
                                            bool shift) {
            if (prefix.size() != 1 || key.size() != 1) return std::string();
            const auto* binding = eqnedit::Equation::resolve_sequence(
                prefix[0], key[0], shift);
            return binding ? std::string(binding->command) : std::string();
        }, py::arg("prefix"), py::arg("key"), py::arg("shift") = false);
}
