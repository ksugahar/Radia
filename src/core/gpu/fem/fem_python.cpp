/* SPDX-License-Identifier: MPL-2.0 */
// Python binding of the FEM device route. It checks every array against the
// declared shape and passes pointers on; all numerics are on the device.
#include "fem_api.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace radia_fem;

namespace {

// Holds a reference to every accepted array for the duration of one call.
class Arrays {
public:
    template <class T> const T* get(py::handle value, long long expected, const char* name) {
        if (value.is_none()) return nullptr;
        if (!py::isinstance<py::array_t<T>>(value))
            throw py::value_error(std::string(name) + " has the wrong type or dtype");
        auto array = py::reinterpret_borrow<py::array>(value);
        if (!array.dtype().is(py::dtype::of<T>()) || !(array.flags() & py::array::c_style))
            throw py::value_error(std::string(name) + " must be a C-contiguous array of the exact dtype");
        if (array.size() != expected)
            throw py::value_error(std::string(name) + " holds " + std::to_string(array.size()) + " values, expected " + std::to_string(expected));
        held_.push_back(array);
        return static_cast<const T*>(array.data());
    }

private:
    std::vector<py::array> held_;
};

FormInput form_input(Arrays& arrays, py::handle value, const GroupInput& group, const char* name) {
    FormInput form;
    if (value.is_none()) return form;
    py::dict d = py::reinterpret_borrow<py::dict>(value);
    form.law = py::cast<int>(d["law"]);
    form.nq = py::cast<int>(d["nq"]);
    form.affine = py::cast<bool>(d["affine"]);
    const long long points = form.affine ? 1 : form.nq;
    form.weights = arrays.get<double>(d["weights"], form.nq, name);
    form.table = arrays.get<double>(d["table"], (long long)(group.ncls) * form.nq * group.nl * 3, name);
    form.dgeo = arrays.get<double>(d["dgeo"], points * group.ngeo * group.dim, name);
    form.jacobian = arrays.get<double>(d["jacobian"], group.ne * points * 3 * group.dim, name);
    form.tensor = arrays.get<double>(d["tensor"], group.ne * points * 6, name);
    return form;
}

std::vector<GroupInput> group_inputs(Arrays& arrays, const py::list& groups) {
    std::vector<GroupInput> inputs;
    for (py::handle item : groups) {
        py::dict d = py::reinterpret_borrow<py::dict>(item);
        GroupInput group;
        group.dim = py::cast<int>(d["dim"]);
        group.ne = py::cast<int>(d["ne"]);
        group.nl = py::cast<int>(d["nl"]);
        group.ncls = py::cast<int>(d["ncls"]);
        group.ngeo = py::cast<int>(d["ngeo"]);
        group.dofs = arrays.get<int>(d["dofs"], (long long)(group.ne) * group.nl, "dofs");
        group.cls = arrays.get<int>(d["cls"], group.ne, "cls");
        group.nodes = arrays.get<double>(d["nodes"], (long long)(group.ne) * group.ngeo * 3, "nodes");
        group.coef = arrays.get<double>(d["coef"], (long long)(group.ne) * COLUMNS, "coef");
        group.stiff = form_input(arrays, d["stiff"], group, "stiff form");
        group.mass = form_input(arrays, d["mass"], group, "mass form");
        inputs.push_back(group);
    }
    return inputs;
}

py::dict amg_dict(const AmgReport& report) {
    py::dict d;
    d["rows"] = report.rows;
    d["nonzeros"] = report.nonzeros;
    d["pinned"] = report.pinned;
    d["coarse"] = report.coarse;
    return d;
}

py::array_t<double> vector_of(int n) { return py::array_t<double>(py::ssize_t(n)); }

class PySystem {
public:
    PySystem(int ndof, bool complex_values, const py::list& groups, const py::array& free_dofs) : n_(ndof), complex_(complex_values) {
        Arrays arrays;
        std::vector<GroupInput> inputs = group_inputs(arrays, groups);
        for (const GroupInput& input : inputs) tensor_values_.push_back(input.stiff.tensor ? (long long)(input.ne) * input.stiff.nq * 6 : 0);
        const int* mask = arrays.get<int>(free_dofs, ndof, "free_dofs");
        py::gil_scoped_release release;
        system_ = std::make_unique<System>(ndof, complex_values, inputs, mask);
    }
    py::dict setup_amg(int min_coarse, double theta) {
        AmgReport report;
        {
            py::gil_scoped_release release;
            report = system_->setup_amg(min_coarse, theta);
        }
        return amg_dict(report);
    }
    py::list setup_ams(int nedge, int nv, const py::array& edges, const py::array& coordinates, int min_coarse, double theta) {
        Arrays arrays;
        const int* edge_vertices = arrays.get<int>(edges, (long long)(nedge) * 2, "edge_vertices");
        const double* coords = arrays.get<double>(coordinates, (long long)(nv) * 3, "coordinates");
        if (!edge_vertices || !coords) throw py::value_error("AMS needs the mesh edges and coordinates");
        std::vector<AmgReport> reports;
        {
            py::gil_scoped_release release;
            reports = system_->setup_ams(nedge, nv, edge_vertices, coords, min_coarse, theta, AmsOptions{});
        }
        py::list out;
        for (const AmgReport& report : reports) out.append(amg_dict(report));
        return out;
    }
    py::list setup_bddc(const py::list& split, const py::array& wb_of_dof, int nw, int nedge, int nv, const py::object& edges,
                        const py::object& coordinates, int cycles, int min_coarse, double theta) {
        Arrays arrays;
        std::vector<BddcGroupInput> inputs;
        for (py::handle item : split) {
            py::dict d = py::reinterpret_borrow<py::dict>(item);
            BddcGroupInput input;
            input.nw = py::cast<int>(d["nw"]);
            input.ni = py::cast<int>(d["ni"]);
            input.wb_local = arrays.get<int>(d["wb_local"], input.nw, "wb_local");
            input.int_local = arrays.get<int>(d["int_local"], input.ni, "int_local");
            inputs.push_back(input);
        }
        const int* wb = arrays.get<int>(wb_of_dof, n_, "wb_of_dof");
        const int* edge_vertices = arrays.get<int>(edges, (long long)(nedge) * 2, "edge_vertices");
        const double* coords = arrays.get<double>(coordinates, (long long)(nv) * 3, "coordinates");
        std::vector<AmgReport> reports;
        {
            py::gil_scoped_release release;
            reports = system_->setup_bddc(inputs, wb, nw, nedge, nv, edge_vertices, coords, cycles, min_coarse, theta);
        }
        py::list out;
        for (const AmgReport& report : reports) out.append(amg_dict(report));
        return out;
    }
    void update_tensors(const py::list& tensors) {
        if (py::len(tensors) != tensor_values_.size()) throw py::value_error("one tensor entry per element group is required");
        Arrays arrays;
        std::vector<const double*> pointers;
        for (size_t g = 0; g < tensor_values_.size(); ++g) {
            if (!tensor_values_[g] && !tensors[g].is_none()) throw py::value_error("this group was not built with a point tensor");
            pointers.push_back(arrays.get<double>(tensors[g], tensor_values_[g], "tensor"));
        }
        py::gil_scoped_release release;
        system_->update_tensors(pointers);
    }
    py::array_t<double> probe_action(int which, const py::array& x) {
        Arrays arrays;
        const double* in = arrays.get<double>(x, n_, "x");
        if (!in) throw py::value_error("probe input is required");
        py::array_t<double> y = vector_of(n_);
        double* out = y.mutable_data();
        {
            py::gil_scoped_release release;
            system_->probe_action(which, in, out);
        }
        return y;
    }
    py::tuple probe_preconditioner(const py::array& x_re, const py::object& x_im) {
        Arrays arrays;
        const double* re = arrays.get<double>(x_re, n_, "x_re");
        const double* im = arrays.get<double>(x_im, n_, "x_im");
        py::array_t<double> y_re = vector_of(n_), y_im = vector_of(complex_ ? n_ : 0);
        double* out_re = y_re.mutable_data();
        double* out_im = y_im.mutable_data();
        {
            py::gil_scoped_release release;
            system_->probe_preconditioner(re, im, out_re, out_im);
        }
        return py::make_tuple(y_re, y_im);
    }
    py::array_t<double> probe_ams_transfer(int space, bool transpose, const py::array& x) {
        const std::pair<int, int> shape = system_->ams_transfer_shape();
        Arrays arrays;
        const double* in = arrays.get<double>(x, transpose ? shape.first : shape.second, "x");
        if (!in) throw py::value_error("probe input is required");
        py::array_t<double> y = vector_of(transpose ? shape.second : shape.first);
        double* out = y.mutable_data();
        {
            py::gil_scoped_release release;
            system_->probe_ams_transfer(space, transpose, in, out);
        }
        return y;
    }
    py::dict solve(const py::array& rhs_re, const py::object& rhs_im, const py::object& lift_re, const py::object& lift_im,
                   double tolerance, int max_iterations, double relative_limit) {
        Arrays arrays;
        const double* f_re = arrays.get<double>(rhs_re, n_, "rhs_re");
        const double* f_im = arrays.get<double>(rhs_im, n_, "rhs_im");
        const double* g_re = arrays.get<double>(lift_re, n_, "lift_re");
        const double* g_im = arrays.get<double>(lift_im, n_, "lift_im");
        if (!f_re) throw py::value_error("a right-hand side is required");
        py::array_t<double> x_re = vector_of(n_), x_im = vector_of(complex_ ? n_ : 0);
        double* out_re = x_re.mutable_data();
        double* out_im = x_im.mutable_data();
        SolveReport report;
        {
            py::gil_scoped_release release;
            report = system_->solve(f_re, f_im, g_re, g_im, tolerance, max_iterations, relative_limit, out_re, out_im);
        }
        py::dict out;
        out["x_re"] = x_re;
        out["x_im"] = x_im;
        out["iterations"] = report.iterations;
        out["true_residual"] = report.true_residual;
        out["initial_norm"] = report.initial_norm;
        out["final_norm"] = report.final_norm;
        out["history"] = report.history;
        return out;
    }
    py::list ledger() const {
        py::list out;
        for (const LedgerRow& row : system_->ledger()) {
            py::dict d;
            d["payload"] = row.payload;
            d["direction"] = row.direction;
            d["bytes"] = row.bytes;
            out.append(d);
        }
        return out;
    }
    int nonzeros() const { return system_->nonzeros(); }
    uint64_t peak_device_bytes() const { return system_->peak_device_bytes(); }
    py::dict timings() const {
        py::dict d;
        for (const auto& item : system_->timings()) d[py::str(item.first)] = item.second;
        return d;
    }

private:
    int n_ = 0;
    bool complex_ = false;
    std::vector<long long> tensor_values_;
    std::unique_ptr<System> system_;
};

}  // namespace

PYBIND11_MODULE(_fem_cuda, m) {
    m.doc() = "Device-built FEM operators, AMG, AMS and BDDC. No matrix crosses this boundary.";
    m.def("device_name", &device_name);
    m.attr("LAW_NONE") = int(LAW_NONE);
    m.attr("LAW_GRAD") = int(LAW_GRAD);
    m.attr("LAW_VALUE") = int(LAW_VALUE);
    m.attr("LAW_CURL") = int(LAW_CURL);
    m.attr("LAW_COVARIANT") = int(LAW_COVARIANT);
    py::class_<PySystem>(m, "System")
        .def(py::init<int, bool, const py::list&, const py::array&>())
        .def("setup_amg", &PySystem::setup_amg)
        .def("setup_ams", &PySystem::setup_ams)
        .def("setup_bddc", &PySystem::setup_bddc)
        .def("update_tensors", &PySystem::update_tensors)
        .def("probe_action", &PySystem::probe_action)
        .def("probe_preconditioner", &PySystem::probe_preconditioner)
        .def("probe_ams_transfer", &PySystem::probe_ams_transfer)
        .def("solve", &PySystem::solve)
        .def("ledger", &PySystem::ledger)
        .def("nonzeros", &PySystem::nonzeros)
        .def("peak_device_bytes", &PySystem::peak_device_bytes)
        .def("timings", &PySystem::timings);
}
