/************************************************************************//**
 * File: rad_ngsolve.cpp
 * Description: NGSolve CoefficientFunction binding for Radia
 * Project: Radia
 * First release: October 2025
 *
 * IMPORTANT:
 * - NGSolve uses meters, Radia uses millimeters
 * - Automatic unit conversion: m -> mm (multiply by 1000)
 *
 * Field types supported:
 * - 'b': Magnetic flux density (Tesla)
 * - 'h': Magnetic field (A/m)
 * - 'a': Vector potential (T*m)
 * - 'm': Magnetization (A/m)
 *
 * Coordinate transformation (v0.07):
 * - origin: Translation vector (meters)
 * - u_axis, v_axis, w_axis: Local coordinate system (auto-normalized)
 *
 * Batch evaluation (v0.08):
 * - Implements efficient batch evaluation for multiple points
 * - Reduces Python call overhead from O(M) to O(1)
 * - Enables OpenMP parallelization in Radia core
 *
 * @version 0.08
 ***************************************************************************/

#include <iostream>
#include <string>
#include <vector>
#include <fem.hpp>
#include <python_ngstd.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <cmath>

namespace py = pybind11;

namespace ngfem
{

class RadiaFieldCF : public CoefficientFunction
{
public:
	int radia_obj;
	std::string field_type;

	// Coordinate transformation
	double origin[3];      // Translation vector [m]
	double u_axis[3];      // Local u-axis (normalized)
	double v_axis[3];      // Local v-axis (normalized)
	double w_axis[3];      // Local w-axis (normalized)
	bool use_transform;    // Whether to apply coordinate transformation

	// Computation settings
	py::object precision;  // Computation precision (None = use Radia default)
	py::object use_hmatrix; // Use H-matrix acceleration (None = keep current setting)

	RadiaFieldCF(int obj, const std::string& ftype = "b",
	             py::object py_origin = py::none(),
	             py::object py_u = py::none(),
	             py::object py_v = py::none(),
	             py::object py_w = py::none(),
	             py::object py_precision = py::none(),
	             py::object py_use_hmatrix = py::none())
	    : CoefficientFunction(3), radia_obj(obj), field_type(ftype), use_transform(false),
	      precision(py_precision), use_hmatrix(py_use_hmatrix)
	{
		// Validate field type
		if (field_type != "b" && field_type != "h" &&
		    field_type != "a" && field_type != "m") {
			throw std::invalid_argument(
				"Invalid field_type. Must be 'b' (flux density), "
				"'h' (magnetic field), 'a' (vector potential), or 'm' (magnetization)");
		}

		// Default: identity transformation
		origin[0] = 0.0; origin[1] = 0.0; origin[2] = 0.0;
		u_axis[0] = 1.0; u_axis[1] = 0.0; u_axis[2] = 0.0;
		v_axis[0] = 0.0; v_axis[1] = 1.0; v_axis[2] = 0.0;
		w_axis[0] = 0.0; w_axis[1] = 0.0; w_axis[2] = 1.0;

		// Parse origin
		if (!py_origin.is_none()) {
			parse_vector(py_origin, origin);
			use_transform = true;
		}

		// Parse and normalize u-axis
		if (!py_u.is_none()) {
			parse_vector(py_u, u_axis);
			normalize(u_axis);
			use_transform = true;
		}

		// Parse and normalize v-axis
		if (!py_v.is_none()) {
			parse_vector(py_v, v_axis);
			normalize(v_axis);
			use_transform = true;
		}

		// Parse and normalize w-axis
		if (!py_w.is_none()) {
			parse_vector(py_w, w_axis);
			normalize(w_axis);
			use_transform = true;
		}

		// Apply computation settings
		py::gil_scoped_acquire acquire;
		py::module_ rad = py::module_::import("radia");

		// Set H-matrix usage if specified
		if (!py_use_hmatrix.is_none()) {
			bool enable_hmatrix = py_use_hmatrix.cast<bool>();
			if (enable_hmatrix) {
				rad.attr("SolverHMatrixEnable")();
			} else {
				rad.attr("SolverHMatrixDisable")();
			}
		}

		// Set precision if specified
		if (!py_precision.is_none()) {
			double prec = py_precision.cast<double>();
			// Set precision for all field computation types
			std::string prec_str = "PrcB->" + std::to_string(prec) +
			                       ",PrcA->" + std::to_string(prec) +
			                       ",PrcH->" + std::to_string(prec) +
			                       ",PrcM->" + std::to_string(prec);
			rad.attr("FldCmpPrc")(prec_str);
		}
	}

private:
	void parse_vector(py::object py_vec, double vec[3]) {
		if (py::isinstance<py::list>(py_vec)) {
			py::list lst = py_vec.cast<py::list>();
			if (lst.size() != 3) {
				throw std::invalid_argument("Vector must have 3 components");
			}
			vec[0] = lst[0].cast<double>();
			vec[1] = lst[1].cast<double>();
			vec[2] = lst[2].cast<double>();
		} else if (py::isinstance<py::tuple>(py_vec)) {
			py::tuple tup = py_vec.cast<py::tuple>();
			if (tup.size() != 3) {
				throw std::invalid_argument("Vector must have 3 components");
			}
			vec[0] = tup[0].cast<double>();
			vec[1] = tup[1].cast<double>();
			vec[2] = tup[2].cast<double>();
		} else {
			throw std::invalid_argument("Vector must be a list or tuple");
		}
	}

	void normalize(double vec[3]) {
		double norm = std::sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
		if (norm < 1e-12) {
			throw std::invalid_argument("Cannot normalize zero vector");
		}
		vec[0] /= norm;
		vec[1] /= norm;
		vec[2] /= norm;
	}

	double dot(const double a[3], const double b[3]) const {
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

public:

	virtual ~RadiaFieldCF() {
		// Acquire GIL before releasing Python objects to prevent memory leaks
		// When NGSolve destroys this CoefficientFunction, we must ensure
		// Python reference counting is done safely
		py::gil_scoped_acquire acquire;
		precision.release();
		use_hmatrix.release();
	}

	// Scalar evaluation - not used for vector fields, return 0
	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
	    return 0.0;
	}

	// Single point evaluation (for direct calls like cf(mip))
	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                     FlatVector<> result) const override
	{
	    py::gil_scoped_acquire acquire;

	    try {
	        auto pnt = mip.GetPoint();
	        int dim = pnt.Size();

	        // Get global coordinates (NGSolve, in meters)
	        double p_global[3];
	        p_global[0] = pnt[0];
	        p_global[1] = (dim >= 2) ? pnt[1] : 0.0;
	        p_global[2] = (dim >= 3) ? pnt[2] : 0.0;

	        // Apply coordinate transformation if enabled
	        double p_local[3];
	        if (use_transform) {
	            double p_translated[3];
	            p_translated[0] = p_global[0] - origin[0];
	            p_translated[1] = p_global[1] - origin[1];
	            p_translated[2] = p_global[2] - origin[2];

	            p_local[0] = dot(u_axis, p_translated);
	            p_local[1] = dot(v_axis, p_translated);
	            p_local[2] = dot(w_axis, p_translated);
	        } else {
	            p_local[0] = p_global[0];
	            p_local[1] = p_global[1];
	            p_local[2] = p_global[2];
	        }

	        // Convert m -> mm
	        py::list coords;
	        coords.append(p_local[0] * 1000.0);
	        coords.append(p_local[1] * 1000.0);
	        coords.append(p_local[2] * 1000.0);

	        py::module_ rad = py::module_::import("radia");
	        py::object field_result = rad.attr("Fld")(radia_obj, field_type, coords);

	        py::list field_list = field_result.cast<py::list>();

	        double f_local[3];
	        f_local[0] = field_list[0].cast<double>();
	        f_local[1] = field_list[1].cast<double>();
	        f_local[2] = field_list[2].cast<double>();

	        // Transform field back to global coordinate system
	        double f_global[3];
	        if (use_transform) {
	            f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
	            f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
	            f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
	        } else {
	            f_global[0] = f_local[0];
	            f_global[1] = f_local[1];
	            f_global[2] = f_local[2];
	        }

	        double scale = (field_type == "a") ? 0.001 : 1.0;

	        result(0) = f_global[0] * scale;
	        result(1) = f_global[1] * scale;
	        result(2) = f_global[2] * scale;

	    } catch (std::exception &e) {
	        std::cerr << "[RadiaField] Single point error (" << field_type << "): "
	                  << e.what() << std::endl;
	        result(0) = 0.0;
	        result(1) = 0.0;
	        result(2) = 0.0;
	    }
	}

	// Batch evaluation: evaluate field at multiple points in one call
	// This significantly reduces Python call overhead and enables OpenMP parallelization
	virtual void Evaluate(const BaseMappedIntegrationRule& mir,
	                     BareSliceMatrix<> result) const override
	{
	    py::gil_scoped_acquire acquire;

	    try {
	        size_t npts = mir.Size();

	        // Collect all points in a Python list of lists
	        py::list points_list;

	        for (size_t i = 0; i < npts; i++) {
	            auto pnt = mir[i].GetPoint();
	            int dim = pnt.Size();

	            // Get global coordinates (NGSolve, in meters)
	            double p_global[3];
	            p_global[0] = pnt[0];
	            p_global[1] = (dim >= 2) ? pnt[1] : 0.0;
	            p_global[2] = (dim >= 3) ? pnt[2] : 0.0;

	            // Apply coordinate transformation if enabled
	            double p_local[3];
	            if (use_transform) {
	                double p_translated[3];
	                p_translated[0] = p_global[0] - origin[0];
	                p_translated[1] = p_global[1] - origin[1];
	                p_translated[2] = p_global[2] - origin[2];

	                p_local[0] = dot(u_axis, p_translated);
	                p_local[1] = dot(v_axis, p_translated);
	                p_local[2] = dot(w_axis, p_translated);
	            } else {
	                p_local[0] = p_global[0];
	                p_local[1] = p_global[1];
	                p_local[2] = p_global[2];
	            }

	            // Convert m -> mm
	            py::list coords;
	            coords.append(p_local[0] * 1000.0);
	            coords.append(p_local[1] * 1000.0);
	            coords.append(p_local[2] * 1000.0);

	            points_list.append(coords);
	        }

	        // Single Python call for all points!
	        py::module_ rad = py::module_::import("radia");
	        py::object field_results = rad.attr("Fld")(radia_obj, field_type, points_list);

	        // field_results is a list of [Bx, By, Bz] for each point
	        py::list results_list = field_results.cast<py::list>();

	        // Extract results and apply transformations
	        double scale = (field_type == "a") ? 0.001 : 1.0;

	        for (size_t i = 0; i < npts; i++) {
	            py::list field_list = results_list[i].cast<py::list>();

	            double f_local[3];
	            f_local[0] = field_list[0].cast<double>();
	            f_local[1] = field_list[1].cast<double>();
	            f_local[2] = field_list[2].cast<double>();

	            // Transform field back to global coordinate system
	            double f_global[3];
	            if (use_transform) {
	                f_global[0] = u_axis[0]*f_local[0] + v_axis[0]*f_local[1] + w_axis[0]*f_local[2];
	                f_global[1] = u_axis[1]*f_local[0] + v_axis[1]*f_local[1] + w_axis[1]*f_local[2];
	                f_global[2] = u_axis[2]*f_local[0] + v_axis[2]*f_local[1] + w_axis[2]*f_local[2];
	            } else {
	                f_global[0] = f_local[0];
	                f_global[1] = f_local[1];
	                f_global[2] = f_local[2];
	            }

	            result(0, i) = f_global[0] * scale;
	            result(1, i) = f_global[1] * scale;
	            result(2, i) = f_global[2] * scale;
	        }

	    } catch (std::exception &e) {
	        std::cerr << "[RadiaField] Batch evaluation error (" << field_type << "): "
	                  << e.what() << std::endl;
	        // Fill with zeros on error
	        size_t npts = mir.Size();
	        for (size_t i = 0; i < npts; i++) {
	            result(i, 0) = 0.0;
	            result(i, 1) = 0.0;
	            result(i, 2) = 0.0;
	        }
	    }
	}
};

} // namespace ngfem

PYBIND11_MODULE(rad_ngsolve, m) {
	m.doc() = "NGSolve CoefficientFunction interface for Radia (with m->mm conversion and coordinate transformation)";

	// Unified interface with coordinate transformation
	py::class_<ngfem::RadiaFieldCF,
	           std::shared_ptr<ngfem::RadiaFieldCF>,
	           ngfem::CoefficientFunction>(m, "RadiaField")
	    .def(py::init<int>(), py::arg("radia_obj"),
	         "Create Radia field CoefficientFunction (default: magnetic flux density)")
	    .def(py::init<int, const std::string&>(),
	         py::arg("radia_obj"), py::arg("field_type"),
	         "Create Radia field CoefficientFunction\n"
	         "field_type: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)")
	    .def(py::init<int, const std::string&, py::object, py::object, py::object, py::object, py::object, py::object>(),
	         py::arg("radia_obj"),
	         py::arg("field_type") = "b",
	         py::arg("origin") = py::none(),
	         py::arg("u_axis") = py::none(),
	         py::arg("v_axis") = py::none(),
	         py::arg("w_axis") = py::none(),
	         py::arg("precision") = py::none(),
	         py::arg("use_hmatrix") = py::none(),
	         "Create Radia field CoefficientFunction with full control\n\n"
	         "Parameters:\n"
	         "  radia_obj: Radia object ID\n"
	         "  field_type: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)\n"
	         "  origin: Translation vector [x, y, z] in meters (default: [0, 0, 0])\n"
	         "  u_axis: Local u-axis [ux, uy, uz] (default: [1, 0, 0]) - will be normalized\n"
	         "  v_axis: Local v-axis [vx, vy, vz] (default: [0, 1, 0]) - will be normalized\n"
	         "  w_axis: Local w-axis [wx, wy, wz] (default: [0, 0, 1]) - will be normalized\n"
	         "  precision: Computation precision in Tesla (default: None = Radia default)\n"
	         "  use_hmatrix: Enable H-matrix acceleration (default: None = keep current setting)\n\n"
	         "Coordinate transformation:\n"
	         "  1. Global point p is translated by origin: p' = p - origin\n"
	         "  2. p' is projected onto local axes: p_local = [u*p', v*p', w*p']\n"
	         "  3. Field is calculated in Radia's coordinate system\n"
	         "  4. Field is transformed back to global: F = u*F_local[0] + v*F_local[1] + w*F_local[2]\n\n"
	         "Performance control:\n"
	         "  - precision: Controls accuracy vs speed trade-off (smaller = more accurate)\n"
	         "  - use_hmatrix: True = fast for large N, False = accurate for small N\n\n"
	         "Example:\n"
	         "  # High accuracy, no H-matrix\n"
	         "  B_cf = rad_ngsolve.RadiaField(magnet, 'b', precision=1e-6, use_hmatrix=False)\n\n"
	         "  # Fast evaluation with H-matrix\n"
	         "  B_cf = rad_ngsolve.RadiaField(magnet, 'b', use_hmatrix=True)")
	    .def_readonly("radia_obj", &ngfem::RadiaFieldCF::radia_obj)
	    .def_readonly("field_type", &ngfem::RadiaFieldCF::field_type)
	    .def_readonly("use_transform", &ngfem::RadiaFieldCF::use_transform)
	    .def_readonly("precision", &ngfem::RadiaFieldCF::precision)
	    .def_readonly("use_hmatrix", &ngfem::RadiaFieldCF::use_hmatrix);
}
