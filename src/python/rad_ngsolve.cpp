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
 * - 'a': Vector potential (T·m)
 * - 'm': Magnetization (A/m)
 *
 * @version 0.06
 ***************************************************************************/

#include <iostream>
#include <string>
#include <fem.hpp>
#include <python_ngstd.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <cmath>

namespace py = pybind11;

namespace ngfem
{

class RadiaFieldCF : public CoefficientFunction
{
public:
	int radia_obj;
	std::string field_type;

	RadiaFieldCF(int obj, const std::string& ftype = "b")
	    : CoefficientFunction(3), radia_obj(obj), field_type(ftype)
	{
		// Validate field type
		if (field_type != "b" && field_type != "h" &&
		    field_type != "a" && field_type != "m") {
			throw std::invalid_argument(
				"Invalid field_type. Must be 'b' (flux density), "
				"'h' (magnetic field), 'a' (vector potential), or 'm' (magnetization)");
		}
	}

	virtual ~RadiaFieldCF() {}

	virtual double Evaluate(const BaseMappedIntegrationPoint& mip) const override
	{
	    return 0.0;
	}

	virtual void Evaluate(const BaseMappedIntegrationPoint& mip,
	                     FlatVector<> result) const override
	{
	    py::gil_scoped_acquire acquire;

	    try {
	        auto pnt = mip.GetPoint();
	        int dim = pnt.Size();

	        // NGSolve uses meters, Radia uses millimeters
	        // Convert m -> mm by multiplying by 1000
	        py::list coords;
	        coords.append(pnt[0] * 1000.0);
	        coords.append((dim >= 2) ? pnt[1] * 1000.0 : 0.0);
	        coords.append((dim >= 3) ? pnt[2] * 1000.0 : 0.0);

	        py::module_ rad = py::module_::import("radia");
	        py::object field_result = rad.attr("Fld")(radia_obj, field_type, coords);

	        py::list field_list = field_result.cast<py::list>();

	        // Coordinate system conversion for vector potential 'a':
	        // NGSolve differentiates in meters, Radia uses millimeters
	        // curl(A) = ∂A/∂x_m = (∂A/∂x_mm) × (∂x_mm/∂x_m) = (∂A/∂x_mm) × 1000
	        // To get correct B = curl(A), we need to divide A by 1000
	        double scale = (field_type == "a") ? 0.001 : 1.0;

	        result(0) = field_list[0].cast<double>() * scale;
	        result(1) = field_list[1].cast<double>() * scale;
	        result(2) = field_list[2].cast<double>() * scale;

	    } catch (std::exception &e) {
	        std::cerr << "[RadiaField] Error (" << field_type << "): "
	                  << e.what() << std::endl;
	        result(0) = 0.0;
	        result(1) = 0.0;
	        result(2) = 0.0;
	    }
	}
};

} // namespace ngfem

PYBIND11_MODULE(rad_ngsolve, m) {
	m.doc() = "NGSolve CoefficientFunction interface for Radia (with m->mm conversion)";

	// Unified interface
	py::class_<ngfem::RadiaFieldCF,
	           std::shared_ptr<ngfem::RadiaFieldCF>,
	           ngfem::CoefficientFunction>(m, "RadiaField")
	    .def(py::init<int>(), py::arg("radia_obj"),
	         "Create Radia field CoefficientFunction (default: magnetic flux density)")
	    .def(py::init<int, const std::string&>(),
	         py::arg("radia_obj"), py::arg("field_type"),
	         "Create Radia field CoefficientFunction\n"
	         "field_type: 'b' (flux density), 'h' (field), 'a' (vector potential), 'm' (magnetization)")
	    .def_readonly("radia_obj", &ngfem::RadiaFieldCF::radia_obj)
	    .def_readonly("field_type", &ngfem::RadiaFieldCF::field_type);
}
