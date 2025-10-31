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
 * @authors Claude Code
 * @version 0.04
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

class RadiaBFieldCF : public CoefficientFunction
{
public:
    int radia_obj;
    std::string field_comp;

    RadiaBFieldCF(int obj, const std::string& comp = "b")
        : CoefficientFunction(3), radia_obj(obj), field_comp(comp)
    {
    }

    virtual ~RadiaBFieldCF() {}

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
            py::object B_result = rad.attr("Fld")(radia_obj, field_comp, coords);

            py::list B_list = B_result.cast<py::list>();
            result(0) = B_list[0].cast<double>();
            result(1) = B_list[1].cast<double>();
            result(2) = B_list[2].cast<double>();

        } catch (std::exception &e) {
            std::cerr << "[RadBfield] Error: " << e.what() << std::endl;
            result(0) = 0.0;
            result(1) = 0.0;
            result(2) = 0.0;
        }
    }
};

class RadiaHFieldCF : public CoefficientFunction
{
public:
    int radia_obj;

    RadiaHFieldCF(int obj) : CoefficientFunction(3), radia_obj(obj) {}
    virtual ~RadiaHFieldCF() {}

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

            py::list coords;
            coords.append(pnt[0] * 1000.0);
            coords.append((dim >= 2) ? pnt[1] * 1000.0 : 0.0);
            coords.append((dim >= 3) ? pnt[2] * 1000.0 : 0.0);

            py::module_ rad = py::module_::import("radia");
            py::object H_result = rad.attr("Fld")(radia_obj, "h", coords);
            py::list H_list = H_result.cast<py::list>();

            result(0) = H_list[0].cast<double>();
            result(1) = H_list[1].cast<double>();
            result(2) = H_list[2].cast<double>();

        } catch (std::exception &e) {
            std::cerr << "[RadHfield] Error: " << e.what() << std::endl;
            result(0) = 0.0;
            result(1) = 0.0;
            result(2) = 0.0;
        }
    }
};

class RadiaAFieldCF : public CoefficientFunction
{
public:
    int radia_obj;

    RadiaAFieldCF(int obj) : CoefficientFunction(3), radia_obj(obj) {}
    virtual ~RadiaAFieldCF() {}

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

            py::list coords;
            coords.append(pnt[0] * 1000.0);
            coords.append((dim >= 2) ? pnt[1] * 1000.0 : 0.0);
            coords.append((dim >= 3) ? pnt[2] * 1000.0 : 0.0);

            py::module_ rad = py::module_::import("radia");
            py::object A_result = rad.attr("Fld")(radia_obj, "a", coords);
            py::list A_list = A_result.cast<py::list>();

            result(0) = A_list[0].cast<double>();
            result(1) = A_list[1].cast<double>();
            result(2) = A_list[2].cast<double>();

        } catch (std::exception &e) {
            std::cerr << "[RadAfield] Error: " << e.what() << std::endl;
            result(0) = 0.0;
            result(1) = 0.0;
            result(2) = 0.0;
        }
    }
};

} // namespace ngfem

PYBIND11_MODULE(rad_ngsolve, m) {
    m.doc() = "NGSolve CoefficientFunction interface for Radia (with m->mm conversion)";

    py::class_<ngfem::RadiaBFieldCF,
               std::shared_ptr<ngfem::RadiaBFieldCF>,
               ngfem::CoefficientFunction>(m, "RadBfield")
        .def(py::init<int>(), py::arg("radia_obj"))
        .def(py::init<int, const std::string&>(),
             py::arg("radia_obj"), py::arg("field_comp"))
        .def_readonly("radia_obj", &ngfem::RadiaBFieldCF::radia_obj)
        .def_readonly("field_comp", &ngfem::RadiaBFieldCF::field_comp);

    py::class_<ngfem::RadiaHFieldCF,
               std::shared_ptr<ngfem::RadiaHFieldCF>,
               ngfem::CoefficientFunction>(m, "RadHfield")
        .def(py::init<int>(), py::arg("radia_obj"))
        .def_readonly("radia_obj", &ngfem::RadiaHFieldCF::radia_obj);

    py::class_<ngfem::RadiaAFieldCF,
               std::shared_ptr<ngfem::RadiaAFieldCF>,
               ngfem::CoefficientFunction>(m, "RadAfield")
        .def(py::init<int>(), py::arg("radia_obj"))
        .def_readonly("radia_obj", &ngfem::RadiaAFieldCF::radia_obj);
}
