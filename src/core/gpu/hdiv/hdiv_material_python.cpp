// Thin Python owner/publisher for the shared NGSolve material staging boundary.
#include "hdiv_material_staging.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <mutex>
namespace py=pybind11;
namespace {
using namespace radia::hdiv_cuda;
class MaterialOwner {
    std::mutex mutex;
    std::unique_ptr<MaterialStagingContext> context;
public:
    MaterialOwner(const std::string& path,int order):context(std::make_unique<MaterialStagingContext>(path,order)){}
    MaterialQuadraturePayload quadrature(std::size_t first,std::size_t count){
        std::lock_guard<std::mutex> lock(mutex);
        if(!context)throw std::invalid_argument("Material context is closed");
        return stageMaterialQuadrature(*context,first,count);
    }
    void close(){std::lock_guard<std::mutex> lock(mutex);context.reset();}
};
template<class T> py::array_t<T> ownedArray(const std::vector<T>& values){
    py::array_t<T> result(values.size());
    std::copy(values.begin(),values.end(),result.mutable_data());return result;
}
py::dict publish(const MaterialQuadraturePayload& p){
    py::dict out;
    out["schema"]="radia.hdiv-material-quadrature.v1";
    out["evaluation_column_offsets"]=ownedArray(p.columnOffsets);
    out["evaluation_row_indices"]=ownedArray(p.rowIndices);
    out["evaluation_values"]=ownedArray(p.values);
    out["evaluation_rows"]=3*p.weights.size();out["dofs"]=p.dofs;
    out["weights"]=ownedArray(p.weights);out["element"]=ownedArray(p.elements);
    auto points=ownedArray(p.points);points.resize({py::ssize_t(p.weights.size()),py::ssize_t(3)});out["points"]=points;
    out["global_point_indices"]=ownedArray(p.globalPointIndices);
    out["integration_order"]=p.integrationOrder;out["total_points"]=p.totalPoints;
    out["total_elements"]=p.totalElements;out["first_element"]=p.firstElement;out["element_count"]=p.elementCount;
    out["potential_shape_entries"]=p.potentialShapeEntries;out["tangent_nonzero_bound"]=p.tangentNonzeroBound;
    out["epoch"]=1; // Immutable mesh/space; no in-place update operation exists.
    return out;
}
}
void registerMaterialStaging(py::module_& module){
    py::class_<MaterialOwner>(module,"_MaterialStagingContext")
        .def(py::init<const std::string&,int>(),py::arg("mesh_path"),py::arg("order"),py::call_guard<py::gil_scoped_release>())
        .def("quadrature",[](MaterialOwner& self,std::size_t first,std::size_t count){
            MaterialQuadraturePayload payload;
            {py::gil_scoped_release release;payload=self.quadrature(first,count);}
            return publish(payload);
        },py::arg("first_element")=0,py::arg("element_count")=std::numeric_limits<std::size_t>::max())
        .def("close",&MaterialOwner::close,py::call_guard<py::gil_scoped_release>());
}
