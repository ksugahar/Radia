// Thin Python owner/publisher for the shared NGSolve material staging boundary.
#include "hdiv_charge_staging.hpp"
#include "hdiv_mixed_staging.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <mutex>
namespace py=pybind11;
namespace {
using namespace radia::hdiv_cuda;
class MaterialOwner {
    std::mutex mutex;
    std::unique_ptr<MaterialStagingContext> context;
public:
    MaterialOwner(const std::string& path,int order,int curveOrder,const std::vector<std::string>& periodicBoundaries):context(std::make_unique<MaterialStagingContext>(path,order,curveOrder,periodicBoundaries)){
        if(!periodicBoundaries.empty()&&!isHexMesh(*context))throw std::invalid_argument("Periodic material staging requires a pure HEX mesh");}
    MaterialQuadraturePayload quadrature(std::size_t first,std::size_t count){
        std::lock_guard<std::mutex> lock(mutex);
        if(!context)throw std::invalid_argument("Material context is closed");
        return stageMaterialQuadrature(*context,first,count);
    }
    FlatEvaluationTables flatTables(){
        std::lock_guard<std::mutex> lock(mutex);
        if(!context)throw std::invalid_argument("Material context is closed");
        return isHexMesh(*context)?stageHexEvaluationTables(*context):stageFlatEvaluationTables(*context);
    }
    FlatChargeTables chargeTables(){
        std::lock_guard<std::mutex> lock(mutex);
        if(!context)throw std::invalid_argument("Material context is closed");
        if(isMixedMesh(*context))return stageMixedChargeTables(*context);
        return isHexMesh(*context)?stageHexChargeTables(*context):stageFlatChargeTables(*context);
    }
    std::vector<FlatEvaluationTables> mixedTables(){
        std::lock_guard<std::mutex> lock(mutex);
        if(!context)throw std::invalid_argument("Material context is closed");
        return stageMixedEvaluationTables(*context);
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
template<class T> py::array_t<T> ownedRows(const std::vector<T>& values,std::size_t rows){
    auto result=ownedArray(values);result.resize({py::ssize_t(rows),py::ssize_t(rows?values.size()/rows:0)});return result;
}
py::dict publishFlat(const FlatEvaluationTables& t){
    py::dict out;
    out["schema"]="radia.hdiv-flat-evaluation-tables.v1";
    out["family"]=t.nodes==27?"hex":t.nodes==10?"curved-tet":"affine-tet";
    out["dofs"]=t.dofs;out["total_points"]=t.totalPoints;out["integration_order"]=t.integrationOrder;
    out["rule"]=ownedRows(t.rule,t.rulePoints);
    auto shapes=ownedArray(t.shapes);
    shapes.resize({py::ssize_t(t.shapes.size()/(3*t.local*t.rulePoints)),py::ssize_t(t.rulePoints),py::ssize_t(t.local),py::ssize_t(3)});
    out["reference_shapes"]=shapes;
    out["vertices"]=ownedRows(t.vertices,t.elements);out["element_class"]=ownedArray(t.elementClass);
    out["element_dofs"]=ownedRows(t.elementDofs,t.elements);out["element_points"]=ownedRows(t.elementPoints,t.elements);
    return out;
}
py::dict publishCharge(const FlatChargeTables& t){
    py::dict out;
    out["schema"]="radia.hdiv-flat-charge-tables.v1";
    out["family"]=t.mixed?"mixed":t.hex?"hex":t.curved?"curved-tet":"affine-tet";
    out["dofs"]=t.dofs;
    if(t.mixed){
        out["host_types"]=ownedArray(t.hostTypes);out["host_class"]=ownedArray(t.hostClass);out["class_info"]=ownedRows(t.classInfo,t.classInfo.size()/3);
        out["host_dofs"]=ownedArray(t.hostDofs);out["class_tables"]=ownedArray(t.classTables);
        out["charge_host"]=ownedArray(t.chargeHost);out["charge_kind"]=ownedArray(t.chargeKind);out["charge_exponents"]=ownedArray(t.chargeExponents);
        out["cell_vertices"]=ownedRows(t.cellVertices,t.cells);out["face_vertices"]=ownedRows(t.faceVertices,t.faces);
        return out;
    }
    auto tables=[](const std::vector<double>& v,std::size_t monomials,std::size_t width){auto a=ownedArray(v);
        a.resize({py::ssize_t(monomials*width?v.size()/(monomials*width):0),py::ssize_t(monomials),py::ssize_t(width)});return a;};
    out["cell_tables"]=tables(t.cellTables,t.cellMonomials,t.cellWidth);out["face_tables"]=tables(t.faceTables,t.faceMonomials,t.faceWidth);
    out["cell_class"]=ownedArray(t.cellClass);out["cell_dofs"]=ownedRows(t.cellDofs,t.cells);
    out["face_class"]=ownedArray(t.faceClass);out["face_dofs"]=ownedRows(t.faceDofs,t.faces);
    out["cell_exponents"]=ownedRows(t.cellExponents,t.cellMonomials);out["face_exponents"]=ownedRows(t.faceExponents,t.faceMonomials);
    out["cell_vertices"]=ownedRows(t.cellVertices,t.cells);out["face_vertices"]=ownedRows(t.faceVertices,t.faces);
    if(t.curved){out["cell_numbers"]=ownedRows(t.cellNumbers,t.cells);out["face_numbers"]=ownedRows(t.faceNumbers,t.faces);
        out["cell_mass_rule"]=ownedRows(t.cellRule,t.cellRule.size()/4);out["face_mass_rule"]=ownedRows(t.faceRule,t.faceRule.size()/4);}
    return out;
}
}
void registerMaterialStaging(py::module_& module){
    py::class_<MaterialOwner>(module,"_MaterialStagingContext")
        .def(py::init<const std::string&,int,int,const std::vector<std::string>&>(),py::arg("mesh_path"),py::arg("order"),py::arg("curve_order")=1,
             py::arg("periodic_boundaries")=std::vector<std::string>{},py::call_guard<py::gil_scoped_release>())
        .def("quadrature",[](MaterialOwner& self,std::size_t first,std::size_t count){
            MaterialQuadraturePayload payload;
            {py::gil_scoped_release release;payload=self.quadrature(first,count);}
            return publish(payload);
        },py::arg("first_element")=0,py::arg("element_count")=std::numeric_limits<std::size_t>::max())
        .def("flat_tables",[](MaterialOwner& self){
            FlatEvaluationTables tables;
            {py::gil_scoped_release release;tables=self.flatTables();}
            return publishFlat(tables);
        },"NGSolve reference HDiv shapes per vertex-rank class plus primitive geometry (affine TET vertices, curved P2 TET nodes or HEX Q2 lattices); no evaluation matrix.")
        .def("mixed_tables",[](MaterialOwner& self){
            std::vector<FlatEvaluationTables> groups;
            {py::gil_scoped_release release;groups=self.mixedTables();}
            py::list out;for(const auto& t:groups){auto d=publishFlat(t);d["family"]="mixed";out.append(d);}
            return out;
        },"One flat_tables payload per element type of a straight tet, pyramid and hex mesh, with global DOF and point numbers")
        .def("charge_tables",[](MaterialOwner& self){
            FlatChargeTables tables;
            {py::gil_scoped_release release;tables=self.chargeTables();}
            return publishCharge(tables);
        },"NGSolve charge tables per vertex-rank class (affine TET projections, curved TET moments or HEX reference projections) plus cell/exterior-face DOFs and geometry; no charge map.")
        .def("close",&MaterialOwner::close,py::call_guard<py::gil_scoped_release>());
}
