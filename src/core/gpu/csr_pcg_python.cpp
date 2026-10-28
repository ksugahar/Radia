#ifndef RADIA_CSR_PYTHON_MODULE
#define RADIA_CSR_PYTHON_MODULE _gpu_csr
#endif
#include "csr_pcg_host.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py=pybind11;
template<class T> using Array=py::array_t<T,py::array::c_style>;
static py::tuple solve(Array<int> row,Array<int> col,Array<double> val,
                       Array<double> rhs,Array<double> diagonal,double tol,int maxiter){
    if(row.ndim()!=1 || col.ndim()!=1 || val.ndim()!=1 || rhs.ndim()!=1 || diagonal.ndim()!=1)
        throw py::value_error("CSR inputs must be one-dimensional C-contiguous arrays");
    const size_t n=rhs.size(),nnz=val.size();
    radia_csr::validate_dimensions(n,nnz,tol,maxiter);
    if(row.size()!=n+1 || col.size()!=nnz || diagonal.size()!=n)
        throw py::value_error("CSR array lengths disagree");
    // Own snapshots while the GIL is released; another Python thread cannot
    // mutate the buffers used for validation or asynchronous device transfers.
    std::vector<int> rp(row.data(),row.data()+row.size()),ci(col.data(),col.data()+col.size());
    std::vector<double> av(val.data(),val.data()+val.size()),b(rhs.data(),rhs.data()+rhs.size()),di(diagonal.data(),diagonal.data()+diagonal.size());
    py::array_t<double> result(n);
    double* output=result.mutable_data();
    radia_csr::Result info{};
    {py::gil_scoped_release release;
     info=radia_csr::solve_gpu_host(static_cast<int>(n),static_cast<int>(nnz),rp.data(),ci.data(),av.data(),b.data(),di.data(),tol,maxiter,output);}
    py::dict metadata;metadata["backend"]="gpu";metadata["precision"]="float64";
    metadata["iterations"]=info.iterations;metadata["true_relative_residual"]=info.true_relative_residual;
    metadata["preconditioner"]="Jacobi; not IC or AMS";
    metadata["transfer_included"]=true;
    return py::make_tuple(result,metadata);
}
// Private phase used while assembling the shared AMG binding. It makes no
// convergence promise and must never be substituted for a complete solve.
static py::array_t<double> amg_smooth(Array<int> row,Array<int> col,Array<double> val,
                                     Array<double> rhs,Array<double> l1,int sweeps) {
    if(row.ndim()!=1||col.ndim()!=1||val.ndim()!=1||rhs.ndim()!=1||l1.ndim()!=1)
        throw py::value_error("AMG arrays must be one-dimensional and C-contiguous");
    const size_t n=rhs.size(),nnz=val.size();
    radia_csr::validate_dimensions(n,nnz,1e-6,sweeps);
    if(row.size()!=n+1||col.size()!=nnz||l1.size()!=n)
        throw py::value_error("AMG array lengths disagree");
    std::vector<int> rp(row.data(),row.data()+row.size()),ci(col.data(),col.data()+col.size());
    std::vector<double> av(val.data(),val.data()+val.size()),b(rhs.data(),rhs.data()+rhs.size()),weights(l1.data(),l1.data()+l1.size());
    py::array_t<double> result(n);double* output=result.mutable_data();
    {py::gil_scoped_release release;
     radia_csr::amg_smooth_gpu_host(int(n),int(nnz),rp.data(),ci.data(),av.data(),b.data(),weights.data(),sweeps,output);}
    return result;
}
PYBIND11_MODULE(RADIA_CSR_PYTHON_MODULE,m){
    m.doc()="Explicit CSR CUDA solve; no MATLAB ABI or import-time CUDA calls";
    m.def("_amg_l1_smooth",&amg_smooth,py::arg("row_pointers").noconvert(),py::arg("columns").noconvert(),
          py::arg("values").noconvert(),py::arg("rhs").noconvert(),py::arg("l1_weights").noconvert(),py::arg("sweeps"),
          "Internal shared AMS/AMG smoothing phase; not a complete preconditioner or solver");
    m.def("solve",&solve,py::arg("row_pointers").noconvert(),py::arg("columns").noconvert(),
          py::arg("values").noconvert(),py::arg("rhs").noconvert(),py::arg("inverse_diagonal").noconvert(),
          py::arg("tolerance"),py::arg("max_iterations"));
}
