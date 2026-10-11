/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include "sparsesolv/preconditioners/compact_ams.hpp"
#include <sparsecholesky.hpp>
#include <limits>

namespace py = pybind11;
namespace ngla {
template<class T, class Range> py::array_t<T> AcceleratorArray(const Range& values) {
    py::array_t<T> out(values.Size());
    auto dst = out.mutable_unchecked<1>();
    for (size_t i = 0; i < size_t(values.Size()); ++i) dst(i) = values[i];
    return out;
}
inline py::dict AcceleratorCSR(const SparseMatrix<double>& a) {
    if (a.NZE() > size_t(std::numeric_limits<int>::max()))
        throw std::overflow_error("AMS CUDA bridge requires int32 CSR capacity");
    py::array_t<int> rows(a.Height()+1), cols(a.NZE());
    py::array_t<double> vals(a.NZE());
    auto r=rows.mutable_unchecked<1>(); auto c=cols.mutable_unchecked<1>();
    auto v=vals.mutable_unchecked<1>();
    int k=0;
    for (int i=0; i<a.Height(); ++i) {
        r(i)=k;
        auto ci=a.GetRowIndices(i); auto av=a.GetRowValues(i);
        for (int j=0; j<ci.Size(); ++j) { c(k)=ci[j]; v(k++)=av[j]; }
    }
    r(a.Height())=k;
    py::dict out; out["rows"]=a.Height(); out["cols"]=a.Width();
    out["rowptr"]=rows; out["columns"]=cols; out["values"]=vals;
    return out;
}
inline py::dict AcceleratorCoarse(const BaseMatrix* inv,
                                  const std::vector<double>& dense, int n) {
    py::dict out;
    if (!dense.empty()) {
        py::array_t<double> data(dense.size());
        std::copy(dense.begin(),dense.end(),data.mutable_data());
        out["kind"]="dense_inverse"; out["values"]=data;
    } else if (inv) {
        auto factor=dynamic_cast<const SparseCholesky<double>*>(inv);
        if (!factor) throw std::runtime_error("AMS accelerator requires sparsecholesky coarse factors");
        // NGSolve owns factorization, ordering, null rows and compressed storage.
        // Expand only the strictly upper factor's row indices, never its inverse.
        auto first=factor->GetFirstInRow(); auto ri=factor->GetFirstInRowRI();
        auto indices=factor->GetRowIndex2(); auto lf=factor->GetLFact();
        if (lf.Size()>size_t(std::numeric_limits<int>::max()))
            throw std::overflow_error("AMS coarse factor exceeds int32 capacity");
        py::array_t<int> rp(first.Size()), ci(lf.Size());
        for (int i=0; i<first.Size(); ++i) rp.mutable_data()[i]=int(first[i]);
        for (int i=0; i<factor->GetNUsed(); ++i)
            for (size_t j=first[i]; j<first[i+1]; ++j)
                ci.mutable_data()[j]=indices[ri[i]+j-first[i]];
        out["kind"]="ngsolve_sparse_ldlt";
        out["rowptr"]=rp; out["columns"]=ci;
        out["values"]=AcceleratorArray<double>(lf);
        out["order"]=AcceleratorArray<int>(factor->GetOrder());
        out["diagonal_inverse"]=AcceleratorArray<double>(factor->GetDiag());
    } else out["kind"]="ten_l1_jacobi_sweeps";
    out["rows"]=n;
    return out;
}
inline py::list AcceleratorAMG(const CompactAMG& amg) {
    py::list levels;
    amg.VisitAcceleratorLevels([&](size_t index, const auto& a, const auto* p,
        const auto* pt, const auto& l1, const auto* inverse, const auto& dense, int smooth) {
        py::dict l; l["A"]=AcceleratorCSR(a); l["num_smooth"]=smooth;
        l["l1"]=py::array_t<double>(l1.size(),l1.data());
        if (p) l["P"]=AcceleratorCSR(*p);
        if (pt) l["Pt"]=AcceleratorCSR(*pt);
        if (index+1==size_t(amg.NumLevels()))
            l["coarse"]=AcceleratorCoarse(inverse,dense,a.Height());
        levels.append(l);
    });
    return levels;
}
inline py::dict AcceleratorAMS(const HypreBasedAMS& ams) {
    if (ams.GetMixedPrecision()) throw std::runtime_error("AMS CUDA snapshot requires FP64");
    if (ams.GetNumSmooth()<1 || (ams.GetCycleType()!=1 && ams.GetCycleType()!=7))
        throw std::runtime_error("AMS CUDA snapshot requires cycle 1 or 7 and positive smoothing count");
    py::dict out; out["schema"]="radia.ams-accelerator.v1";
    out["A"]=AcceleratorCSR(ams.GetAbc());
    const auto& norms=ams.GetL1Norms();
    out["l1"]=py::array_t<double>(norms.size(),norms.data());
    out["cycle_type"]=ams.GetCycleType(); out["num_smooth"]=ams.GetNumSmooth();
    out["correction_weight"]=ams.GetCorrectionWeight(); out["beta_zero"]=ams.GetBetaZero();
    py::array_t<int> free(ams.GetNdofHC());
    auto freedofs=ams.GetFreeDofs();
    for (int i=0;i<ams.GetNdofHC();++i)
        free.mutable_data()[i]=!freedofs || freedofs->Test(i);
    out["free"]=free;
    py::list spaces;
    auto append=[&](const auto& p,const auto& pt, const CompactAMG* amg) {
        if (!amg) throw std::runtime_error("AMS CUDA snapshot requires CompactAMG auxiliary solvers");
        py::dict s; s["P"]=AcceleratorCSR(p); s["Pt"]=AcceleratorCSR(pt);
        s["levels"]=AcceleratorAMG(*amg); spaces.append(s);
    };
    if (!ams.GetBetaZero()) append(ams.GetGrad(),ams.GetGradT(),ams.GetBGAsAMG());
    append(ams.GetPix(),ams.GetPixT(),ams.GetBPixAsAMG());
    append(ams.GetPiy(),ams.GetPiyT(),ams.GetBPiyAsAMG());
    append(ams.GetPiz(),ams.GetPizT(),ams.GetBPizAsAMG());
    out["spaces"]=spaces;
    return out;
}
}
