// HACApK H-matrix of the tetrahedral volume PEEC partial inductance (rad_peec_tet).
//
// One DOF per interior face (the lowest-order H(div) face current); the cluster coordinate is the
// face centroid and the support radius covers both owning tetrahedra, so touching supports stay in
// dense leaves.  Entries come from radia::peec_tet::InductanceKernel, the same rules as the dense
// assembly, and the operator is applied with the exactly symmetric H-matrix product.
#pragma once

#include <memory>
#include <vector>

#include "rad_hacapk.h"
#include "rad_peec_tet.h"

class RadHACApKTetPEEC : public RadHACApKBase {
public:
    explicit RadHACApKTetPEEC(std::unique_ptr<radia::peec_tet::InductanceKernel> kernel);
    ~RadHACApKTetPEEC() override = default;

    double GetInteractionMatrixElement(int dof_i, int dof_j) const override;
    const radia::peec_tet::InductanceKernel& Kernel() const { return *m_kernel; }

    // Y = L X for nrhs real right-hand sides stored row-major [nrhs][ndof].
    void Apply(const std::vector<double>& x, int nrhs, std::vector<double>& y);

protected:
    void ExtractCoordinates() override;
    void OnBeforeBuild() override {}
    void InitializeInvChi() override;
    bool IsVariableDOF() const override { return false; }
    int GetUniformNFFC() const override { return 1; }
    bool UseSymmetricFill() const override { return true; }

private:
    std::unique_ptr<radia::peec_tet::InductanceKernel> m_kernel;
};

// Parameters for the 1/R face-current kernel: dense near field, ranks up to 400.
inline RadHACApKParams RadHACApKTetPEECDefaultParams() {
    RadHACApKParams p;
    p.aca_eps = 1.0e-8;
    p.leaf_size = 64;
    p.eta = 2.0;
    p.max_rank = 400;
    p.print_level = 0;
    return p;
}
