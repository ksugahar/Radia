/*
 * rad_mmm_matrices.cpp
 *
 * MMM (Magnetic Moment Method) standalone matrix construction
 *
 * Part of Radia project - pybind11 module: mmm_core
 */

#include "rad_mmm_matrices.h"

#include <cstring>
#include <algorithm>
#include <stdexcept>
#include <cmath>

// Intel MKL for BLAS/LAPACK
#include <mkl.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radia {
namespace mmm {

//=========================================================================
// TriangleFace Implementation
//=========================================================================

void TriangleFace::ComputeDerivedQuantities() {
    edge01 = v1 - v0;
    edge02 = v2 - v0;
    Vec3d cross = edge01.cross(edge02);
    area = cross.norm() * 0.5;
    if (area > 1e-30) {
        normal = cross / (2.0 * area);
    } else {
        normal = Vec3d(0, 0, 1);
    }
}

//=========================================================================
// MMMBuilder Implementation
//=========================================================================

MMMBuilder::MMMBuilder() {}

MMMBuilder::~MMMBuilder() {}

void MMMBuilder::Clear() {
    elements_.clear();
}

int MMMBuilder::TotalDOF() const {
    int total = 0;
    for (const auto& elem : elements_) {
        total += elem.dof;
    }
    return total;
}

void MMMBuilder::AddTetrahedron(const std::vector<Vec3d>& vertices) {
    if (vertices.size() != 4) {
        throw std::runtime_error("Tetrahedron requires exactly 4 vertices");
    }

    MMMElement elem;
    elem.dof = 3;  // MMM: 3 DOF (Mx, My, Mz)

    // Compute centroid
    elem.center = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / 4.0;

    // Compute volume using scalar triple product
    Vec3d a = vertices[1] - vertices[0];
    Vec3d b = vertices[2] - vertices[0];
    Vec3d c = vertices[3] - vertices[0];
    elem.volume = std::abs(a.dot(b.cross(c))) / 6.0;

    // Create 4 triangular faces with outward normals
    // Standard tetrahedron face ordering (opposite vertex indices)
    // Face 0: v1, v3, v2 (opposite v0)
    // Face 1: v0, v2, v3 (opposite v1)
    // Face 2: v0, v3, v1 (opposite v2)
    // Face 3: v0, v1, v2 (opposite v3)

    int face_indices[4][3] = {
        {1, 3, 2},  // Face opposite v0
        {0, 2, 3},  // Face opposite v1
        {0, 3, 1},  // Face opposite v2
        {0, 1, 2}   // Face opposite v3
    };

    elem.triangles.resize(4);
    for (int f = 0; f < 4; f++) {
        TriangleFace& tri = elem.triangles[f];
        tri.v0 = vertices[face_indices[f][0]];
        tri.v1 = vertices[face_indices[f][1]];
        tri.v2 = vertices[face_indices[f][2]];
        tri.ComputeDerivedQuantities();

        // Ensure normal points outward (away from centroid)
        Vec3d to_center = elem.center - tri.v0;
        if (tri.normal.dot(to_center) > 0) {
            // Normal points inward, flip it
            std::swap(tri.v1, tri.v2);
            tri.ComputeDerivedQuantities();
        }
    }

    elements_.push_back(std::move(elem));
}

void MMMBuilder::AddTetrahedron(const double* vertices) {
    std::vector<Vec3d> verts(4);
    for (int i = 0; i < 4; i++) {
        verts[i] = Vec3d(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
    }
    AddTetrahedron(verts);
}

void MMMBuilder::AddHexahedron(const std::vector<Vec3d>& vertices) {
    if (vertices.size() != 8) {
        throw std::runtime_error("Hexahedron requires exactly 8 vertices");
    }

    MMMElement elem;
    elem.dof = 6;  // MSC: 6 DOF (sigma per face)

    // Compute centroid (average of 8 vertices)
    elem.center = Vec3d(0, 0, 0);
    for (const auto& v : vertices) {
        elem.center = elem.center + v;
    }
    elem.center = elem.center / 8.0;

    // Approximate volume (assumes nearly rectangular hex)
    // For accurate volume, we'd need to split into tetrahedra
    // This is sufficient for interaction matrix computation
    Vec3d diag = vertices[7] - vertices[0];
    Vec3d dx = vertices[1] - vertices[0];
    Vec3d dy = vertices[3] - vertices[0];
    Vec3d dz = vertices[4] - vertices[0];
    elem.volume = std::abs(dx.dot(dy.cross(dz)));

    // Standard hexahedron face ordering:
    // Face 0: bottom (v0, v3, v2, v1) - normal -Z
    // Face 1: top    (v4, v5, v6, v7) - normal +Z
    // Face 2: front  (v0, v1, v5, v4) - normal -Y
    // Face 3: back   (v2, v3, v7, v6) - normal +Y
    // Face 4: left   (v0, v4, v7, v3) - normal -X
    // Face 5: right  (v1, v2, v6, v5) - normal +X

    int face_indices[6][4] = {
        {0, 3, 2, 1},  // Bottom
        {4, 5, 6, 7},  // Top
        {0, 1, 5, 4},  // Front
        {2, 3, 7, 6},  // Back
        {0, 4, 7, 3},  // Left
        {1, 2, 6, 5}   // Right
    };

    // Each quad face -> 2 triangles
    elem.triangles.resize(12);
    elem.face_normals.resize(6);
    elem.face_areas.resize(6);
    elem.eval_points.resize(6);

    for (int f = 0; f < 6; f++) {
        Vec3d v0 = vertices[face_indices[f][0]];
        Vec3d v1 = vertices[face_indices[f][1]];
        Vec3d v2 = vertices[face_indices[f][2]];
        Vec3d v3 = vertices[face_indices[f][3]];

        // Triangle 1: v0, v1, v2
        TriangleFace& tri1 = elem.triangles[f * 2];
        tri1.v0 = v0;
        tri1.v1 = v1;
        tri1.v2 = v2;
        tri1.ComputeDerivedQuantities();

        // Triangle 2: v0, v2, v3
        TriangleFace& tri2 = elem.triangles[f * 2 + 1];
        tri2.v0 = v0;
        tri2.v1 = v2;
        tri2.v2 = v3;
        tri2.ComputeDerivedQuantities();

        // Face normal (average of two triangles)
        elem.face_normals[f] = (tri1.normal + tri2.normal).normalized();

        // Total face area
        elem.face_areas[f] = tri1.area + tri2.area;

        // Evaluation point (Yano-Sugahara: face centroid)
        elem.eval_points[f] = (v0 + v1 + v2 + v3) / 4.0;

        // Ensure normals point outward
        Vec3d to_center = elem.center - elem.eval_points[f];
        if (elem.face_normals[f].dot(to_center) > 0) {
            elem.face_normals[f] = elem.face_normals[f] * (-1.0);
            // Also flip triangle normals
            tri1.normal = tri1.normal * (-1.0);
            tri2.normal = tri2.normal * (-1.0);
        }
    }

    elements_.push_back(std::move(elem));
}

void MMMBuilder::AddHexahedron(const double* vertices) {
    std::vector<Vec3d> verts(8);
    for (int i = 0; i < 8; i++) {
        verts[i] = Vec3d(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
    }
    AddHexahedron(verts);
}

//=========================================================================
// Permanent Magnet Element Methods
//=========================================================================

void MMMBuilder::AddPermanentMagnetTetrahedron(const std::vector<Vec3d>& vertices,
                                                 const Vec3d& Mr) {
    if (vertices.size() != 4) {
        throw std::runtime_error("Tetrahedron requires exactly 4 vertices");
    }

    MMMElement elem;
    elem.dof = 0;  // PM elements have 0 DOF in solve (fixed magnetization)
    elem.material_type = MMMMaterialType::FixedPM;
    elem.Mr = Mr;

    // Compute centroid
    elem.center = (vertices[0] + vertices[1] + vertices[2] + vertices[3]) / 4.0;

    // Compute volume
    Vec3d a = vertices[1] - vertices[0];
    Vec3d b = vertices[2] - vertices[0];
    Vec3d c = vertices[3] - vertices[0];
    elem.volume = std::abs(a.dot(b.cross(c))) / 6.0;

    // Create 4 triangular faces (same as soft tetrahedron)
    int face_indices[4][3] = {
        {1, 3, 2}, {0, 2, 3}, {0, 3, 1}, {0, 1, 2}
    };

    elem.triangles.resize(4);
    for (int f = 0; f < 4; f++) {
        TriangleFace& tri = elem.triangles[f];
        tri.v0 = vertices[face_indices[f][0]];
        tri.v1 = vertices[face_indices[f][1]];
        tri.v2 = vertices[face_indices[f][2]];
        tri.ComputeDerivedQuantities();

        // Ensure normal points outward
        Vec3d to_center = elem.center - tri.v0;
        if (tri.normal.dot(to_center) > 0) {
            std::swap(tri.v1, tri.v2);
            tri.ComputeDerivedQuantities();
        }
    }

    elements_.push_back(std::move(elem));
}

void MMMBuilder::AddPermanentMagnetTetrahedron(const double* vertices, const double* Mr) {
    std::vector<Vec3d> verts(4);
    for (int i = 0; i < 4; i++) {
        verts[i] = Vec3d(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
    }
    Vec3d mr(Mr[0], Mr[1], Mr[2]);
    AddPermanentMagnetTetrahedron(verts, mr);
}

void MMMBuilder::AddPermanentMagnetHexahedron(const std::vector<Vec3d>& vertices,
                                                const Vec3d& Mr) {
    if (vertices.size() != 8) {
        throw std::runtime_error("Hexahedron requires exactly 8 vertices");
    }

    MMMElement elem;
    elem.dof = 0;  // PM elements have 0 DOF in solve
    elem.material_type = MMMMaterialType::FixedPM;
    elem.Mr = Mr;

    // Compute centroid
    elem.center = Vec3d(0, 0, 0);
    for (const auto& v : vertices) {
        elem.center = elem.center + v;
    }
    elem.center = elem.center / 8.0;

    // Compute volume
    Vec3d dx = vertices[1] - vertices[0];
    Vec3d dy = vertices[3] - vertices[0];
    Vec3d dz = vertices[4] - vertices[0];
    elem.volume = std::abs(dx.dot(dy.cross(dz)));

    // Create 6 quad faces (same as soft hexahedron)
    int face_indices[6][4] = {
        {0, 3, 2, 1}, {4, 5, 6, 7}, {0, 1, 5, 4},
        {2, 3, 7, 6}, {0, 4, 7, 3}, {1, 2, 6, 5}
    };

    elem.triangles.resize(12);
    elem.face_normals.resize(6);
    elem.face_areas.resize(6);
    elem.eval_points.resize(6);

    for (int f = 0; f < 6; f++) {
        Vec3d v0 = vertices[face_indices[f][0]];
        Vec3d v1 = vertices[face_indices[f][1]];
        Vec3d v2 = vertices[face_indices[f][2]];
        Vec3d v3 = vertices[face_indices[f][3]];

        TriangleFace& tri1 = elem.triangles[f * 2];
        tri1.v0 = v0; tri1.v1 = v1; tri1.v2 = v2;
        tri1.ComputeDerivedQuantities();

        TriangleFace& tri2 = elem.triangles[f * 2 + 1];
        tri2.v0 = v0; tri2.v1 = v2; tri2.v2 = v3;
        tri2.ComputeDerivedQuantities();

        elem.face_normals[f] = (tri1.normal + tri2.normal).normalized();
        elem.face_areas[f] = tri1.area + tri2.area;
        elem.eval_points[f] = (v0 + v1 + v2 + v3) / 4.0;

        Vec3d to_center = elem.center - elem.eval_points[f];
        if (elem.face_normals[f].dot(to_center) > 0) {
            elem.face_normals[f] = elem.face_normals[f] * (-1.0);
            tri1.normal = tri1.normal * (-1.0);
            tri2.normal = tri2.normal * (-1.0);
        }
    }

    elements_.push_back(std::move(elem));
}

void MMMBuilder::AddPermanentMagnetHexahedron(const double* vertices, const double* Mr) {
    std::vector<Vec3d> verts(8);
    for (int i = 0; i < 8; i++) {
        verts[i] = Vec3d(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
    }
    Vec3d mr(Mr[0], Mr[1], Mr[2]);
    AddPermanentMagnetHexahedron(verts, mr);
}

double MMMBuilder::SolidAngleTriangle(const Vec3d& obs, const Vec3d& v0,
                                        const Vec3d& v1, const Vec3d& v2) const {
    // van Oosterom & Strackee formula for solid angle of a triangle
    Vec3d r0 = v0 - obs;
    Vec3d r1 = v1 - obs;
    Vec3d r2 = v2 - obs;

    double d0 = r0.norm();
    double d1 = r1.norm();
    double d2 = r2.norm();

    // Handle degenerate cases
    if (d0 < 1e-30 || d1 < 1e-30 || d2 < 1e-30) {
        return 0.0;  // Point on vertex
    }

    double numer = r0.dot(r1.cross(r2));
    double denom = d0 * d1 * d2
                 + d0 * r1.dot(r2)
                 + d1 * r0.dot(r2)
                 + d2 * r0.dot(r1);

    return 2.0 * std::atan2(numer, denom);
}

Vec3d MMMBuilder::FieldFromChargedTriangle(const TriangleFace& tri, const Vec3d& obs,
                                            double sigma) const {
    // H field from uniformly charged triangle using solid angle
    // H = -sigma / (4*pi) * grad(Omega)
    // For uniform charge: H = sigma * Omega / (4*pi) * n (pointing toward surface)

    // Compute solid angle
    double omega = SolidAngleTriangle(obs, tri.v0, tri.v1, tri.v2);

    // H field contribution (using solid angle gradient approach)
    // This is the simplified formula for field from a uniformly charged triangle
    Vec3d H = tri.normal * (sigma * omega * constants::INV_FOUR_PI);

    return H;
}

void MMMBuilder::Compute3x3Block(int elem_i, int elem_j, double* N_block) const {
    // Compute 3x3 interaction block between two tetrahedra
    // N[i,j] relates magnetization M_j to field at element i
    // H_i = N[i,j] * M_j

    const MMMElement& elem_src = elements_[elem_j];  // Source element
    const MMMElement& elem_obs = elements_[elem_i];  // Observer element

    // Initialize block to zero
    std::memset(N_block, 0, 9 * sizeof(double));

    // For tetrahedra: use surface charge method
    // Each face has charge sigma = M . n
    // Field at observer = sum over faces of field from charged triangle

    Vec3d obs = elem_obs.center;

    // For each face of source element
    for (const auto& tri : elem_src.triangles) {
        // Field contribution when M = ex (unit x-magnetization)
        // sigma = M . n = n_x
        Vec3d H_x = FieldFromChargedTriangle(tri, obs, tri.normal.x);

        // Field contribution when M = ey
        Vec3d H_y = FieldFromChargedTriangle(tri, obs, tri.normal.y);

        // Field contribution when M = ez
        Vec3d H_z = FieldFromChargedTriangle(tri, obs, tri.normal.z);

        // Accumulate into N matrix (row-major)
        // N[0,:] = dH_x/dM (row 0)
        N_block[0] += H_x.x;  N_block[1] += H_y.x;  N_block[2] += H_z.x;
        // N[1,:] = dH_y/dM (row 1)
        N_block[3] += H_x.y;  N_block[4] += H_y.y;  N_block[5] += H_z.y;
        // N[2,:] = dH_z/dM (row 2)
        N_block[6] += H_x.z;  N_block[7] += H_y.z;  N_block[8] += H_z.z;
    }
}

void MMMBuilder::Compute6x6Block(int elem_i, int elem_j, double* K_block) const {
    // Compute 6x6 MSC interaction block between two hexahedra
    // K[i,j] relates face charges sigma_j to field at evaluation points of element i

    const MMMElement& elem_src = elements_[elem_j];  // Source element
    const MMMElement& elem_obs = elements_[elem_i];  // Observer element

    // Initialize block to zero
    std::memset(K_block, 0, 36 * sizeof(double));

    // For MSC hexahedra: 6 evaluation points, 6 source faces
    // K[f_obs, f_src] = n_obs . H(eval_point_obs, sigma=1 on face f_src)

    for (int f_obs = 0; f_obs < 6; f_obs++) {
        Vec3d obs = elem_obs.eval_points[f_obs];
        Vec3d n_obs = elem_obs.face_normals[f_obs];

        for (int f_src = 0; f_src < 6; f_src++) {
            // Two triangles per face
            Vec3d H_total(0, 0, 0);

            // Triangle 1
            const TriangleFace& tri1 = elem_src.triangles[f_src * 2];
            H_total = H_total + FieldFromChargedTriangle(tri1, obs, 1.0);

            // Triangle 2
            const TriangleFace& tri2 = elem_src.triangles[f_src * 2 + 1];
            H_total = H_total + FieldFromChargedTriangle(tri2, obs, 1.0);

            // K[f_obs, f_src] = n_obs . H_total
            K_block[f_obs * 6 + f_src] = n_obs.dot(H_total);
        }
    }
}

void MMMBuilder::Compute3x6Block(int tetra_idx, int hexa_idx, double* N_block) const {
    // Tetrahedron observer, hexahedron source
    // 3 rows (observer DOF), 6 columns (source faces)

    const MMMElement& elem_src = elements_[hexa_idx];  // Source hexahedron
    const MMMElement& elem_obs = elements_[tetra_idx]; // Observer tetrahedron

    std::memset(N_block, 0, 18 * sizeof(double));

    Vec3d obs = elem_obs.center;

    for (int f_src = 0; f_src < 6; f_src++) {
        Vec3d H_total(0, 0, 0);

        // Two triangles per face
        const TriangleFace& tri1 = elem_src.triangles[f_src * 2];
        H_total = H_total + FieldFromChargedTriangle(tri1, obs, 1.0);

        const TriangleFace& tri2 = elem_src.triangles[f_src * 2 + 1];
        H_total = H_total + FieldFromChargedTriangle(tri2, obs, 1.0);

        // N[Hx, sigma_f] = H_total.x, etc.
        N_block[0 * 6 + f_src] = H_total.x;
        N_block[1 * 6 + f_src] = H_total.y;
        N_block[2 * 6 + f_src] = H_total.z;
    }
}

void MMMBuilder::Compute6x3Block(int hexa_idx, int tetra_idx, double* N_block) const {
    // Hexahedron observer, tetrahedron source
    // 6 rows (observer faces), 3 columns (source magnetization)

    const MMMElement& elem_src = elements_[tetra_idx]; // Source tetrahedron
    const MMMElement& elem_obs = elements_[hexa_idx];  // Observer hexahedron

    std::memset(N_block, 0, 18 * sizeof(double));

    for (int f_obs = 0; f_obs < 6; f_obs++) {
        Vec3d obs = elem_obs.eval_points[f_obs];
        Vec3d n_obs = elem_obs.face_normals[f_obs];

        // Compute H at evaluation point from tetra with unit magnetization in each direction
        for (const auto& tri : elem_src.triangles) {
            Vec3d H_x = FieldFromChargedTriangle(tri, obs, tri.normal.x);
            Vec3d H_y = FieldFromChargedTriangle(tri, obs, tri.normal.y);
            Vec3d H_z = FieldFromChargedTriangle(tri, obs, tri.normal.z);

            // N[f_obs, Mx] += n_obs . H_x, etc.
            N_block[f_obs * 3 + 0] += n_obs.dot(H_x);
            N_block[f_obs * 3 + 1] += n_obs.dot(H_y);
            N_block[f_obs * 3 + 2] += n_obs.dot(H_z);
        }
    }
}

MMMMatrices MMMBuilder::Build() {
    MMMMatrices result;
    int n_elem = static_cast<int>(elements_.size());

    if (n_elem == 0) {
        return result;
    }

    // Compute DOF offsets and count element types
    result.n_elements = n_elem;
    result.dof_offset.resize(n_elem + 1);
    result.material_types.resize(n_elem);
    result.Mr.resize(n_elem * 3, 0.0);
    result.n_pm_elements = 0;
    result.n_soft_elements = 0;

    result.dof_offset[0] = 0;
    for (int i = 0; i < n_elem; i++) {
        result.dof_offset[i + 1] = result.dof_offset[i] + elements_[i].dof;
        result.material_types[i] = elements_[i].material_type;

        // Store Mr for all elements (0 for soft materials)
        result.Mr[i * 3 + 0] = elements_[i].Mr.x;
        result.Mr[i * 3 + 1] = elements_[i].Mr.y;
        result.Mr[i * 3 + 2] = elements_[i].Mr.z;

        if (elements_[i].IsPermanentMagnet()) {
            result.n_pm_elements++;
        } else {
            result.n_soft_elements++;
        }
    }
    result.total_dof = result.dof_offset[n_elem];

    // Allocate interaction matrix (column-major for LAPACK)
    size_t matrix_size = static_cast<size_t>(result.total_dof) * result.total_dof;
    result.N.resize(matrix_size, 0.0);

    // Build matrix blocks
    #pragma omp parallel for schedule(dynamic) if(n_elem > 10)
    for (int col = 0; col < n_elem; col++) {
        int dof_col = elements_[col].dof;
        int offset_col = result.dof_offset[col];

        for (int row = 0; row < n_elem; row++) {
            int dof_row = elements_[row].dof;
            int offset_row = result.dof_offset[row];

            // Allocate temporary block
            std::vector<double> block(dof_row * dof_col, 0.0);

            // Compute block based on element types
            if (dof_row == 3 && dof_col == 3) {
                Compute3x3Block(row, col, block.data());
            } else if (dof_row == 6 && dof_col == 6) {
                Compute6x6Block(row, col, block.data());
            } else if (dof_row == 3 && dof_col == 6) {
                Compute3x6Block(row, col, block.data());
            } else if (dof_row == 6 && dof_col == 3) {
                Compute6x3Block(row, col, block.data());
            }

            // Copy block to matrix (transpose for column-major storage)
            for (int j = 0; j < dof_col; j++) {
                for (int i = 0; i < dof_row; i++) {
                    // block is row-major: block[i][j] at i*dof_col+j
                    // result.N is column-major: N[row+i, col+j] at (col+j)*total_dof + (row+i)
                    size_t idx = static_cast<size_t>(offset_col + j) * result.total_dof + (offset_row + i);
                    result.N[idx] = block[i * dof_col + j];
                }
            }
        }
    }

    return result;
}

//=========================================================================
// MMMSolver Implementation
//=========================================================================

MMMSolver::MMMSolver()
    : total_dof_(0), n_elements_(0), matrix_set_(false) {}

MMMSolver::~MMMSolver() {}

void MMMSolver::SetMatrix(const MMMMatrices& matrices) {
    N_ = matrices.N;
    dof_offset_ = matrices.dof_offset;
    total_dof_ = matrices.total_dof;
    n_elements_ = matrices.n_elements;
    matrix_set_ = true;
}

void MMMSolver::SetMatrix(const double* N, const int* dof_offset, int total_dof, int n_elem) {
    size_t matrix_size = static_cast<size_t>(total_dof) * total_dof;
    N_.assign(N, N + matrix_size);
    dof_offset_.assign(dof_offset, dof_offset + n_elem + 1);
    total_dof_ = total_dof;
    n_elements_ = n_elem;
    matrix_set_ = true;
}

void MMMSolver::BuildSystemMatrix(std::vector<double>& A,
                                   const std::vector<double>& inv_chi,
                                   bool chi_per_element) const {
    // Build A = -diag(inv_chi) - N (ELF-compatible)
    size_t n = total_dof_;
    A = N_;  // Copy N

    // Negate: A = -N
    for (size_t i = 0; i < A.size(); i++) {
        A[i] = -A[i];
    }

    // Subtract diagonal: A -= diag(inv_chi) (ELF-compatible: -1/chi)
    if (chi_per_element) {
        // Broadcast chi from elements to DOF
        for (int e = 0; e < n_elements_; e++) {
            double inv_chi_e = inv_chi[e];
            int start = dof_offset_[e];
            int end = dof_offset_[e + 1];
            for (int i = start; i < end; i++) {
                // Column-major: A[i,i] at i*n + i
                A[static_cast<size_t>(i) * n + i] -= inv_chi_e;  // ELF-compatible
            }
        }
    } else {
        // inv_chi is already per-DOF
        for (size_t i = 0; i < n; i++) {
            A[i * n + i] -= inv_chi[i];  // ELF-compatible
        }
    }
}

void MMMSolver::MatVec(const std::vector<double>& A,
                        const std::vector<double>& x,
                        std::vector<double>& y) const {
    // y = A * x using BLAS dgemv
    cblas_dgemv(CblasColMajor, CblasNoTrans,
                total_dof_, total_dof_,
                1.0, A.data(), total_dof_,
                x.data(), 1,
                0.0, y.data(), 1);
}

int MMMSolver::SolveLU(const std::vector<double>& inv_chi,
                        const std::vector<double>& H_ext,
                        std::vector<double>& M,
                        bool chi_per_element) {
    if (!matrix_set_) {
        return -1;
    }

    // Build system matrix A = diag(inv_chi) - N
    std::vector<double> A;
    BuildSystemMatrix(A, inv_chi, chi_per_element);

    // Build RHS: b = H_ext (copy)
    M = H_ext;

    // Solve A * M = H_ext using LAPACK dgesv
    std::vector<int> ipiv(total_dof_);
    int info = LAPACKE_dgesv(LAPACK_COL_MAJOR, total_dof_, 1,
                              A.data(), total_dof_,
                              ipiv.data(),
                              M.data(), total_dof_);

    return info;
}

int MMMSolver::SolveBiCGSTAB(const std::vector<double>& inv_chi,
                              const std::vector<double>& H_ext,
                              std::vector<double>& M,
                              double tol,
                              int max_iter,
                              bool chi_per_element) {
    if (!matrix_set_) {
        return -1;
    }

    int n = total_dof_;

    // Build system matrix A = diag(inv_chi) - N
    std::vector<double> A;
    BuildSystemMatrix(A, inv_chi, chi_per_element);

    // Initialize M = 0
    M.assign(n, 0.0);

    // BiCGSTAB algorithm
    std::vector<double> r(n), r0(n), p(n), v(n), s(n), t(n);

    // r = b - A*x = H_ext - 0 = H_ext
    r = H_ext;
    r0 = r;

    double rho_old = 1.0, alpha = 1.0, omega = 1.0;
    double rho_new, beta;

    double b_norm = cblas_dnrm2(n, H_ext.data(), 1);
    if (b_norm < 1e-30) b_norm = 1.0;

    for (int iter = 0; iter < max_iter; iter++) {
        rho_new = cblas_ddot(n, r0.data(), 1, r.data(), 1);

        if (std::abs(rho_new) < 1e-30) {
            return iter;  // Breakdown
        }

        if (iter == 0) {
            p = r;
        } else {
            beta = (rho_new / rho_old) * (alpha / omega);
            // p = r + beta * (p - omega * v)
            for (int i = 0; i < n; i++) {
                p[i] = r[i] + beta * (p[i] - omega * v[i]);
            }
        }

        // v = A * p
        MatVec(A, p, v);

        double r0v = cblas_ddot(n, r0.data(), 1, v.data(), 1);
        if (std::abs(r0v) < 1e-30) {
            return iter;  // Breakdown
        }

        alpha = rho_new / r0v;

        // s = r - alpha * v
        for (int i = 0; i < n; i++) {
            s[i] = r[i] - alpha * v[i];
        }

        double s_norm = cblas_dnrm2(n, s.data(), 1);
        if (s_norm / b_norm < tol) {
            // x = x + alpha * p
            cblas_daxpy(n, alpha, p.data(), 1, M.data(), 1);
            return iter + 1;
        }

        // t = A * s
        MatVec(A, s, t);

        double tt = cblas_ddot(n, t.data(), 1, t.data(), 1);
        double ts = cblas_ddot(n, t.data(), 1, s.data(), 1);

        if (std::abs(tt) < 1e-30) {
            omega = 0.0;
        } else {
            omega = ts / tt;
        }

        // x = x + alpha * p + omega * s
        cblas_daxpy(n, alpha, p.data(), 1, M.data(), 1);
        cblas_daxpy(n, omega, s.data(), 1, M.data(), 1);

        // r = s - omega * t
        for (int i = 0; i < n; i++) {
            r[i] = s[i] - omega * t[i];
        }

        double r_norm = cblas_dnrm2(n, r.data(), 1);
        if (r_norm / b_norm < tol) {
            return iter + 1;
        }

        rho_old = rho_new;
    }

    return max_iter;  // Did not converge
}

//=========================================================================
// Nonlinear Iteration Helpers
//=========================================================================

std::vector<double> ComputeChiFromBH(
    const std::vector<double>& M_flat,
    const std::vector<double>& H_flat,
    const std::vector<std::array<double, 2>>& bh_curve,
    int n_elem) {

    std::vector<double> chi(n_elem);

    for (int e = 0; e < n_elem; e++) {
        // Compute |H| and |M| for this element
        int offset = e * 3;
        double Hx = H_flat[offset], Hy = H_flat[offset + 1], Hz = H_flat[offset + 2];
        double Mx = M_flat[offset], My = M_flat[offset + 1], Mz = M_flat[offset + 2];

        double H_mag = std::sqrt(Hx*Hx + Hy*Hy + Hz*Hz);
        double M_mag = std::sqrt(Mx*Mx + My*My + Mz*Mz);

        if (H_mag < 1e-10) {
            // Use initial slope
            if (!bh_curve.empty() && bh_curve.size() > 1) {
                double H1 = bh_curve[1][0];
                double B1 = bh_curve[1][1];
                if (H1 > 1e-10) {
                    chi[e] = (B1 / constants::MU_0 - H1) / H1;
                } else {
                    chi[e] = 999.0;  // High permeability default
                }
            } else {
                chi[e] = 999.0;
            }
            continue;
        }

        // Interpolate B-H curve
        double B_interp = 0.0;
        int n_pts = static_cast<int>(bh_curve.size());

        if (H_mag <= bh_curve[0][0]) {
            B_interp = bh_curve[0][1];
        } else if (H_mag >= bh_curve[n_pts - 1][0]) {
            B_interp = bh_curve[n_pts - 1][1];
        } else {
            // Linear interpolation
            for (int i = 1; i < n_pts; i++) {
                if (H_mag <= bh_curve[i][0]) {
                    double H0 = bh_curve[i - 1][0];
                    double H1 = bh_curve[i][0];
                    double B0 = bh_curve[i - 1][1];
                    double B1 = bh_curve[i][1];
                    double t = (H_mag - H0) / (H1 - H0 + 1e-30);
                    B_interp = B0 + t * (B1 - B0);
                    break;
                }
            }
        }

        // chi = M / H = (B/mu_0 - H) / H
        double M_from_B = B_interp / constants::MU_0 - H_mag;
        chi[e] = M_from_B / (H_mag + 1e-30);

        // Clamp to reasonable range
        if (chi[e] < 0.0) chi[e] = 0.0;
        if (chi[e] > 1e6) chi[e] = 1e6;
    }

    return chi;
}

double CheckConvergence(
    const std::vector<double>& B_new,
    const std::vector<double>& B_old,
    double B_sat,
    int n_elem) {

    double max_rel_change = 0.0;

    for (int e = 0; e < n_elem; e++) {
        int offset = e * 3;

        double dBx = B_new[offset] - B_old[offset];
        double dBy = B_new[offset + 1] - B_old[offset + 1];
        double dBz = B_new[offset + 2] - B_old[offset + 2];

        double dB_mag = std::sqrt(dBx*dBx + dBy*dBy + dBz*dBz);
        double rel_change = dB_mag / (B_sat + 1e-30);

        if (rel_change > max_rel_change) {
            max_rel_change = rel_change;
        }
    }

    return max_rel_change;
}

//=========================================================================
// MMMFieldComputer Implementation
//=========================================================================

MMMFieldComputer::MMMFieldComputer() {}

MMMFieldComputer::~MMMFieldComputer() {}

void MMMFieldComputer::SetElements(const std::vector<MMMElement>& elements) {
    elements_ = elements;

    // Build DOF offset array
    int n_elem = static_cast<int>(elements_.size());
    dof_offset_.resize(n_elem + 1);
    dof_offset_[0] = 0;
    for (int i = 0; i < n_elem; i++) {
        dof_offset_[i + 1] = dof_offset_[i] + elements_[i].dof;
    }

    // Initialize material data from elements
    material_types_.resize(n_elem);
    Mr_.resize(n_elem * 3, 0.0);
    for (int i = 0; i < n_elem; i++) {
        material_types_[i] = elements_[i].material_type;
        Mr_[i * 3 + 0] = elements_[i].Mr.x;
        Mr_[i * 3 + 1] = elements_[i].Mr.y;
        Mr_[i * 3 + 2] = elements_[i].Mr.z;
    }
}

void MMMFieldComputer::SetMaterialData(const std::vector<MMMMaterialType>& material_types,
                                        const std::vector<double>& Mr) {
    material_types_ = material_types;
    Mr_ = Mr;
}

int MMMFieldComputer::TotalDOF() const {
    if (dof_offset_.empty()) return 0;
    return dof_offset_.back();
}

std::vector<double> MMMFieldComputer::GetElementCenters() const {
    std::vector<double> centers(elements_.size() * 3);
    for (size_t i = 0; i < elements_.size(); i++) {
        centers[i * 3 + 0] = elements_[i].center.x;
        centers[i * 3 + 1] = elements_[i].center.y;
        centers[i * 3 + 2] = elements_[i].center.z;
    }
    return centers;
}

std::vector<double> MMMFieldComputer::GetElementVolumes() const {
    std::vector<double> volumes(elements_.size());
    for (size_t i = 0; i < elements_.size(); i++) {
        volumes[i] = elements_[i].volume;
    }
    return volumes;
}

Vec3d MMMFieldComputer::DipoleField(const Vec3d& M, const Vec3d& elem_center,
                                     double volume, const Vec3d& obs) const {
    // Dipole field: H = (1/4pi) * [3(m.r)r/r^5 - m/r^3]
    // where m = M * volume (magnetic moment)
    Vec3d r = obs - elem_center;
    double r_mag = r.norm();

    if (r_mag < 1e-12) {
        return Vec3d(0, 0, 0);  // At source point
    }

    Vec3d m = M * volume;
    double r3 = r_mag * r_mag * r_mag;
    double r5 = r3 * r_mag * r_mag;

    double m_dot_r = m.dot(r);

    Vec3d H = r * (3.0 * m_dot_r / r5) - m * (1.0 / r3);
    H = H * constants::INV_FOUR_PI;

    return H;
}

Vec3d MMMFieldComputer::SurfaceChargeField(const MMMElement& elem,
                                            const double* sigma,
                                            const Vec3d& obs) const {
    // For MSC hexahedron: field from 6 surface charges
    // H = sum_i (sigma_i / 4pi) * Omega_i * n_i
    // where Omega_i is solid angle of face i

    Vec3d H(0, 0, 0);

    // Iterate over triangular faces
    int n_tri = static_cast<int>(elem.triangles.size());
    for (int t = 0; t < n_tri; t++) {
        const TriangleFace& tri = elem.triangles[t];

        // Solid angle of triangle
        Vec3d r0 = tri.v0 - obs;
        Vec3d r1 = tri.v1 - obs;
        Vec3d r2 = tri.v2 - obs;

        double d0 = r0.norm();
        double d1 = r1.norm();
        double d2 = r2.norm();

        if (d0 < 1e-12 || d1 < 1e-12 || d2 < 1e-12) {
            continue;  // Observer on triangle
        }

        // Solid angle using van Oosterom & Strackee formula
        double numer = r0.dot(r1.cross(r2));
        double denom = d0 * d1 * d2 +
                       d0 * r1.dot(r2) +
                       d1 * r0.dot(r2) +
                       d2 * r0.dot(r1);

        double omega = 2.0 * std::atan2(numer, denom);

        // For hexahedra, map triangle to face index
        // Each quad face has 2 triangles: face_idx = t / 2
        int face_idx = t / 2;
        if (face_idx < 6) {
            double sigma_f = sigma[face_idx];
            H = H + tri.normal * (sigma_f * omega * constants::INV_FOUR_PI);
        }
    }

    return H;
}

std::vector<double> MMMFieldComputer::ComputeHField(
    const std::vector<double>& M,
    const std::vector<double>& obs_points,
    int n_points) const {

    std::vector<double> H_field(n_points * 3, 0.0);

    int n_elem = static_cast<int>(elements_.size());
    bool has_material_data = !material_types_.empty();

    #pragma omp parallel for schedule(dynamic) if(n_points > 100)
    for (int p = 0; p < n_points; p++) {
        Vec3d obs(obs_points[p * 3], obs_points[p * 3 + 1], obs_points[p * 3 + 2]);
        Vec3d H_total(0, 0, 0);

        for (int e = 0; e < n_elem; e++) {
            const MMMElement& elem = elements_[e];
            int dof_start = dof_offset_[e];

            // Check if this is a permanent magnet element
            bool is_pm = false;
            if (has_material_data && e < static_cast<int>(material_types_.size())) {
                is_pm = (material_types_[e] == MMMMaterialType::FixedPM ||
                         material_types_[e] == MMMMaterialType::DemagPM);
            }

            if (is_pm) {
                // PM element: use fixed Mr (stored in Mr_ array)
                Vec3d M_pm(Mr_[e * 3], Mr_[e * 3 + 1], Mr_[e * 3 + 2]);
                H_total = H_total + DipoleField(M_pm, elem.center, elem.volume, obs);
            } else if (elem.dof == 3) {
                // Soft tetrahedron: use solved M
                Vec3d M_elem(M[dof_start], M[dof_start + 1], M[dof_start + 2]);
                H_total = H_total + DipoleField(M_elem, elem.center, elem.volume, obs);
            } else if (elem.dof == 6) {
                // Soft hexahedron: surface charge method
                H_total = H_total + SurfaceChargeField(elem, &M[dof_start], obs);
            } else if (elem.dof == 0 && !is_pm) {
                // Element with 0 DOF but not marked as PM - use stored Mr from element
                H_total = H_total + DipoleField(elem.Mr, elem.center, elem.volume, obs);
            }
        }

        H_field[p * 3 + 0] = H_total.x;
        H_field[p * 3 + 1] = H_total.y;
        H_field[p * 3 + 2] = H_total.z;
    }

    return H_field;
}

std::vector<double> MMMFieldComputer::ComputeBField(
    const std::vector<double>& M,
    const std::vector<double>& obs_points,
    int n_points) const {

    // B = mu_0 * H (in air)
    std::vector<double> H_field = ComputeHField(M, obs_points, n_points);

    for (size_t i = 0; i < H_field.size(); i++) {
        H_field[i] *= constants::MU_0;
    }

    return H_field;
}

//=========================================================================
// MMMNonlinearSolver Implementation
//=========================================================================

MMMNonlinearSolver::MMMNonlinearSolver()
    : total_dof_(0), n_elements_(0),
      tol_(1e-4), max_iter_(100), relax_(1.0), linear_solver_(0),
      has_chi_init_(false) {}

MMMNonlinearSolver::~MMMNonlinearSolver() {}

void MMMNonlinearSolver::SetMatrix(const double* N, const int* dof_offset,
                                    int total_dof, int n_elem) {
    total_dof_ = total_dof;
    n_elements_ = n_elem;

    // Copy interaction matrix (column-major)
    N_.assign(N, N + static_cast<size_t>(total_dof) * total_dof);

    // Copy DOF offsets
    dof_offset_.assign(dof_offset, dof_offset + n_elem + 1);
}

void MMMNonlinearSolver::SetBHCurve(const std::vector<std::array<double, 2>>& bh_curve) {
    bh_curve_ = bh_curve;
}

void MMMNonlinearSolver::SetParameters(double tol, int max_iter,
                                        double relax, int linear_solver) {
    tol_ = tol;
    max_iter_ = max_iter;
    relax_ = relax;
    linear_solver_ = linear_solver;
}

void MMMNonlinearSolver::SetInitialChi(const std::vector<double>& chi_init) {
    chi_init_ = chi_init;
    has_chi_init_ = true;
}

double MMMNonlinearSolver::InterpolateBH(double H_mag) const {
    // Linear interpolation of B-H curve
    if (bh_curve_.empty()) return 0.0;

    if (H_mag <= bh_curve_[0][0]) {
        return bh_curve_[0][1];
    }
    if (H_mag >= bh_curve_.back()[0]) {
        return bh_curve_.back()[1];
    }

    // Binary search for interval
    for (size_t i = 1; i < bh_curve_.size(); i++) {
        if (H_mag <= bh_curve_[i][0]) {
            double H0 = bh_curve_[i-1][0];
            double H1 = bh_curve_[i][0];
            double B0 = bh_curve_[i-1][1];
            double B1 = bh_curve_[i][1];

            double t = (H_mag - H0) / (H1 - H0 + 1e-30);
            return B0 + t * (B1 - B0);
        }
    }

    return bh_curve_.back()[1];
}

double MMMNonlinearSolver::ComputeChiFromH(double H_mag) const {
    // chi = M / H = (B/mu_0 - H) / H = B/(mu_0*H) - 1
    double B = InterpolateBH(H_mag);
    if (H_mag < 1e-30) {
        // At zero field, use initial slope
        if (bh_curve_.size() >= 2 && bh_curve_[0][0] < 1e-30 && bh_curve_[1][0] > 1e-30) {
            double dB = bh_curve_[1][1] - bh_curve_[0][1];
            double dH = bh_curve_[1][0] - bh_curve_[0][0];
            double mu_r = dB / (constants::MU_0 * dH);
            return mu_r - 1.0;
        }
        return 0.0;  // Default: chi=0 (air)
    }

    double mu_r = B / (constants::MU_0 * H_mag);
    double chi = mu_r - 1.0;

    // Clamp to reasonable range
    if (chi < 0.0) chi = 0.0;
    if (chi > 1e6) chi = 1e6;

    return chi;
}

void MMMNonlinearSolver::ComputeLocalH(const std::vector<double>& H_ext,
                                        const std::vector<double>& M) {
    // H_local = H_ext + H_demag = H_ext + N * M
    // But actually for internal H, we need: H_local = H_ext - N * M (demagnetizing field opposes M)
    // Standard formulation: H_internal = H_ext - N_demag * M

    H_local_.resize(n_elements_ * 3);

    for (int e = 0; e < n_elements_; e++) {
        int dof_start = dof_offset_[e];
        int dof_end = dof_offset_[e + 1];
        int n_dof = dof_end - dof_start;

        // Start with external field
        double Hx = 0, Hy = 0, Hz = 0;
        if (n_dof >= 3) {
            // For tetrahedra: H_ext is directly the external field
            // For simplicity, assume H_ext is expanded to total_dof
            Hx = H_ext[dof_start];
            Hy = H_ext[dof_start + 1];
            Hz = H_ext[dof_start + 2];
        }

        // Compute demagnetizing field: H_demag = -N * M
        // The matrix N stores the demagnetization tensor
        // H_internal_i = H_ext_i - sum_j N_ij * M_j

        double Hd_x = 0, Hd_y = 0, Hd_z = 0;
        for (int i = 0; i < total_dof_; i++) {
            // Column-major: N_[i + j * total_dof_]
            Hd_x += N_[dof_start + 0 + i * total_dof_] * M[i];
            Hd_y += N_[dof_start + 1 + i * total_dof_] * M[i];
            Hd_z += N_[dof_start + 2 + i * total_dof_] * M[i];
        }

        // H_internal = H_ext + N * M (N is interaction, includes demagnetization)
        // Note: The sign depends on how N is defined
        H_local_[e * 3 + 0] = Hx + Hd_x;
        H_local_[e * 3 + 1] = Hy + Hd_y;
        H_local_[e * 3 + 2] = Hz + Hd_z;
    }
}

std::vector<double> MMMNonlinearSolver::ExpandH_ext(const std::vector<double>& H_ext) const {
    std::vector<double> H_expanded(total_dof_);

    if (H_ext.size() == 3) {
        // Uniform field: broadcast to all DOF
        for (int e = 0; e < n_elements_; e++) {
            int dof_start = dof_offset_[e];
            int dof_end = dof_offset_[e + 1];
            int n_dof = dof_end - dof_start;

            if (n_dof == 3) {
                H_expanded[dof_start + 0] = H_ext[0];
                H_expanded[dof_start + 1] = H_ext[1];
                H_expanded[dof_start + 2] = H_ext[2];
            } else if (n_dof == 6) {
                // For MSC hexahedra, set H.n for each face
                // Simplified: assume face normals align with axes
                // Face order: -z, +z, -y, +y, -x, +x
                H_expanded[dof_start + 0] = -H_ext[2];  // -z face
                H_expanded[dof_start + 1] = H_ext[2];   // +z face
                H_expanded[dof_start + 2] = -H_ext[1];  // -y face
                H_expanded[dof_start + 3] = H_ext[1];   // +y face
                H_expanded[dof_start + 4] = -H_ext[0];  // -x face
                H_expanded[dof_start + 5] = H_ext[0];   // +x face
            }
        }
    } else if (H_ext.size() == static_cast<size_t>(total_dof_)) {
        H_expanded = H_ext;
    } else {
        throw std::runtime_error("H_ext must have size 3 or total_dof");
    }

    return H_expanded;
}

std::vector<double> MMMNonlinearSolver::GetLocalB() const {
    // B = mu_0 * (H + M)
    // For tetrahedra: B = mu_0 * H + mu_0 * M
    std::vector<double> B_local(n_elements_ * 3);

    for (int e = 0; e < n_elements_; e++) {
        int dof_start = dof_offset_[e];
        int n_dof = dof_offset_[e + 1] - dof_start;

        if (n_dof == 3) {
            // Tetrahedra: M is directly magnetization vector
            B_local[e * 3 + 0] = constants::MU_0 * (H_local_[e * 3 + 0] + M_[dof_start + 0]);
            B_local[e * 3 + 1] = constants::MU_0 * (H_local_[e * 3 + 1] + M_[dof_start + 1]);
            B_local[e * 3 + 2] = constants::MU_0 * (H_local_[e * 3 + 2] + M_[dof_start + 2]);
        } else {
            // MSC: more complex, use H_local directly
            B_local[e * 3 + 0] = constants::MU_0 * H_local_[e * 3 + 0];
            B_local[e * 3 + 1] = constants::MU_0 * H_local_[e * 3 + 1];
            B_local[e * 3 + 2] = constants::MU_0 * H_local_[e * 3 + 2];
        }
    }

    return B_local;
}

MMMNonlinearResult MMMNonlinearSolver::Solve(const std::vector<double>& H_ext) {
    MMMNonlinearResult result;
    result.converged = false;
    result.iterations = 0;
    result.final_error = 1e30;

    if (bh_curve_.empty()) {
        throw std::runtime_error("B-H curve not set");
    }

    // Expand H_ext to full DOF size
    std::vector<double> H_ext_full = ExpandH_ext(H_ext);

    // Initialize chi
    chi_.resize(n_elements_);
    if (has_chi_init_) {
        chi_ = chi_init_;
    } else {
        // Use chi at first B-H point (low-field permeability)
        double chi0 = ComputeChiFromH(1.0);  // Evaluate at H=1 A/m
        for (int e = 0; e < n_elements_; e++) {
            chi_[e] = chi0;
        }
    }

    // Initialize M
    M_.resize(total_dof_, 0.0);

    // Previous B for convergence check
    std::vector<double> B_old(n_elements_ * 3, 0.0);

    // Create linear solver
    MMMSolver linear_solver;
    linear_solver.SetMatrix(N_.data(), dof_offset_.data(), total_dof_, n_elements_);

    // Newton-Raphson iteration
    for (int iter = 0; iter < max_iter_; iter++) {
        result.iterations = iter + 1;

        // Compute 1/chi for solver
        std::vector<double> inv_chi(n_elements_);
        for (int e = 0; e < n_elements_; e++) {
            inv_chi[e] = 1.0 / (chi_[e] + 1e-30);
        }

        // Solve linear system
        std::vector<double> M_new;
        if (linear_solver_ == 0) {
            int info = linear_solver.SolveLU(inv_chi, H_ext_full, M_new, true);
            if (info != 0) {
                throw std::runtime_error("LU solve failed");
            }
        } else {
            int n_iter = linear_solver.SolveBiCGSTAB(inv_chi, H_ext_full, M_new, 1e-10, 1000, true);
            if (n_iter < 0) {
                throw std::runtime_error("BiCGSTAB solve failed");
            }
        }

        // Apply relaxation
        if (relax_ < 1.0 && iter > 0) {
            for (int i = 0; i < total_dof_; i++) {
                M_[i] = (1.0 - relax_) * M_[i] + relax_ * M_new[i];
            }
        } else {
            M_ = M_new;
        }

        // Compute local H field
        ComputeLocalH(H_ext_full, M_);

        // Compute new B field
        std::vector<double> B_new = GetLocalB();

        // Check convergence
        double B_sat = 0.0;
        for (size_t i = 0; i < B_new.size(); i++) {
            B_sat = std::max(B_sat, std::abs(B_new[i]));
        }
        if (B_sat < 1e-10) B_sat = 1.0;  // Avoid division by zero

        result.final_error = CheckConvergence(B_new, B_old, B_sat, n_elements_);

        if (result.final_error < tol_) {
            result.converged = true;
            break;
        }

        // Update chi from local H
        std::vector<double> chi_new(n_elements_);
        for (int e = 0; e < n_elements_; e++) {
            double Hx = H_local_[e * 3 + 0];
            double Hy = H_local_[e * 3 + 1];
            double Hz = H_local_[e * 3 + 2];
            double H_mag = std::sqrt(Hx*Hx + Hy*Hy + Hz*Hz);

            chi_new[e] = ComputeChiFromH(H_mag);
        }

        // Apply relaxation to chi
        if (relax_ < 1.0) {
            for (int e = 0; e < n_elements_; e++) {
                chi_[e] = (1.0 - relax_) * chi_[e] + relax_ * chi_new[e];
            }
        } else {
            chi_ = chi_new;
        }

        B_old = B_new;
    }

    result.M = M_;
    result.chi = chi_;

    return result;
}

}  // namespace mmm
}  // namespace radia
