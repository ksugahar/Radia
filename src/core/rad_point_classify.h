#ifndef RAD_POINT_CLASSIFY_H
#define RAD_POINT_CLASSIFY_H

#include "gmvect.h"
#include <vector>

namespace RadPointClassify {

enum PointClass {
    INSIDE = 0,
    NEAR = 1,
    FAR = 2
};

struct ClassifyResult {
    PointClass classification;
    int nearest_elem_id;
    double distance;
};

struct ElementData {
    TVector3d center;
    double size;
    std::vector<TVector3d> vertices;
    std::vector<std::vector<int>> faces;  // Face topology (vertex indices)
    int num_faces;
};

// Solid Angle Method (Rigorous, ELF-compatible)
double ComputeTriangleSolidAngle(const TVector3d& obs,
                                  const TVector3d& v0,
                                  const TVector3d& v1,
                                  const TVector3d& v2);

bool PointInPolyhedronSolidAngle(const TVector3d& pt,
                                  const std::vector<TVector3d>& vertices,
                                  const std::vector<std::vector<int>>& faces);

bool PointInTetrahedronSolidAngle(const TVector3d& pt,
                                   const TVector3d& v0, const TVector3d& v1,
                                   const TVector3d& v2, const TVector3d& v3);

bool PointInHexahedronSolidAngle(const TVector3d& pt, const TVector3d verts[8]);

// Fast Approximate Methods
bool PointInTetrahedron(const TVector3d& pt,
                        const TVector3d& v0, const TVector3d& v1,
                        const TVector3d& v2, const TVector3d& v3);

bool PointInHexahedron(const TVector3d& pt, const TVector3d verts[8]);

// Utility Functions
double ComputeElementSize(const TVector3d* vertices, int n_verts);

void ComputeElementAABB(const TVector3d* vertices, int n_verts,
                        TVector3d& aabb_min, TVector3d& aabb_max);

void ClassifyPoints(int n_points,
                    const double* points,
                    const std::vector<ElementData>& elements,
                    double near_threshold,
                    std::vector<ClassifyResult>& results);

void ClassifyPointsFromHandle(int n_points,
                              const double* points,
                              int container_handle,
                              double near_threshold,
                              int* classification,
                              int* nearest_elem);

} // namespace RadPointClassify

#endif // RAD_POINT_CLASSIFY_H
