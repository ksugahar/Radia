// compact_stubs.cpp — Stub implementations for symbols referenced
// by meshclass.cpp/basegeom.cpp but never actually called in our
// compact Netgen workflow (mesh construction + curving + save only).
//
// These stubs resolve linker errors without pulling in the full
// Netgen meshing/optimization/linalg subsystems.

#include <mystdlib.h>
#include "meshing.hpp"

// MarkedTet/MarkedPrism/MarkedIdentification/MarkedTri: complete type
// needed for BisectionInfo destructor (unique_ptr<Array<MarkedTet>>)
namespace netgen {
  class MarkedTet {
  public:
    PointIndex pnums[4];
    int matindex = 0;
    unsigned int marked:2;
    unsigned int flagged:1;
    unsigned int tetedge1:3;
    unsigned int tetedge2:3;
    char faceedges[4];
    bool incorder = false;
    unsigned int order:6;
    int8_t newest_vertex = -1;
    MarkedTet() = default;
  };
  class MarkedPrism {
  public:
    PointIndex pnums[6];
    int matindex = 0;
    unsigned int marked:2;
    MarkedPrism() = default;
  };
  class MarkedIdentification {
  public:
    int np = 0;
    PointIndex pnums[8];
    unsigned int marked:1;
    MarkedIdentification() = default;
  };
  class MarkedTri {
  public:
    PointIndex pnums[3];
    PointGeomInfo IsectP1gi, IsectP2gi;
    unsigned int IsectFlag:1;
    MarkedTri() = default;
  };
}

namespace netgen
{

// DenseMatrix: now using real densemat.cpp (no longer stubbed)
// Mat<3,3> CalcInverse still needed (template specialization, not in densemat.cpp)
void CalcInverse(const Mat<3,3,double>& m, Mat<3,3,double>& inv) {
  double det = m(0,0)*(m(1,1)*m(2,2)-m(1,2)*m(2,1))
              -m(0,1)*(m(1,0)*m(2,2)-m(1,2)*m(2,0))
              +m(0,2)*(m(1,0)*m(2,1)-m(1,1)*m(2,0));
  if (fabs(det) < 1e-40) { inv = 0.0; return; }
  double idet = 1.0/det;
  inv(0,0) =  idet*(m(1,1)*m(2,2)-m(1,2)*m(2,1));
  inv(0,1) = -idet*(m(0,1)*m(2,2)-m(0,2)*m(2,1));
  inv(0,2) =  idet*(m(0,1)*m(1,2)-m(0,2)*m(1,1));
  inv(1,0) = -idet*(m(1,0)*m(2,2)-m(1,2)*m(2,0));
  inv(1,1) =  idet*(m(0,0)*m(2,2)-m(0,2)*m(2,0));
  inv(1,2) = -idet*(m(0,0)*m(1,2)-m(0,2)*m(1,0));
  inv(2,0) =  idet*(m(1,0)*m(2,1)-m(1,1)*m(2,0));
  inv(2,1) = -idet*(m(0,0)*m(2,1)-m(0,1)*m(2,0));
  inv(2,2) =  idet*(m(0,0)*m(1,1)-m(0,1)*m(1,0));
}
double Det(const Mat<2,2,double>& m) {
  return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}

// ============================================================
// Meshing algorithm stubs (from basegeom.cpp virtual methods)
// NetgenGeometry::GenerateMesh/MeshFace/OptimizeSurface are
// virtual methods that we override in CallbackGeometry. The
// base class implementations reference the full mesher.
// ============================================================

// Meshing2: too many virtual methods to stub. Skip entirely.
// Instead, we patch basegeom.cpp to remove mesher dependencies.

// Volume meshing stubs (referenced by basegeom)
MESHING3_RESULT MeshVolume(const MeshingParameters&, Mesh&) { return MESHING3_OK; }
MESHING3_RESULT OptimizeVolume(const MeshingParameters&, Mesh&) { return MESHING3_OK; }
void MeshQuality3d(const Mesh&, NgArray<int,0,int>*) {}

// ============================================================
// Mesh utility stubs (from meshclass.cpp)
// ============================================================
double CalcTetBadness(const Point3d&, const Point3d&, const Point3d&, const Point3d&,
                      double, const MeshingParameters&) { return 0.0; }
int CheckMesh3D(const Mesh&) { return 0; }
int IntersectTriangleTriangle(const Point<3>**, const Point<3>**) { return 0; }
double ComputeCylinderRadius(const Point3d&, const Point3d&, const Point3d&, const Point3d&) { return 0.0; }

// Mesh functions removed from meshclass_patched.cpp (compact build)
void Mesh::BuildElementSearchTree(int) {}
bool Mesh::PointContainedIn3DElementOld(const Point3d&, double* const, ElementIndex, double) const { return false; }
void Mesh::SetLocalH(netgen::Point<3> , netgen::Point<3> , double, int) {}
void Mesh::RestrictLocalH(const Point3d&, double, int) {}
void Mesh::RestrictLocalHLine(const Point3d&, const Point3d&, double, int) {}
void Mesh::LoadLocalMeshSize(const filesystem::path&) {}
void Mesh::SetGlobalH(double h) { hglob = h; }
void Mesh::SetMinimalH(double h) { hmin = h; }
double Mesh::GetH(const Point3d&, int) const { return hglob; }
bool Mesh::PointContainedIn2DElement(const Point3d&, double* const, SurfaceElementIndex, bool) const { return false; }
bool Mesh::PointContainedIn3DElement(const Point3d&, double* const, ElementIndex, double) const { return false; }

// BisectionInfo (mesh refinement)
BisectionInfo::BisectionInfo() {}
BisectionInfo::~BisectionInfo() {}

// AnisotropicClusters (mesh quality) — has reference member
// NOTE: no static Mesh here — static init of Mesh causes issues
AnisotropicClusters::AnisotropicClusters(const Mesh& m) : mesh(m) {}
AnisotropicClusters::~AnisotropicClusters() {}
void AnisotropicClusters::Update() {}

// ============================================================
// Global variables
// ============================================================
DLL_HEADER size_t timestamp = 0;

// ============================================================
// Message handler
// ============================================================
void Ng_PrintDest(const char*) {}

// netrule destructor (needed by Meshing2 destructor via unique_ptr)
netrule::~netrule() {}

// ============================================================
// GeometryRegisterArray (from basegeom.cpp, used by Mesh::Load)
// ============================================================
shared_ptr<NetgenGeometry> GeometryRegisterArray::LoadFromMeshFile(istream&) const { return nullptr; }

// ============================================================
// AdFront stubs (from localh.cpp)
// ============================================================
} // namespace netgen

// Forward-declared in headers but never used in our workflow
namespace netgen {
  bool AdFront2::SameSide(const Point<2>&, const Point<2>&, const FlatArray<int, size_t>*) const { return false; }
  void AdFront3::GetFaceBoundingBox(int, Box3d&) const {}
  int AdFront3::SameSide(const Point<3>&, const Point<3>&, const NgArray<int,0,int>*) const { return 0; }
}
