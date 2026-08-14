#include "rad_hdiv_field_evaluator.h"

#include "rad_hdiv_vim.h"

#include <core/taskmanager.hpp>

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace rad_hdiv {
namespace {

struct CompensatedVec3 {
    double sum[3] = {0.0, 0.0, 0.0};
    double correction[3] = {0.0, 0.0, 0.0};

    void Add(int component, double value) {
        const double next = sum[component] + value;
        correction[component] += std::fabs(sum[component]) >= std::fabs(value)
            ? (sum[component] - next) + value
            : (value - next) + sum[component];
        sum[component] = next;
    }

    void Add(const double value[3]) {
        for (int component = 0; component < 3; ++component) Add(component, value[component]);
    }

    void Store(double out[3]) const {
        for (int component = 0; component < 3; ++component)
            out[component] = sum[component] + correction[component];
    }
};

using Vec = std::array<double, 3>;

double Dot(const Vec& a, const Vec& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

Vec Cross(const Vec& a, const Vec& b) {
    return {a[1]*b[2] - a[2]*b[1],
            a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0]};
}

double Norm(const Vec& a) {
    return std::sqrt(Dot(a, a));
}

double Det(const Vec& a, const Vec& b, const Vec& c) {
    return Dot(a, Cross(b, c));
}

struct TetSource {
    double v[4][3]{};
    double coefficient[20]{};
};

struct TriSource {
    double v[3][3]{};
    double sigma0 = 0.0;
    double slope[3]{};
    double hessian[3][3]{};
};

struct CurvedTetSource {
    double nodes[10][3]{};
    double coefficient[4]{}; // 1, xi, eta, zeta
};

struct CurvedTriSource {
    double nodes[6][3]{};
    double coefficient[6]{}; // 1, eta, eta^2, xi, xi*eta, xi^2
};

struct PointSource {
    Vec position{};
    double strength = 0.0;
};

enum class SourceKind : std::uint8_t { Tet, Triangle, CurvedTet, CurvedTriangle, Point };

struct SourceAtom {
    SourceKind kind = SourceKind::Point;
    std::size_t index = 0;
    Vec center{};
    Vec lower{};
    Vec upper{};
    double charge = 0.0;
    Vec raw_first{};
    std::array<double, 9> raw_second{};
};

struct ImageTerm {
    int mask = 0;
    double sign = 1.0;
    double angle = 0.0;                 // rotation about +z; 0 == pure mirror (the historical path)
};

struct TreeNode {
    std::size_t begin = 0;
    std::size_t end = 0;
    int left = -1;
    int right = -1;
    Vec center{};
    Vec lower{};
    Vec upper{};
    double radius = 0.0;
    double charge = 0.0;
    Vec dipole{};
    std::array<double, 9> second{};

    bool IsLeaf() const { return left < 0; }
};

constexpr std::array<double, 4> GL_X = {
    0.06943184420297371, 0.33000947820757187,
    0.6699905217924281, 0.9305681557970262
};
constexpr std::array<double, 4> GL_W = {
    0.1739274225687269, 0.3260725774312731,
    0.3260725774312731, 0.1739274225687269
};

int PolynomialIndex(int ax, int ay, int az) {
    const int degree=ax+ay+az;
    int index=0;
    for(int d=0;d<degree;++d)index+=(d+1)*(d+2)/2;
    for(int x=0;x<ax;++x)index+=degree-x+1;
    return index+ay;
}

double TetDensity(const TetSource& source, const Vec& x) {
    double powers[3][4]={{1.0,x[0],x[0]*x[0],x[0]*x[0]*x[0]},
                         {1.0,x[1],x[1]*x[1],x[1]*x[1]*x[1]},
                         {1.0,x[2],x[2]*x[2],x[2]*x[2]*x[2]}};
    double value=0.0;
    for(int total=0;total<=3;++total)
        for(int ax=0;ax<=total;++ax)
            for(int ay=0;ay<=total-ax;++ay){
                const int az=total-ax-ay;
                value+=source.coefficient[PolynomialIndex(ax,ay,az)]*
                    powers[0][ax]*powers[1][ay]*powers[2][az];
            }
    return value;
}

void AddMoment(SourceAtom& atom, const Vec& x, double dq) {
    atom.charge += dq;
    for (int i = 0; i < 3; ++i) {
        atom.raw_first[i] += dq * x[i];
        for (int j = 0; j < 3; ++j)
            atom.raw_second[3*i+j] += dq * x[i] * x[j];
    }
}

template <std::size_t N>
void SetBounds(SourceAtom& atom, const double (&vertices)[N][3]) {
    for (int k = 0; k < 3; ++k) {
        atom.lower[k] = atom.upper[k] = vertices[0][k];
        for (std::size_t i = 1; i < N; ++i) {
            atom.lower[k] = std::min(atom.lower[k], vertices[i][k]);
            atom.upper[k] = std::max(atom.upper[k], vertices[i][k]);
        }
        atom.center[k] = 0.5 * (atom.lower[k] + atom.upper[k]);
    }
}

SourceAtom MakeTetAtom(const TetSource& source, std::size_t index) {
    SourceAtom atom;
    atom.kind = SourceKind::Tet;
    atom.index = index;
    SetBounds(atom, source.v);
    const Vec e1 = {source.v[1][0]-source.v[0][0], source.v[1][1]-source.v[0][1], source.v[1][2]-source.v[0][2]};
    const Vec e2 = {source.v[2][0]-source.v[0][0], source.v[2][1]-source.v[0][1], source.v[2][2]-source.v[0][2]};
    const Vec e3 = {source.v[3][0]-source.v[0][0], source.v[3][1]-source.v[0][1], source.v[3][2]-source.v[0][2]};
    const double jacobian = std::fabs(Det(e1, e2, e3));
    for (int iu = 0; iu < 4; ++iu) for (int iv = 0; iv < 4; ++iv) for (int iw = 0; iw < 4; ++iw) {
        const double u = GL_X[iu], v = GL_X[iv], w = GL_X[iw];
        const double l0 = (1.0-u)*(1.0-v)*(1.0-w);
        const double l1 = u;
        const double l2 = (1.0-u)*v;
        const double l3 = (1.0-u)*(1.0-v)*w;
        Vec x{};
        for (int k = 0; k < 3; ++k)
            x[k] = l0*source.v[0][k] + l1*source.v[1][k] + l2*source.v[2][k] + l3*source.v[3][k];
        const double rho = TetDensity(source,x);
        const double weight = GL_W[iu]*GL_W[iv]*GL_W[iw]
                            * jacobian*(1.0-u)*(1.0-u)*(1.0-v);
        AddMoment(atom, x, rho*weight);
    }
    return atom;
}

SourceAtom MakeTriAtom(const TriSource& source, std::size_t index) {
    SourceAtom atom;
    atom.kind = SourceKind::Triangle;
    atom.index = index;
    SetBounds(atom, source.v);
    const Vec e1 = {source.v[1][0]-source.v[0][0], source.v[1][1]-source.v[0][1], source.v[1][2]-source.v[0][2]};
    const Vec e2 = {source.v[2][0]-source.v[0][0], source.v[2][1]-source.v[0][1], source.v[2][2]-source.v[0][2]};
    const double jacobian = Norm(Cross(e1, e2));
    for (int iu = 0; iu < 4; ++iu) for (int iv = 0; iv < 4; ++iv) {
        const double u = GL_X[iu], v = GL_X[iv];
        const double l0 = (1.0-u)*(1.0-v);
        const double l1 = u;
        const double l2 = (1.0-u)*v;
        Vec x{};
        for (int k = 0; k < 3; ++k)
            x[k] = l0*source.v[0][k] + l1*source.v[1][k] + l2*source.v[2][k];
        double sigma = source.sigma0 + source.slope[0]*x[0]
                     + source.slope[1]*x[1] + source.slope[2]*x[2];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j)
            sigma += x[i]*source.hessian[i][j]*x[j];
        const double weight = GL_W[iu]*GL_W[iv]*jacobian*(1.0-u);
        AddMoment(atom, x, sigma*weight);
    }
    return atom;
}

SourceAtom MakePointAtom(const PointSource& source, std::size_t index) {
    SourceAtom atom;
    atom.kind = SourceKind::Point;
    atom.index = index;
    atom.center = atom.lower = atom.upper = source.position;
    AddMoment(atom, source.position, source.strength);
    return atom;
}

double CurvedTetDensity(const CurvedTetSource& source, double xi, double eta, double zeta) {
    return source.coefficient[0] + source.coefficient[1]*xi
         + source.coefficient[2]*eta + source.coefficient[3]*zeta;
}

double CurvedTriDensity(const CurvedTriSource& source, double xi, double eta) {
    return source.coefficient[0] + source.coefficient[1]*eta
         + source.coefficient[2]*eta*eta + source.coefficient[3]*xi
         + source.coefficient[4]*xi*eta + source.coefficient[5]*xi*xi;
}

SourceAtom MakeCurvedTetAtom(const CurvedTetSource& source, std::size_t index) {
    SourceAtom atom;
    atom.kind = SourceKind::CurvedTet;
    atom.index = index;
    SetBounds(atom, source.nodes);
    for (int ia = 0; ia < 4; ++ia) for (int ib = 0; ib < 4; ++ib) for (int ic = 0; ic < 4; ++ic) {
        const double a = GL_X[ia], b = GL_X[ib], c = GL_X[ic];
        const double xi = a;
        const double eta = b*(1.0-a);
        const double zeta = c*(1.0-a)*(1.0-b);
        double x[3], jacobian;
        CurvedTetMapMeasure(source.nodes, xi, eta, zeta, x, jacobian);
        const double weight = GL_W[ia]*GL_W[ib]*GL_W[ic]
                            *(1.0-a)*(1.0-a)*(1.0-b)*jacobian;
        AddMoment(atom, {x[0], x[1], x[2]},
                  CurvedTetDensity(source, xi, eta, zeta)*weight);
    }
    return atom;
}

SourceAtom MakeCurvedTriAtom(const CurvedTriSource& source, std::size_t index) {
    SourceAtom atom;
    atom.kind = SourceKind::CurvedTriangle;
    atom.index = index;
    SetBounds(atom, source.nodes);
    for (int iu = 0; iu < 4; ++iu) for (int iv = 0; iv < 4; ++iv) {
        const double xi = GL_X[iu];
        const double eta = GL_X[iv]*(1.0-xi);
        double x[3], jacobian;
        CurvedTriMapMeasure(source.nodes, xi, eta, x, jacobian);
        const double weight = GL_W[iu]*GL_W[iv]*(1.0-xi)*jacobian;
        AddMoment(atom, {x[0], x[1], x[2]}, CurvedTriDensity(source, xi, eta)*weight);
    }
    return atom;
}

} // namespace

struct HDivFieldEvaluator::Impl {
    FieldEvaluatorOptions options;
    std::vector<TetSource> tets;
    std::vector<TriSource> triangles;
    std::vector<CurvedTetSource> curved_tets;
    std::vector<CurvedTriSource> curved_triangles;
    std::vector<PointSource> points;
    std::vector<double> gauss_points;
    std::vector<double> gauss_weights;
    std::vector<SourceAtom> atoms;
    std::vector<ImageTerm> images;
    std::vector<std::size_t> order;
    std::vector<TreeNode> nodes;
    Vec lower{};
    Vec upper{};
    mutable std::atomic<int> last_algorithm{0};

    void ValidateOptions() {
        if (options.leaf_size < 1) throw std::invalid_argument("HDivFieldEvaluator: leaf_size must be >= 1");
        if (!(options.theta > 0.0 && options.theta < 1.0))
            throw std::invalid_argument("HDivFieldEvaluator: theta must be in (0,1)");
        if (options.tree_min_sources < 2)
            throw std::invalid_argument("HDivFieldEvaluator: tree_min_sources must be >= 2");
        if (options.auto_min_work < 1)
            throw std::invalid_argument("HDivFieldEvaluator: auto_min_work must be >= 1");
        if (!(options.tree_relative_tolerance > 0.0 && options.tree_relative_tolerance < 1.0))
            throw std::invalid_argument("HDivFieldEvaluator: tree_relative_tolerance must be in (0,1)");
        if (options.probe_count < 1)
            throw std::invalid_argument("HDivFieldEvaluator: probe_count must be >= 1");
    }

    void SetImages(std::vector<int> masks, std::vector<double> signs) {
        if (masks.size() != signs.size())
            throw std::invalid_argument("HDivFieldEvaluator: image_masks and image_signs size mismatch");
        images.reserve(masks.size());
        for (std::size_t i = 0; i < masks.size(); ++i) {
            if (masks[i] < 0 || masks[i] > 7)
                throw std::invalid_argument("HDivFieldEvaluator: each image mask must be in [0,7] "
                                            "(0 = pure rotation, requires a non-zero rotation angle)");
            if (!std::isfinite(signs[i]))
                throw std::invalid_argument("HDivFieldEvaluator: image signs must be finite");
            images.push_back({masks[i], signs[i], 0.0});
        }
    }

    void SetImageRotations(std::vector<double> angles) {
        if (angles.empty()) {
            for (ImageTerm& image : images) image.angle = 0.0;
        } else {
            if (angles.size() != images.size())
                throw std::invalid_argument(
                    "HDivFieldEvaluator: image rotation angles must match the image count");
            for (std::size_t i = 0; i < angles.size(); ++i) {
                if (!std::isfinite(angles[i]))
                    throw std::invalid_argument("HDivFieldEvaluator: image rotation angles must be finite");
                images[i].angle = angles[i];
            }
        }
        for (const ImageTerm& image : images)
            if (image.mask == 0 && image.angle == 0.0)
                throw std::invalid_argument(
                    "HDivFieldEvaluator: an image with mask 0 and rotation angle 0 is the IDENTITY -- it "
                    "would double the direct term; give it a mirror mask or a non-zero rotation angle");
    }

    // T^-1 on an eval point: mirror on the mask axes, then rotate by -angle (see the charge Gram's
    // ImageEvalPoint -- a mirror is an involution so the historical code could not distinguish the two;
    // a rotation is not, and using +angle here silently reads the wrong image).
    static void ImageInversePoint(const ImageTerm& image, const double v[3], double o[3]) {
        double t[3] = {v[0], v[1], v[2]};
        for (int axis = 0; axis < 3; ++axis) if (image.mask & (1 << axis)) t[axis] = -t[axis];
        if (image.angle != 0.0) {
            const double a = -image.angle, c = std::cos(a), s = std::sin(a);
            const double x = t[0], y = t[1];
            t[0] = c*x - s*y;  t[1] = s*x + c*y;
        }
        o[0] = t[0]; o[1] = t[1]; o[2] = t[2];
    }

    // T forward on a field vector: rotate by +angle, then mirror.
    static void ImageForwardVector(const ImageTerm& image, const double v[3], double o[3]) {
        double t[3] = {v[0], v[1], v[2]};
        if (image.angle != 0.0) {
            const double c = std::cos(image.angle), s = std::sin(image.angle);
            const double x = t[0], y = t[1];
            t[0] = c*x - s*y;  t[1] = s*x + c*y;
        }
        for (int axis = 0; axis < 3; ++axis) if (image.mask & (1 << axis)) t[axis] = -t[axis];
        o[0] = t[0]; o[1] = t[1]; o[2] = t[2];
    }

    int BuildNode(std::size_t begin, std::size_t end, int depth) {
        TreeNode node;
        node.begin = begin;
        node.end = end;
        const double inf = std::numeric_limits<double>::infinity();
        node.lower = {inf, inf, inf};
        node.upper = {-inf, -inf, -inf};
        Vec raw_first{};
        std::array<double, 9> raw_second{};
        for (std::size_t pos = begin; pos < end; ++pos) {
            const SourceAtom& atom = atoms[order[pos]];
            node.charge += atom.charge;
            for (int k = 0; k < 3; ++k) {
                node.lower[k] = std::min(node.lower[k], atom.lower[k]);
                node.upper[k] = std::max(node.upper[k], atom.upper[k]);
                raw_first[k] += atom.raw_first[k];
            }
            for (int k = 0; k < 9; ++k) raw_second[k] += atom.raw_second[k];
        }
        for (int k = 0; k < 3; ++k) node.center[k] = 0.5*(node.lower[k] + node.upper[k]);
        Vec half = {0.5*(node.upper[0]-node.lower[0]),
                    0.5*(node.upper[1]-node.lower[1]),
                    0.5*(node.upper[2]-node.lower[2])};
        node.radius = Norm(half);
        for (int i = 0; i < 3; ++i) {
            node.dipole[i] = raw_first[i] - node.charge*node.center[i];
            for (int j = 0; j < 3; ++j) {
                node.second[3*i+j] = raw_second[3*i+j]
                    - node.center[i]*raw_first[j] - raw_first[i]*node.center[j]
                    + node.charge*node.center[i]*node.center[j];
            }
        }
        const int index = static_cast<int>(nodes.size());
        nodes.push_back(node);
        const std::size_t count = end - begin;
        if (count <= static_cast<std::size_t>(options.leaf_size) || depth >= 64) return index;

        Vec center_lower = {inf, inf, inf};
        Vec center_upper = {-inf, -inf, -inf};
        for (std::size_t pos = begin; pos < end; ++pos) {
            const Vec& c = atoms[order[pos]].center;
            for (int k = 0; k < 3; ++k) {
                center_lower[k] = std::min(center_lower[k], c[k]);
                center_upper[k] = std::max(center_upper[k], c[k]);
            }
        }
        int axis = 0;
        if (center_upper[1]-center_lower[1] > center_upper[axis]-center_lower[axis]) axis = 1;
        if (center_upper[2]-center_lower[2] > center_upper[axis]-center_lower[axis]) axis = 2;
        if (!(center_upper[axis] > center_lower[axis])) return index;
        const std::size_t mid = begin + count/2;
        std::nth_element(order.begin()+begin, order.begin()+mid, order.begin()+end,
            [&](std::size_t a, std::size_t b) { return atoms[a].center[axis] < atoms[b].center[axis]; });
        const int left = BuildNode(begin, mid, depth+1);
        const int right = BuildNode(mid, end, depth+1);
        nodes[index].left = left;
        nodes[index].right = right;
        return index;
    }

    void BuildTree() {
        if (atoms.empty()) throw std::invalid_argument("HDivFieldEvaluator: at least one source is required");
        order.resize(atoms.size());
        std::iota(order.begin(), order.end(), std::size_t(0));
        nodes.reserve(2*atoms.size());
        BuildNode(0, atoms.size(), 0);
        lower = nodes[0].lower;
        upper = nodes[0].upper;
    }

    void AddCurvedTet(const CurvedTetSource& source, const double r[3], double out[3]) const {
        double xi0[3];
        ClosestRefTet(source.nodes, r, xi0);
        static const double corners[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        static const int faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
        const int nq = static_cast<int>(gauss_points.size());
        for (int face = 0; face < 4; ++face) {
            for (int lead = 0; lead < 3; ++lead) {
                const double* b1 = corners[faces[face][lead]];
                const double* b2 = corners[faces[face][(lead+1)%3]];
                const double* b3 = corners[faces[face][(lead+2)%3]];
                double d1[3], d2[3], d3[3], e21[3], e32[3];
                for (int k = 0; k < 3; ++k) {
                    d1[k] = b1[k]-xi0[k]; d2[k] = b2[k]-xi0[k]; d3[k] = b3[k]-xi0[k];
                    e21[k] = b2[k]-b1[k]; e32[k] = b3[k]-b2[k];
                }
                const double determinant = d1[0]*(d2[1]*d3[2]-d2[2]*d3[1])
                                         + d1[1]*(d2[2]*d3[0]-d2[0]*d3[2])
                                         + d1[2]*(d2[0]*d3[1]-d2[1]*d3[0]);
                if (std::fabs(determinant) < 1e-300) continue;
                for (int ia = 0; ia < nq; ++ia) {
                    const double u = gauss_points[static_cast<std::size_t>(ia)];
                    for (int ib = 0; ib < nq; ++ib) {
                        const double v = gauss_points[static_cast<std::size_t>(ib)];
                        for (int ic = 0; ic < nq; ++ic) {
                            const double w = gauss_points[static_cast<std::size_t>(ic)];
                            double xi[3];
                            for (int k = 0; k < 3; ++k)
                                xi[k] = xi0[k] + u*(d1[k] + v*(e21[k] + w*e32[k]));
                            double x[3], jacobian;
                            CurvedTetMapMeasure(source.nodes, xi[0], xi[1], xi[2], x, jacobian);
                            const double dx = r[0]-x[0], dy = r[1]-x[1], dz = r[2]-x[2];
                            const double distance2 = dx*dx + dy*dy + dz*dz;
                            if (distance2 <= 1e-300) continue;
                            const double weight = gauss_weights[static_cast<std::size_t>(ia)]
                                                * gauss_weights[static_cast<std::size_t>(ib)]
                                                * gauss_weights[static_cast<std::size_t>(ic)] / 3.0
                                                * u*u*v*determinant*jacobian;
                            const double scale = CurvedTetDensity(source, xi[0], xi[1], xi[2])
                                               * weight/(distance2*std::sqrt(distance2));
                            out[0] += scale*dx; out[1] += scale*dy; out[2] += scale*dz;
                        }
                    }
                }
            }
        }
    }

    void AddCurvedTriangle(const CurvedTriSource& source, const double r[3], double out[3]) const {
        double xi0[2];
        ClosestRefTri(source.nodes, r, xi0);
        static const double corners[3][2] = {{0,0},{1,0},{0,1}};
        const int nq = static_cast<int>(gauss_points.size());
        for (int edge = 0; edge < 3; ++edge) {
            const double* a = corners[edge];
            const double* b = corners[(edge+1)%3];
            const double e1x = a[0]-xi0[0], e1y = a[1]-xi0[1];
            const double e2x = b[0]-xi0[0], e2y = b[1]-xi0[1];
            const double determinant = e1x*e2y-e1y*e2x;
            for (int iu = 0; iu < nq; ++iu) {
                const double u = gauss_points[static_cast<std::size_t>(iu)];
                for (int iv = 0; iv < nq; ++iv) {
                    const double v = gauss_points[static_cast<std::size_t>(iv)];
                    const double xi = xi0[0] + u*e1x + u*v*(e2x-e1x);
                    const double eta = xi0[1] + u*e1y + u*v*(e2y-e1y);
                    double x[3], jacobian;
                    CurvedTriMapMeasure(source.nodes, xi, eta, x, jacobian);
                    const double dx = r[0]-x[0], dy = r[1]-x[1], dz = r[2]-x[2];
                    const double distance2 = dx*dx + dy*dy + dz*dz;
                    if (distance2 <= 1e-300) continue;
                    const double weight = gauss_weights[static_cast<std::size_t>(iu)]
                                        * gauss_weights[static_cast<std::size_t>(iv)]
                                        * u*determinant*jacobian;
                    const double scale = CurvedTriDensity(source, xi, eta)
                                       * weight/(distance2*std::sqrt(distance2));
                    out[0] += scale*dx; out[1] += scale*dy; out[2] += scale*dz;
                }
            }
        }
    }

    void AddExact(const SourceAtom& atom, const double r[3], double out[3]) const {
        if (atom.kind == SourceKind::Tet) {
            const TetSource& source = tets[atom.index];
            double value[3];
            TetVolFieldCubic(source.v, r, source.coefficient, value);
            for (int k = 0; k < 3; ++k) out[k] += value[k];
        } else if (atom.kind == SourceKind::Triangle) {
            const TriSource& source = triangles[atom.index];
            double value[3];
            QuadTriField(source.v, r, source.sigma0, source.slope, source.hessian, value);
            for (int k = 0; k < 3; ++k) out[k] += value[k];
        } else if (atom.kind == SourceKind::CurvedTet) {
            AddCurvedTet(curved_tets[atom.index], r, out);
        } else if (atom.kind == SourceKind::CurvedTriangle) {
            AddCurvedTriangle(curved_triangles[atom.index], r, out);
        } else {
            const PointSource& source = points[atom.index];
            const double dx = r[0]-source.position[0];
            const double dy = r[1]-source.position[1];
            const double dz = r[2]-source.position[2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            if (r2 <= 1e-300) return;
            const double scale = source.strength/(r2*std::sqrt(r2));
            out[0] += scale*dx; out[1] += scale*dy; out[2] += scale*dz;
        }
    }

    void AddMultipole(const TreeNode& node, const double r[3], double out[3]) const {
        const Vec R = {r[0]-node.center[0], r[1]-node.center[1], r[2]-node.center[2]};
        const double r2 = Dot(R, R);
        if (r2 <= 1e-300) return;
        const double inv_r = 1.0/std::sqrt(r2);
        const double inv_r3 = inv_r/r2;
        const double inv_r5 = inv_r3/r2;
        const double inv_r7 = inv_r5/r2;
        const double pR = Dot(node.dipole, R);
        Vec second_R{};
        double trace = 0.0;
        for (int i = 0; i < 3; ++i) {
            trace += node.second[3*i+i];
            for (int j = 0; j < 3; ++j) second_R[i] += node.second[3*i+j]*R[j];
        }
        const double R_second_R = Dot(R, second_R);
        for (int k = 0; k < 3; ++k) {
            out[k] += node.charge*R[k]*inv_r3;
            out[k] += 3.0*R[k]*pR*inv_r5 - node.dipole[k]*inv_r3;
            out[k] += 7.5*R[k]*R_second_R*inv_r7
                    - 1.5*(trace*R[k] + 2.0*second_R[k])*inv_r5;
        }
    }

    void AddTree(int node_index, const double r[3], double out[3]) const {
        const TreeNode& node = nodes[node_index];
        const Vec R = {r[0]-node.center[0], r[1]-node.center[1], r[2]-node.center[2]};
        const double distance = Norm(R);
        if (!node.IsLeaf() && distance > 0.0 && node.radius/distance <= options.theta) {
            AddMultipole(node, r, out);
            return;
        }
        if (node.IsLeaf()) {
            for (std::size_t pos = node.begin; pos < node.end; ++pos) AddExact(atoms[order[pos]], r, out);
            return;
        }
        AddTree(node.left, r, out);
        AddTree(node.right, r, out);
    }

    void EvaluateBase(const double r[3], double out[3], HDivFieldEvaluator::Algorithm algorithm) const {
        out[0] = out[1] = out[2] = 0.0;
        if (algorithm == HDivFieldEvaluator::Algorithm::Tree) {
            AddTree(0, r, out);
        } else {
            CompensatedVec3 accumulated;
            for (const TetSource& source : tets) {
                double value[3]; TetVolFieldCubic(source.v, r, source.coefficient, value);
                accumulated.Add(value);
            }
            for (const TriSource& source : triangles) {
                double value[3]; QuadTriField(source.v, r, source.sigma0, source.slope, source.hessian, value);
                accumulated.Add(value);
            }
            for (const CurvedTetSource& source : curved_tets) {
                double value[3] = {0.0, 0.0, 0.0};
                AddCurvedTet(source, r, value);
                accumulated.Add(value);
            }
            for (const CurvedTriSource& source : curved_triangles) {
                double value[3] = {0.0, 0.0, 0.0};
                AddCurvedTriangle(source, r, value);
                accumulated.Add(value);
            }
            for (const PointSource& source : points) {
                const double dx = r[0]-source.position[0];
                const double dy = r[1]-source.position[1];
                const double dz = r[2]-source.position[2];
                const double r2 = dx*dx + dy*dy + dz*dz;
                if (r2 <= 1e-300) continue;
                const double scale = source.strength/(r2*std::sqrt(r2));
                accumulated.Add(0, scale*dx);
                accumulated.Add(1, scale*dy);
                accumulated.Add(2, scale*dz);
            }
            accumulated.Store(out);
        }
    }

    void EvaluatePhysical(const double r[3], double out[3], HDivFieldEvaluator::Algorithm algorithm) const {
        double base[3];
        EvaluateBase(r, base, algorithm);
        CompensatedVec3 accumulated;
        accumulated.Add(base);
        for (const ImageTerm& image : images) {
            double reflected[3];
            ImageInversePoint(image, r, reflected);
            double value[3], mapped[3];
            EvaluateBase(reflected, value, algorithm);
            ImageForwardVector(image, value, mapped);
            for (int axis = 0; axis < 3; ++axis)
                accumulated.Add(axis, image.sign*mapped[axis]);
        }
        accumulated.Store(out);
    }

    bool TreePassesProbe(const double* observations, std::size_t n_observations) const {
        const std::size_t count = std::min<std::size_t>(
            n_observations, static_cast<std::size_t>(options.probe_count));
        if (count == 0) return false;
        std::vector<std::size_t> indices;
        indices.reserve(count);
        auto add_index = [&](std::size_t index) {
            if (indices.size() < count
                    && std::find(indices.begin(), indices.end(), index) == indices.end())
                indices.push_back(index);
        };
        add_index(0);
        add_index(n_observations-1);
        std::size_t min_axis[3] = {0, 0, 0}, max_axis[3] = {0, 0, 0};
        std::size_t nearest = 0, farthest = 0;
        double nearest_d2 = std::numeric_limits<double>::infinity(), farthest_d2 = -1.0;
        for (std::size_t i = 0; i < n_observations; ++i) {
            const double* r = observations+3*i;
            for (int axis = 0; axis < 3; ++axis) {
                if (r[axis] < observations[3*min_axis[axis]+axis]) min_axis[axis] = i;
                if (r[axis] > observations[3*max_axis[axis]+axis]) max_axis[axis] = i;
            }
            double distance2 = 0.0;
            for (int axis = 0; axis < 3; ++axis) {
                const double delta = r[axis] < lower[axis] ? lower[axis]-r[axis]
                                   : r[axis] > upper[axis] ? r[axis]-upper[axis] : 0.0;
                distance2 += delta*delta;
            }
            if (distance2 < nearest_d2) { nearest_d2 = distance2; nearest = i; }
            if (distance2 > farthest_d2) { farthest_d2 = distance2; farthest = i; }
        }
        add_index(nearest);
        add_index(farthest);
        for (int axis = 0; axis < 3; ++axis) { add_index(min_axis[axis]); add_index(max_axis[axis]); }
        const std::size_t grid = std::max<std::size_t>(count*2, 2);
        for (std::size_t i = 0; i < grid && indices.size() < count; ++i)
            add_index(i*(n_observations-1)/(grid-1));
        for (std::size_t i = 0; i < n_observations && indices.size() < count; ++i) add_index(i);
        const std::size_t probe_count = indices.size();
        std::vector<double> direct(3*probe_count), tree(3*probe_count);
        const auto direct_start = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < probe_count; ++i)
            EvaluatePhysical(observations+3*indices[i], direct.data()+3*i,
                             HDivFieldEvaluator::Algorithm::Direct);
        const auto direct_end = std::chrono::steady_clock::now();
        for (std::size_t i = 0; i < probe_count; ++i)
            EvaluatePhysical(observations+3*indices[i], tree.data()+3*i,
                             HDivFieldEvaluator::Algorithm::Tree);
        const auto tree_end = std::chrono::steady_clock::now();
        double direct_scale2 = 0.0, error2 = 0.0;
        for (std::size_t i = 0; i < probe_count; ++i) {
            double value2 = 0.0, delta2 = 0.0;
            for (int k = 0; k < 3; ++k) {
                value2 += direct[3*i+k]*direct[3*i+k];
                const double delta = tree[3*i+k]-direct[3*i+k];
                delta2 += delta*delta;
            }
            direct_scale2 = std::max(direct_scale2, value2);
            error2 = std::max(error2, delta2);
        }
        const double relative_error = std::sqrt(error2)/std::max(std::sqrt(direct_scale2), 1e-300);
        const auto direct_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(direct_end-direct_start).count();
        const auto tree_ns = std::chrono::duration_cast<std::chrono::nanoseconds>(tree_end-direct_end).count();
        return relative_error <= options.tree_relative_tolerance && tree_ns*20 < direct_ns*19;
    }
};

HDivFieldEvaluator::HDivFieldEvaluator(std::unique_ptr<Impl> impl) : m_impl(std::move(impl)) {}
void HDivFieldEvaluator::SetImageRotations(std::vector<double> angles)
{
    m_impl->SetImageRotations(std::move(angles));
}

HDivFieldEvaluator::~HDivFieldEvaluator() = default;

std::shared_ptr<HDivFieldEvaluator> HDivFieldEvaluator::FromTet(
    std::vector<double> volume, std::vector<double> surface,
    std::vector<int> image_masks, std::vector<double> image_signs,
    const FieldEvaluatorOptions& options) {
    if (volume.size()%16 != 0 || surface.size()%22 != 0)
        throw std::invalid_argument("HDivFieldEvaluator.from_tet: volume/surface shape mismatch");
    auto impl = std::make_unique<Impl>();
    impl->options = options;
    impl->ValidateOptions();
    impl->SetImages(std::move(image_masks), std::move(image_signs));
    impl->tets.resize(volume.size()/16);
    for (std::size_t e = 0; e < impl->tets.size(); ++e) {
        TetSource& source = impl->tets[e];
        const double* block = volume.data()+16*e;
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k) source.v[i][k] = block[3*i+k];
        source.coefficient[0] = block[12];
        source.coefficient[PolynomialIndex(1,0,0)] = block[13];
        source.coefficient[PolynomialIndex(0,1,0)] = block[14];
        source.coefficient[PolynomialIndex(0,0,1)] = block[15];
        impl->atoms.push_back(MakeTetAtom(source, e));
    }
    impl->triangles.resize(surface.size()/22);
    for (std::size_t e = 0; e < impl->triangles.size(); ++e) {
        TriSource& source = impl->triangles[e];
        const double* block = surface.data()+22*e;
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) source.v[i][k] = block[3*i+k];
        source.sigma0 = block[9];
        for (int k = 0; k < 3; ++k) source.slope[k] = block[10+k];
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) source.hessian[i][k] = block[13+3*i+k];
        impl->atoms.push_back(MakeTriAtom(source, e));
    }
    impl->BuildTree();
    return std::shared_ptr<HDivFieldEvaluator>(new HDivFieldEvaluator(std::move(impl)));
}

std::shared_ptr<HDivFieldEvaluator> HDivFieldEvaluator::FromPolynomialTet(
    std::vector<double> volume, std::vector<double> surface,
    std::vector<int> image_masks, std::vector<double> image_signs,
    const FieldEvaluatorOptions& options) {
    if (volume.size()%32 != 0 || surface.size()%22 != 0)
        throw std::invalid_argument(
            "HDivFieldEvaluator.from_polynomial_tet: volume/surface shape mismatch");
    auto impl = std::make_unique<Impl>();
    impl->options = options;
    impl->ValidateOptions();
    impl->SetImages(std::move(image_masks), std::move(image_signs));
    impl->tets.resize(volume.size()/32);
    for (std::size_t e = 0; e < impl->tets.size(); ++e) {
        TetSource& source = impl->tets[e];
        const double* block = volume.data()+32*e;
        for (int i = 0; i < 4; ++i) for (int k = 0; k < 3; ++k)
            source.v[i][k] = block[3*i+k];
        std::copy_n(block+12,20,source.coefficient);
        impl->atoms.push_back(MakeTetAtom(source,e));
    }
    impl->triangles.resize(surface.size()/22);
    for (std::size_t e = 0; e < impl->triangles.size(); ++e) {
        TriSource& source = impl->triangles[e];
        const double* block = surface.data()+22*e;
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k)
            source.v[i][k] = block[3*i+k];
        source.sigma0 = block[9];
        for (int k = 0; k < 3; ++k) source.slope[k] = block[10+k];
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k)
            source.hessian[i][k] = block[13+3*i+k];
        impl->atoms.push_back(MakeTriAtom(source,e));
    }
    impl->BuildTree();
    return std::shared_ptr<HDivFieldEvaluator>(new HDivFieldEvaluator(std::move(impl)));
}

std::shared_ptr<HDivFieldEvaluator> HDivFieldEvaluator::FromCloud(
    std::vector<double> xyz, std::vector<double> strength,
    std::vector<int> image_masks, std::vector<double> image_signs,
    const FieldEvaluatorOptions& options) {
    if (xyz.size() != 3*strength.size())
        throw std::invalid_argument("HDivFieldEvaluator.from_cloud: xyz/strength shape mismatch");
    auto impl = std::make_unique<Impl>();
    impl->options = options;
    impl->ValidateOptions();
    impl->SetImages(std::move(image_masks), std::move(image_signs));
    impl->points.resize(strength.size());
    for (std::size_t i = 0; i < strength.size(); ++i) {
        impl->points[i].position = {xyz[3*i], xyz[3*i+1], xyz[3*i+2]};
        impl->points[i].strength = strength[i];
        impl->atoms.push_back(MakePointAtom(impl->points[i], i));
    }
    impl->BuildTree();
    return std::shared_ptr<HDivFieldEvaluator>(new HDivFieldEvaluator(std::move(impl)));
}

std::shared_ptr<HDivFieldEvaluator> HDivFieldEvaluator::FromCurvedTet(
    std::vector<double> volume, std::vector<double> surface,
    std::vector<double> gauss_points, std::vector<double> gauss_weights,
    std::vector<int> image_masks, std::vector<double> image_signs,
    const FieldEvaluatorOptions& options) {
    if (volume.size()%34 != 0 || surface.size()%24 != 0)
        throw std::invalid_argument("HDivFieldEvaluator.from_curved_tet: volume/surface shape mismatch");
    if (gauss_points.empty() || gauss_points.size() != gauss_weights.size())
        throw std::invalid_argument("HDivFieldEvaluator.from_curved_tet: invalid Gauss rule");
    auto impl = std::make_unique<Impl>();
    impl->options = options;
    impl->ValidateOptions();
    impl->SetImages(std::move(image_masks), std::move(image_signs));
    impl->gauss_points = std::move(gauss_points);
    impl->gauss_weights = std::move(gauss_weights);
    impl->curved_tets.resize(volume.size()/34);
    for (std::size_t e = 0; e < impl->curved_tets.size(); ++e) {
        CurvedTetSource& source = impl->curved_tets[e];
        const double* block = volume.data()+34*e;
        for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k)
            source.nodes[i][k] = block[3*i+k];
        for (int i = 0; i < 4; ++i) source.coefficient[i] = block[30+i];
        impl->atoms.push_back(MakeCurvedTetAtom(source, e));
    }
    impl->curved_triangles.resize(surface.size()/24);
    for (std::size_t e = 0; e < impl->curved_triangles.size(); ++e) {
        CurvedTriSource& source = impl->curved_triangles[e];
        const double* block = surface.data()+24*e;
        for (int i = 0; i < 6; ++i) for (int k = 0; k < 3; ++k)
            source.nodes[i][k] = block[3*i+k];
        for (int i = 0; i < 6; ++i) source.coefficient[i] = block[18+i];
        impl->atoms.push_back(MakeCurvedTriAtom(source, e));
    }
    impl->BuildTree();
    return std::shared_ptr<HDivFieldEvaluator>(new HDivFieldEvaluator(std::move(impl)));
}

void HDivFieldEvaluator::Evaluate(const double* observations, std::size_t n_observations,
                                  double* output, Algorithm algorithm) const {
    if (!observations || !output) {
        if (n_observations == 0) return;
        throw std::invalid_argument("HDivFieldEvaluator.field: null observation/output buffer");
    }
    const bool automatic = algorithm == Algorithm::Auto;
    if (automatic) algorithm = AlgorithmFor(n_observations);
    if (automatic && algorithm == Algorithm::Tree && !m_impl->TreePassesProbe(observations, n_observations))
        algorithm = Algorithm::Direct;
    m_impl->last_algorithm.store(algorithm == Algorithm::Tree ? 1 : 0, std::memory_order_relaxed);

    // HEX/WEDGE are retained as a point cloud.  For a small target batch, the
    // point-source loop is cheaper than opening a fresh TaskManager region,
    // especially on the first field call after a large solve.  Keep large field
    // maps parallel; express the threshold in physical+IMA interactions rather
    // than target count so reduced-domain solves follow the same cost model.
    constexpr std::size_t kSerialCloudInteractions = 4'000'000;
    const long double cloud_interactions = static_cast<long double>(m_impl->points.size())
        * static_cast<long double>(n_observations)
        * static_cast<long double>(1+m_impl->images.size());
    if (!m_impl->points.empty()
            && cloud_interactions <= static_cast<long double>(kSerialCloudInteractions)) {
        for (std::size_t index = 0; index < n_observations; ++index) {
            const double* r = observations+3*index;
            double total[3];
            m_impl->EvaluatePhysical(r, total, algorithm);
            output[3*index] = total[0];
            output[3*index+1] = total[1];
            output[3*index+2] = total[2];
        }
        return;
    }
    ngcore::RegionTaskManager task_manager;
    ngcore::ParallelFor(ngcore::IntRange(n_observations), [&](std::size_t index) {
        const double* r = observations+3*index;
        double total[3];
        m_impl->EvaluatePhysical(r, total, algorithm);
        output[3*index] = total[0];
        output[3*index+1] = total[1];
        output[3*index+2] = total[2];
    });
}

void HDivFieldEvaluator::EvaluateSerial(const double* observations, std::size_t n_observations,
                                        double* output, Algorithm algorithm) const {
    if (!observations || !output) {
        if (n_observations == 0) return;
        throw std::invalid_argument("HDivFieldEvaluator.field: null observation/output buffer");
    }
    const bool automatic = algorithm == Algorithm::Auto;
    if (automatic) algorithm = AlgorithmFor(n_observations);
    if (automatic && algorithm == Algorithm::Tree && !m_impl->TreePassesProbe(observations, n_observations))
        algorithm = Algorithm::Direct;
    m_impl->last_algorithm.store(algorithm == Algorithm::Tree ? 1 : 0, std::memory_order_relaxed);
    for (std::size_t index = 0; index < n_observations; ++index) {
        const double* r = observations+3*index;
        double total[3];
        m_impl->EvaluatePhysical(r, total, algorithm);
        output[3*index] = total[0];
        output[3*index+1] = total[1];
        output[3*index+2] = total[2];
    }
}

HDivFieldEvaluator::Algorithm HDivFieldEvaluator::AlgorithmFor(std::size_t n_observations) const {
    // IMA is a roundoff-level reduced/full contract.  Different full/reduced
    // source trees need not share the same low-rank truncation, so automatic
    // field evaluation stays on the exact analytic/cloud sum when images exist.
    if (!m_impl->images.empty()) return Algorithm::Direct;
    const std::size_t source_count = SourceCount();
    if (source_count < m_impl->options.tree_min_sources) return Algorithm::Direct;
    const long double work = static_cast<long double>(source_count)
                           * static_cast<long double>(n_observations)
                           * static_cast<long double>(1+m_impl->images.size());
    return work >= static_cast<long double>(m_impl->options.auto_min_work)
        ? Algorithm::Tree : Algorithm::Direct;
}

HDivFieldEvaluator::Algorithm HDivFieldEvaluator::ParseAlgorithm(const std::string& name) {
    if (name == "auto") return Algorithm::Auto;
    if (name == "direct") return Algorithm::Direct;
    if (name == "tree") return Algorithm::Tree;
    throw std::invalid_argument("HDivFieldEvaluator.field: algorithm must be 'auto', 'direct', or 'tree'");
}

const char* HDivFieldEvaluator::AlgorithmName(Algorithm algorithm) {
    if (algorithm == Algorithm::Auto) return "auto";
    if (algorithm == Algorithm::Direct) return "direct";
    return "tree";
}

std::size_t HDivFieldEvaluator::SourceCount() const { return m_impl->atoms.size(); }
std::size_t HDivFieldEvaluator::ImageCount() const { return m_impl->images.size(); }
std::size_t HDivFieldEvaluator::TreeNodeCount() const { return m_impl->nodes.size(); }
int HDivFieldEvaluator::LeafSize() const { return m_impl->options.leaf_size; }
double HDivFieldEvaluator::Theta() const { return m_impl->options.theta; }
std::size_t HDivFieldEvaluator::TreeMinSources() const { return m_impl->options.tree_min_sources; }
std::size_t HDivFieldEvaluator::AutoMinWork() const { return m_impl->options.auto_min_work; }
double HDivFieldEvaluator::TreeRelativeTolerance() const { return m_impl->options.tree_relative_tolerance; }
int HDivFieldEvaluator::ProbeCount() const { return m_impl->options.probe_count; }
const char* HDivFieldEvaluator::SourceRepresentation() const {
    if (!m_impl->curved_tets.empty() || !m_impl->curved_triangles.empty()) return "curved-element-exact";
    if (!m_impl->tets.empty() || !m_impl->triangles.empty()) return "analytic-tet";
    return "element-cloud";
}
HDivFieldEvaluator::Algorithm HDivFieldEvaluator::LastAlgorithm() const {
    return m_impl->last_algorithm.load(std::memory_order_relaxed) == 1 ? Algorithm::Tree : Algorithm::Direct;
}
void HDivFieldEvaluator::Bounds(double lower[3], double upper[3]) const {
    for (int k = 0; k < 3; ++k) { lower[k] = m_impl->lower[k]; upper[k] = m_impl->upper[k]; }
}

} // namespace rad_hdiv
