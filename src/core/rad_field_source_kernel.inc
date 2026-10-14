// Shared exact source evaluation. Include after rad_hdiv_vim.h (CPU) or
// rad_hdiv_field_math.inc (CUDA); integration and summation order are preserved.
namespace RAD_FIELD_NAMESPACE {
using rad_hdiv::AddCloudPoint;
using TetSource=rad_hdiv::FieldTetSource;
using TriSource=rad_hdiv::FieldTriSource;
using CubicTriSource=rad_hdiv::FieldCubicTriSource;
using CurvedTetSource=rad_hdiv::FieldCurvedTetSource;
using CurvedTriSource=rad_hdiv::FieldCurvedTriSource;
using SourceKind=rad_hdiv::FieldSourceKind;
template <std::size_t N>
RAD_FIELD_KERNEL inline void Translate(const double (&v)[N][3], const double origin[3], const double r[3],
               double (&local_v)[N][3], double local_r[3]) {
    for (std::size_t i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k) local_v[i][k] = v[i][k] - origin[k];
    for (int k = 0; k < 3; ++k) local_r[k] = r[k] - origin[k];
}

RAD_FIELD_KERNEL inline double CurvedTetDensity(const CurvedTetSource& source, double xi, double eta, double zeta) {
    return source.coefficient[0] + source.coefficient[1]*xi
         + source.coefficient[2]*eta + source.coefficient[3]*zeta;
}

RAD_FIELD_KERNEL inline double CurvedTriDensity(const CurvedTriSource& source, double xi, double eta) {
    return source.coefficient[0] + source.coefficient[1]*eta
         + source.coefficient[2]*eta*eta + source.coefficient[3]*xi
         + source.coefficient[4]*xi*eta + source.coefficient[5]*xi*xi;
}

RAD_FIELD_KERNEL inline void AddCurvedTet(const CurvedTetSource& source, const double r[3], double out[3], const double* gauss_points,
    const double* gauss_weights, int nq) {
        double xi0[3];
        ClosestRefTet(source.nodes, r, xi0);
        const double corners[4][3] = {{0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        const int faces[4][3] = {{1,2,3},{0,3,2},{0,1,3},{2,1,0}};
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

RAD_FIELD_KERNEL inline void AddCurvedTriangle(const CurvedTriSource& source, const double r[3], double out[3], const double* gauss_points,
    const double* gauss_weights, int nq) {
        double xi0[2];
        ClosestRefTri(source.nodes, r, xi0);
        const double corners[3][2] = {{0,0},{1,0},{0,1}};
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

RAD_FIELD_KERNEL inline void AddFieldSource(const rad_hdiv::FieldSourceView& view,
    const rad_hdiv::FieldSourceRef& atom, const double r[3], double out[3]) {
        if (atom.kind == SourceKind::Tet) {
            const TetSource& source = view.tets[atom.index];
            double v[4][3], local_r[3];
            Translate(source.v, source.origin, r, v, local_r);
            double value[3]{};
            if (source.degree <= 3) {
                TetVolFieldCubic(v, local_r, source.coefficient, value);
            } else {
                double basis[84][3];
                TetVolFieldBasisUpTo6(v, local_r, basis);
                for (int index = 0; index < 84; ++index)
                    for (int k = 0; k < 3; ++k)
                        value[k] += source.coefficient[index]*basis[index][k];
            }
            for (int k = 0; k < 3; ++k) out[k] += value[k];
        } else if (atom.kind == SourceKind::Triangle) {
            const TriSource& source = view.triangles[atom.index];
            double v[3][3], local_r[3];
            Translate(source.v, source.origin, r, v, local_r);
            double value[3];
            QuadTriField(v, local_r, source.sigma0, source.slope, source.hessian, value);
            for (int k = 0; k < 3; ++k) out[k] += value[k];
        } else if (atom.kind == SourceKind::CubicTriangle) {
            const CubicTriSource& source = view.cubic_triangles[atom.index];
            double v[3][3], local_r[3];
            Translate(source.v, source.origin, r, v, local_r);
            double value[3]{};
            if (source.degree <= 3) {
                CubicTriField(v, local_r, source.coefficient, value);
            } else {
                double basis[35][3];
                TriFieldBasisUpTo4(v, local_r, basis);
                for (int index = 0; index < 35; ++index)
                    for (int k = 0; k < 3; ++k)
                        value[k] += source.coefficient[index]*basis[index][k];
            }
            for (int k = 0; k < 3; ++k) out[k] += value[k];
        } else if (atom.kind == SourceKind::CurvedTet) {
            AddCurvedTet(view.curved_tets[atom.index], r, out, view.gauss_points, view.gauss_weights, view.nq);
        } else if (atom.kind == SourceKind::CurvedTriangle) {
            AddCurvedTriangle(view.curved_triangles[atom.index], r, out, view.gauss_points, view.gauss_weights, view.nq);
        } else {
            const rad_hdiv::FieldCloudPoint& source = view.points[atom.index];
            AddCloudPoint(source.position, source.strength, r, out);
        }
    
}
} // namespace RAD_FIELD_NAMESPACE
