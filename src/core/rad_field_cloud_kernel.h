/* Shared CPU/CUDA point-source and quadrupole field kernels. No FE semantics,
 * scheduling, approximation substitution, or 1/(4*pi) scaling lives here. */
#pragma once
#include <cmath>
#include <cstddef>
#include <vector>
#include "rad_field_sources.h"

#ifdef __CUDACC__
#define RAD_FIELD_HD __host__ __device__
#else
#define RAD_FIELD_HD
#endif

namespace rad_hdiv {
struct FieldCloudPoint { double position[3]; double strength; };
struct FieldSourceView {
    const FieldTetSource* tets;
    const FieldTriSource* triangles;
    const FieldCubicTriSource* cubic_triangles;
    const FieldCurvedTetSource* curved_tets;
    const FieldCurvedTriSource* curved_triangles;
    const FieldCloudPoint* points;
    const double* gauss_points;
    const double* gauss_weights;
    int nq;
};
struct FieldCloudNode {
    std::size_t begin, end;
    int left, escape;
    double center[3], radius, charge, dipole[3], second[9];
};
struct FieldCloudImage { int mask; double sign, angle; };
struct FieldCloudSnapshot {
    std::vector<FieldCloudPoint> points;
    std::vector<FieldSourceRef> sources;
    std::vector<FieldTetSource> tets;
    std::vector<FieldTriSource> triangles;
    std::vector<FieldCubicTriSource> cubic_triangles;
    std::vector<FieldCurvedTetSource> curved_tets;
    std::vector<FieldCurvedTriSource> curved_triangles;
    std::vector<double> gauss_points, gauss_weights;
    std::vector<std::size_t> order;
    std::vector<FieldCloudNode> nodes;
    std::vector<FieldCloudImage> images;
    double theta;
};

struct FieldCompensated3 {
    double sum[3]{}, correction[3]{};
    RAD_FIELD_HD void Add(int k, double x) {
        const double next = sum[k]+x;
        correction[k] += std::fabs(sum[k]) >= std::fabs(x)
            ? (sum[k]-next)+x : (x-next)+sum[k];
        sum[k] = next;
    }
    RAD_FIELD_HD void Add(const double v[3]) { for(int k=0;k<3;++k) Add(k,v[k]); }
    RAD_FIELD_HD void Store(double v[3]) const {
        for(int k=0;k<3;++k) v[k]=sum[k]+correction[k];
    }
};

RAD_FIELD_HD inline void AddCloudPoint(const double position[3], double strength,
                                      const double r[3], double out[3]) {
    const double dx=r[0]-position[0], dy=r[1]-position[1], dz=r[2]-position[2];
    const double r2=dx*dx+dy*dy+dz*dz;
    if(r2<=1e-300) return; // Preserve the existing retained-cloud convention.
    const double scale=strength/(r2*std::sqrt(r2));
    out[0]+=scale*dx; out[1]+=scale*dy; out[2]+=scale*dz;
}

RAD_FIELD_HD inline void AddFieldMultipole(const double center[3], double charge,
        const double dipole[3], const double second[9], const double r[3], double out[3]) {
    const double R[3]={r[0]-center[0],r[1]-center[1],r[2]-center[2]};
    const double r2=R[0]*R[0]+R[1]*R[1]+R[2]*R[2];
    if(r2<=1e-300) return;
    const double inv_r=1.0/std::sqrt(r2), inv_r3=inv_r/r2;
    const double inv_r5=inv_r3/r2, inv_r7=inv_r5/r2;
    const double pR=dipole[0]*R[0]+dipole[1]*R[1]+dipole[2]*R[2];
    double second_R[3]{}, trace=0;
    for(int i=0;i<3;++i) {
        trace+=second[3*i+i];
        for(int j=0;j<3;++j) second_R[i]+=second[3*i+j]*R[j];
    }
    const double R_second_R=R[0]*second_R[0]+R[1]*second_R[1]+R[2]*second_R[2];
    for(int k=0;k<3;++k) {
        out[k]+=charge*R[k]*inv_r3;
        out[k]+=3.0*R[k]*pR*inv_r5-dipole[k]*inv_r3;
        out[k]+=7.5*R[k]*R_second_R*inv_r7
                -1.5*(trace*R[k]+2.0*second_R[k])*inv_r5;
    }
}
} // namespace rad_hdiv
#undef RAD_FIELD_HD
