// Device declarations matching the pure scalar API in rad_hdiv_vim.h.
RAD_FIELD_KERNEL double TriPotential(const double V[3][3], const double r[3]);
RAD_FIELD_KERNEL double PhiTet(const double V[4][3], const double P[3]);
RAD_FIELD_KERNEL void TriField(const double V[3][3], const double r[3], double out[3]);
RAD_FIELD_KERNEL void TetField(const double V[4][3], const double P[3], double out[3]);
RAD_FIELD_KERNEL void TriMoment1(const double V[3][3], const double r[3], double out[3]);
RAD_FIELD_KERNEL void TriMoment2(const double V[3][3], const double r[3], double out[3][3]);
RAD_FIELD_KERNEL void TetMoment1(const double V[4][3], const double r[3], double out[3]);
RAD_FIELD_KERNEL void TetPotentialMomentsDirectionalUpTo1(
    const double V[4][3], const double dV[4][3], const double r[3], const double dr[3],
    double value[4], double direction[4]);
RAD_FIELD_KERNEL void TetPotentialMomentsUpTo3(const double V[4][3], const double r[3],
                              double out[20]);
RAD_FIELD_KERNEL void TetPotentialMomentsUpTo6(const double V[4][3], const double r[3],
                              double out[84]);
RAD_FIELD_KERNEL void TetPotentialMomentsDirectionalUpTo3(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[20],double direction[20]);
RAD_FIELD_KERNEL void TetPotentialMomentsDirectionalUpTo6(const double V[4][3],const double dV[4][3],
    const double r[3],const double dr[3],double value[84],double direction[84]);
RAD_FIELD_KERNEL void TriPotentialMomentsUpTo4(const double V[3][3], const double r[3],
                              double out[35]);
RAD_FIELD_KERNEL void TriPotentialMomentsUpTo2(const double V[3][3], const double r[3],
                              double out[10]);
RAD_FIELD_KERNEL void TriPotentialMomentsDirectionalUpTo4(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[35],double direction[35]);
RAD_FIELD_KERNEL void TriPotentialMomentsDirectionalUpTo2(
    const double V[3][3],const double dV[3][3],const double r[3],const double dr[3],
    double value[10],double direction[10]);
RAD_FIELD_KERNEL void TetVolFieldLinear(const double V[4][3], const double r[3], double rho0,
                       const double g[3], double out[3]);
RAD_FIELD_KERNEL void TetVolFieldQuadratic(const double V[4][3], const double r[3], double rho0,
                          const double g[3], const double Q[3][3], double out[3]);
RAD_FIELD_KERNEL void TetVolFieldCubic(const double V[4][3], const double r[3],
                      const double coefficient[20], double out[3]);
RAD_FIELD_KERNEL void TetVolFieldCubicBasis(const double V[4][3], const double r[3],
                           double out[20][3]);
RAD_FIELD_KERNEL void TetVolFieldBasisUpTo6(const double V[4][3], const double r[3],
                           double out[84][3]);
RAD_FIELD_KERNEL void LinTriField(const double V[3][3], const double r[3], double sigma0,
                 const double s[3], double out[3]);
RAD_FIELD_KERNEL void QuadTriField(const double V[3][3], const double r[3], double sigma0,
                  const double s[3], const double S[3][3], double out[3]);
RAD_FIELD_KERNEL void QuadTriFieldBasis(const double V[3][3], const double r[3],
                       double out[10][3]);
RAD_FIELD_KERNEL void TriFieldBasisUpTo4(const double V[3][3], const double r[3],
                        double out[35][3]);
RAD_FIELD_KERNEL void CubicTriField(const double V[3][3], const double r[3],
                   const double coefficient[20], double out[3]);
RAD_FIELD_KERNEL void TetVolFieldLinearDirectional(
    const double V[4][3], const double dV[4][3],
    const double r[3], const double dr[3],
    double rho0, double drho0, const double g[3], const double dg[3],
    double value[3], double direction[3]);
RAD_FIELD_KERNEL void QuadTriFieldDirectional(
    const double V[3][3], const double dV[3][3],
    const double r[3], const double dr[3],
    double sigma0, double dsigma0,
    const double s[3], const double ds[3],
    const double S[3][3], const double dS[3][3],
    double value[3], double direction[3]);
RAD_FIELD_KERNEL void ClosestPointTriangle(const double p[3], const double a[3], const double b[3], const double c[3],
                          double out[3]);
RAD_FIELD_KERNEL void ClosestPointTet(const double V[4][3], const double p[3], double out[3]);
RAD_FIELD_KERNEL void ClosestRefTri(const double nodes[6][3], const double p[3], double xi0[2]);
RAD_FIELD_KERNEL double CurvedTriPotential(const double nodes[6][3], int e0, int e1, const double p[3],
                          const double* gl, const double* gw, int nq,
                          bool include_measure = true);
RAD_FIELD_KERNEL void ClosestRefTet(const double nodes[10][3], const double p[3], double xi0[3]);
RAD_FIELD_KERNEL double CurvedTetPotential(const double nodes[10][3], int e0, int e1, int e2, const double p[3],
                          const double* gl, const double* gw, int nq,
                          bool include_measure = true);
RAD_FIELD_KERNEL void CurvedTriMapMeasure(const double nodes[6][3], double xi, double eta, double X[3], double& dA);
RAD_FIELD_KERNEL void CurvedTetMapMeasure(const double nodes[10][3], double xi, double eta, double zeta, double X[3], double& dV);
