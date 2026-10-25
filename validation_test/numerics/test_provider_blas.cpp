// Provider-neutral LP64/complex ABI analytic checks; no Radia output oracle.
#include "rad_blas.h"
// Fortran declarations stay private to TUs without NGSolve la.hpp/clapack.h.
#if defined(RADIA_BLAS_MKL)
#include <mkl_lapack.h>
#elif defined(RADIA_BLAS_NVPL)
#include <nvpl_lapack.h>
#endif
#include <cmath>
#include <complex>
#include <limits>
#include <cstdio>
#include <stdexcept>
int main() {
    double a[]={1,2,3,4,5,6}, b[]={1,2,3,4,5,6}, c[4]={0};
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,2,2,3,1,a,3,b,2,0,c,2);
    const double expected[]={22,28,49,64};
    for(int i=0;i<4;++i) if(c[i]!=expected[i]) throw std::runtime_error("CBLAS row-major GEMM");
    // Complex 1x1 solve with nonzero imaginary parts exercises pair layout and
    // Fortran LAPACK symbol/calling convention, not merely a real-only cast.
    radia_blas_complex_double z={2,1},rhs={1,3};
    int n=1,nrhs=1,pivot=0,info=-1;
    zgesv_(&n,&nrhs,&z,&n,&pivot,&rhs,&n,&info);
    if(info!=0 || std::abs(rhs.real-1.)>8*std::numeric_limits<double>::epsilon() ||
       std::abs(rhs.imag-1.)>8*std::numeric_limits<double>::epsilon())
        throw std::runtime_error("LP64 complex Fortran ABI");
    double matrix[]={2,1,1,3}, y[]={4,7};
    lapack_int ipiv[2];
    if(LAPACKE_dgesv(LAPACK_ROW_MAJOR,2,1,matrix,2,ipiv,y,1)!=0 ||
       std::abs(y[0]-1)>8*std::numeric_limits<double>::epsilon() ||
       std::abs(y[1]-2)>8*std::numeric_limits<double>::epsilon())
        throw std::runtime_error("LAPACKE pivot/layout ABI");
    std::printf("{\"analytic_blas\":\"PASS\",\"provider\":\"%s\"}\n",RADIA_BLAS_NAME);
}
