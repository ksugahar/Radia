# Explicit Chebyshev mass inverse action

Research option `mass-chebyshev`; existing defaults unchanged. The source
implements eight fixed left-Jacobi-preconditioned Chebyshev steps from zero.
No data-dependent early exit occurs inside this inverse action. The outer
PCG original-equation residual and nonlinear residual gates remain unchanged.

Let D = invChi diag(M), S = D^(-1/2) invChi M D^(-1/2), and U be the
infinity norm of the similar matrix D^(-1) invChi M, inflated by 1e-12.
For SPD M, 0 < lambda(S) <= U. Set L=U/10, theta=(U+L)/2,
delta=(U-L)/2. L is a smoothing target, not a certified minimum eigenvalue.
The degree-eight residual polynomial is

    p(lambda)=T8((theta-lambda)/delta)/T8(theta/delta).

For lambda in [L,U], |p|<1. For 0<lambda<L, the positive numerator decreases
from the denominator to one, so 0<p<1 as well. Thus q=(1-p)/lambda>0
throughout the actual positive spectrum. D^(-1/2) q(S) D^(-1/2) is SPD.
The recurrence is fixed and independent of the RHS, preserving a linear
symmetric inverse action for PCG. Floating-point and original residual gates
still matter; an indefinite material is rejected, never repaired or replaced.
One extra length-n FP64 direction buffer is included in native accounting.

Established context: [PETSc Chebyshev documentation](https://petsc.org/main/manualpages/KSP/KSPCHEBYSHEV/)
and [eigenvalue interval guidance](https://petsc.org/main/manualpages/KSP/KSPChebyshevSetEigenvalues/).
This implementation is independently written; it does not link PETSc.

Required acceptance: scalar residual-polynomial check, identical frozen
operator/RHS positive/negative/zero solves and invalid-mass recovery, independent
CPU nonlinear comparison on the RTX, then same artifacts on the compute GPU.
Performance requires fresh transfer-inclusive whole and prepared timings.
Compilation or reduced iteration counts alone do not establish a speedup.
