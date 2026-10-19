# Nonlinear HDiv material scaling boundary

The accepted 288-DoF tests do not certify larger nonlinear systems. The 2412-DoF public sphere has 39102 material points and 3519180 mapped shape entries. Its conservative current MATLAB material bound is 537830768 bytes; explicit budgets are required. The independent CPU case at H=(0,0,100000) A/m takes 50 Newton iterations at nonlinear tolerance 1e-10, so ten iterations must not be advertised as typical.

The 6594-DoF sphere has 112161 material points and 10094490 potential shape entries, exceeding the existing native 65536-point and 4000000-entry limits. These limits remain unchanged. A capacity change needs a chunked native interface, not an unchecked increase.

Proposed implementation, not yet built or accepted:

- Select an explicit bounded element interval; let NGSolve own IntegrationRuleSpace numbering, mapped shapes, orientation, and weights. Return checked global sample indices rather than assuming element sample numbering is contiguous.
- Limit every native map to 65536 points and 4000000 potential shape entries. Bound integer products and sparse allocation before construction. Compare each chunk with the corresponding independent CPU samples.
- Accumulate material energy, gradient, and consistent symmetric tangent on the GPU using the same quadrature. Gather only the assembled gradient/tangent at the existing solver boundary. Keep the geometry H-matrix through all Newton iterations.
- Bound sparse tangent structure from element connectivity; the existing nine-times-dense quadratic bound prevents scaling and cannot simply be removed. Record sparse temporary bounds separately from unknown CUDA-library workspace and sampled available memory.
- Sequence H-matrix construction and material-map upload so transient geometry potentials do not overlap all material maps unnecessarily. Measure actual owned lifetimes; do not add independent peaks or claim a global peak from buffer counters.
- First test chunk versus unchunked energy/gradient/tangent on 288 DoFs, then independent CPU full solves at 2412, 6594 and above 10000 DoFs. Preserve fixed-source nonlinear residual, original-equation inner residual and Armijo energy acceptance. Use a declared SPD sparse preconditioner before large dense factorization becomes a memory barrier.

No numerical threshold relaxation, geometry substitution, CPU fallback, or large-case speed claim follows from this plan. Transfer-inclusive cold/warm measurements and same-artifact acceptance remain separate gates.

The first 2412-DoF comparison against CPU production defaults failed: coefficient relative difference 2.0015e-4 and independently reconstructed nonlinear residual 5.2321e-6. Decomposition found RHS error 1.7443e-16, material-gradient error 1.7644e-14 relative to source, and demagnetizing-action error 5.2321e-6 relative to source. Existing CPU `ho_far_factor=inf` produces an independent residual 9.8025e-11 for the same GPU coefficients, consistent with the GPU internal 9.3226e-11. The GPU curved path uses the all-high-order integration profile, while CPU defaults use far_quad=3/ho_far_factor=2.0. This setting difference matters for nonlinear acceptance. Record the explicit CPU integration profile; do not relax parity tolerances or call it parity with CPU default split integration. A fresh nonlinear all-high-order CPU fixture is required.
