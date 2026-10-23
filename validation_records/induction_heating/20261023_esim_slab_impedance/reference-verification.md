# Analytic verification of the adopted slab reference

This is an algebraic derivation independent of the discrete production solver.
Its companion verify_reference.py checks the identities with exact symbolic
algebra; symbolic-verification.json binds the executed script and this text by
SHA256. Scientific review of assumptions and applicability remains required.
The source of the diffusion equation is Haus and Melcher, *Electromagnetic
Fields and Energy* (1989), ISBN 9780132490207, section 10.7:
https://web.mit.edu/6.013_book/www/chapter10/10.7.html .
The canonical bibliography key is haus1989electromagnetic.
That section supplies the differential equation and exponential solutions;
the symmetric boundary conditions below specify this adopted case, rather than
claiming its worked example has the same boundary conditions. A DOI for this
book has not been verified. No DOI belonging to a review or unrelated paper is
substituted for the actual reference.

## Equation, boundary conditions and uniqueness

Let x run from the slab centre (0) to its surface (a), with constant positive
real conductivity sigma and permeability mu. Phasors use exp(+j omega t).
Combining Ampere, Ohm and Faraday gives H'' = j omega mu sigma H.
Choose k with k^2 = j omega mu sigma and Re(k)>0. Define
H(x) = H0 cosh(k x)/cosh(k a).
Direct differentiation gives H'' = k^2 H, H'(0)=0 and H(a)=H0.
In surface-depth coordinates z=a-x the same expression is
H(z)=H0 cosh(k(a-z))/cosh(k a).
For positive omega, mu and sigma, conjugating H changes j to -j; the conjugate
does not satisfy the declared phasor equation. It must not be an accepted oracle.
For a difference u between two solutions with homogeneous boundary data,
integration by parts gives integral |u'|^2 + j omega mu sigma |u|^2 = 0.
Its real part makes u constant; u(a)=0 then gives uniqueness.

## Surface impedance and complex power

With consistent tangential orientation, J=H' and E=rho H', rho=1/sigma.
Zs=E(a)/H(a)=rho k tanh(k a).
Multiplying the differential equation by rho H* and integrating gives
rho H'(a) H*(a) = integral_0^a [rho |H'|^2 + j omega mu |H|^2] dx,
because H'(0)=0. Division by |H0|^2 identifies Zs with the complex power
integral normalization. The common time-average factor 1/2 cancels.
Its real part is the nonnegative Joule loss and its imaginary part is the
positive inductive term for exp(+j omega t). This verifies the reference used
for the impedance criterion independently of the discrete solver.

## Phase quantity, convergence and limits of inference

The phase criterion is sqrt(integral |H_FE-H_exact|^2 /
integral |H_exact|^2), where H_FE is the piecewise-linear FE interpolant,
integrated element by element with 32-point Gauss-Legendre quadrature.
The existing measured 16/32-point comparison is in the original results;
this text does not introduce a new quadrature or validation run.
Smooth constant-coefficient fields and uniformly refined smooth grading give
an asymptotic second-order expectation for this P1 approximation. The approved
q>=1.5 is an empirical acceptance lower bound, not a universal theorem for
surface impedance or an unstructured adaptive mesh. Rates use actual h_max/a.
Here a/delta is about 105. h_max lies deep in negligible field, so h_min/delta
is also recorded and halves under this fixed-grid refinement. These data
validate deep-skin diffusion, not finite-thickness or midplane accuracy.
