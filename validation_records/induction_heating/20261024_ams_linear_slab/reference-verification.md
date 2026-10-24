# Analytic reference and 3D embedding (manual derivation)

Reference: canonical bibliography `haus1989electromagnetic`, Haus and Melcher,
Electromagnetic Fields and Energy (1989), ISBN-13 9780132490207, section 10.7.
The slab symmetry specialization and this vector-potential embedding must be
reviewed; this document is not an executed computer-algebra verification.

Let k^2 = i omega mu sigma, k = (1+i)/delta, with constant positive mu and sigma.
For -a <= x <= a,

    A = (0, mu H0 sinh(k x)/(k cosh(k a)), 0),
    H = curl(A)/mu = (0, 0, H0 cosh(k x)/cosh(k a)),
    B = mu H, E = -i omega A, J = sigma E.

Since A_y'' = k^2 A_y, curl(mu^-1 curl(A)) + i omega sigma A = 0.
The tangent trace of this exact A is prescribed on all six extrusion faces.
This is a driven extrusion, not an isolated coil. The field quantities use
relative complex volume L2 errors, not magnitude-only comparisons.

At x = s a, s = +/-1, -s E_y/H_z = (k/sigma) tanh(k a). The outward sign and
volume-to-boundary traces must be retained. Mean Joule power is
(1/(2 sigma)) integral |J|^2 dV = Ly Lz H0^2 Re[(k/sigma) tanh(k a)], including
both depth faces, for peak phasors. A same-matrix direct solve checks solver
parity, not independent formulation accuracy; the analytic reference is the
independent formulation check. H and B are proportional, not two independent
detectors in this constant-mu case.
