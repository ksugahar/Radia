# Reduced eddy-current models

Radia uses Foster modal models, PRIMA projection, and snapshot POD for
electromagnetic model reduction. Retired CLN/Cauer construction guides are
available in Git history; they are not supported implementation recipes.

## Foster modal dynamics

`radia.vim.HCurlEddyFosterModel` diagonalizes the reduced pair `R, L` using
`R v = lambda L v` and `V^T L V = I`. Its harmonic modal coordinates satisfy
`z_k = (V^T P u)_k / (s + lambda_k)`. This preserves the reduced direct
frequency-domain solution when all modes are retained. Truncation requires
an error check over the intended operating band.

`radia.maglev.position_foster.MovingHCurlFosterFamily` shares this operator
between positions and interpolates the source port matrix. Models with
different reduced operators cannot be combined as this family.

A nonzero surface-impedance block is generally non-rational in frequency.
The Foster implementation retains the direct harmonic solve for that block
and rejects unsupported state-space export. It does not replace the block
with an undocumented finite circuit approximation.

## PRIMA projection

`radia.prima_hacapk.PRIMAHACApKModel` projects the PEEC system `R + s L`
onto a Krylov basis. HACApK accelerates inductance matrix-vector products;
the model evaluates the projected matrices and port vector. Choose the
reduced dimension from comparison with the full model over the required
frequency band. A continued-fraction construction is not this projection.

## Modal-bulk SIBC

The eigenmode-bulk and surface-impedance composition lives under the legacy
module path `radia.maglev.mixed_galerkin`. Its public method name is
**Modal-bulk SIBC**. Geometry, frequency range, modal truncation, and
surface corrections need their own validation; historical CLN comparisons
do not establish acceptance of a current implementation.

Keep full-model, reduced-model, and force/thermal validation separate.
Agreement in impedance alone does not establish force or loss accuracy.
