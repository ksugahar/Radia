# Reduced eddy-current models

Radia uses Foster modal models, PRIMA projection, and snapshot POD for
electromagnetic model reduction. Retired ladder construction guides are not
supported implementation recipes.

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

The descriptor branch projects the MNA form `G = [[R, A^T], [-A, 0]]`,
whose symmetric part is positive semidefinite, so every reduced pole is
stable. Projecting the symmetric saddle `[[R, A^T], [A, 0]]` on the same
Krylov basis gives the same in-band impedance but left unstable
out-of-band poles.

### Circuit and Simulink export

`radia.prima_export` turns a PRIMA result into a descriptor port model
`E x' = A x + B u`, `y = C x + D u` (`from_prima_hacapk`, or
`DescriptorPortModel` from explicit matrices). The input is the port
current (impedance) or the port voltage (admittance). `inverse()` switches
orientation exactly.

- **LTspice:** `write_ltspice_subckt` writes a `.subckt` with one
  `Laplace=` E or G source per real first- or second-order pole section. It
  uses the orientation whose transfer function is proper. The element fixes
  the port relation, so either drive works. The pole-residue form is
  checked against the descriptor response (relative 1e-6). Unstable or
  defective poles, and a model that is improper in both orientations, are
  rejected. AC analysis reproduces the model to round-off. A transient run
  inherits LTspice's numerical Laplace convolution: a coupled RL step
  response was within 1e-3.
- **Simulink:** `export_prima_lti_json` writes `radia.prima.port_model.v1`.
  `radia.simulink.loadPrimaPortModel` builds the `dss` object and refuses
  the file unless `freqresp` reproduces the stored check response.
  `radia.simulink.addPrimaLTIBlock` places it in the Control System
  Toolbox LTI System block. That block needs a proper model, so export an
  inductive port as an admittance.

Tests: `tests/test_prima_export.py`,
`tests/ltspice/test_prima_laplace_ltspice.py` and
`tests/matlab/test_prima_port_model.m`.

## Modal-bulk SIBC

The eigenmode-bulk and surface-impedance composition lives under the legacy
module path `radia.maglev.mixed_galerkin`. Its public method name is
**Modal-bulk SIBC**. Geometry, frequency range, modal truncation, and
surface corrections need their own validation; historical reduced-model comparisons
do not establish acceptance of a current implementation.

Keep full-model, reduced-model, and force/thermal validation separate.
Agreement in impedance alone does not establish force or loss accuracy.
