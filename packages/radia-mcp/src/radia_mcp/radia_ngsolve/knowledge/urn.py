r"""Universal Relaxation Network (URN) -- causal/passive rational fitting of a

SHOWCASE NOTEBOOK: docs/universal_relaxation_network/urn_showcase.ipynb -- legacy validation figures (URN-vs-VF, NASA/TDK fits, older attention ablation).
frequency response, with direct time-domain (relaxation-network / SPICE / ADE)
synthesis.

URN is a KAN-inspired network that decomposes a measured or computed frequency
response Z(omega) (impedance, dispersive permittivity/permeability, or an
open-boundary DtN symbol G_n(omega)) into a sparse sum of PHYSICAL relaxation
mechanisms -- Debye, Cole-Cole, Cole-Davidson, Havriliak-Negami, CPE, Warburg,
Gerischer, RLC, skin-effect -- in both a series and a parallel (admittance)
branch, with KAN-style adaptive relaxation-time (tau) refinement.  Because every basis is a passive,
causal relaxation, the fitted model is GUARANTEED causal/passive and maps
directly to a stable time-domain realisation: an equivalent SPICE circuit
(generate_spice_netlist) or, equivalently, a set of first-order auxiliary
differential equations (one per relaxation pole).

This is the right tool to turn a frequency-domain absorbing-BC / dispersive
material response into a TIME-DOMAIN representation for FETD / Newmark-beta:
fit G_n(omega) (or eps(omega), mu(omega)) with URN, then each relaxation
mechanism becomes a local-in-time auxiliary ODE -- no history convolution,
and passivity (hence time-stepping stability) is built in.  Naively dropping a
complex eps/mu sampled at one frequency into the time domain is non-causal;
URN provides the causal broadband surrogate.

Reference: K. Sugahara and Y. Sato, "KAN-inspired Universal Relaxation Network
for Automatic Discovery of Physical Relaxation Mechanisms with Direct Circuit
Synthesis," IEEE Access, 2026.  Canonical runtime implementation:
src/radia/urn; paper assets and benchmarks:
docs/universal_relaxation_network/.

API (this module):
  get_urn_documentation(topic) -- knowledge text (this module)
  run_urn_fit(freqs, Z, ...)   -- fit and return mechanisms + NRMSE + SPICE
  urn_fit_from_csv(path, ...)  -- read a (freq, Re, Im) CSV, fit, return a report
"""

import os

URN_OVERVIEW = r"""
# Universal Relaxation Network (URN)

URN fits a complex frequency response  Z(omega)  by a SPARSE sum of physical
relaxation basis functions and synthesises a causal time-domain (circuit / ADE)
realisation directly from the fit.

Why it exists (vs Vector Fitting):
- Vector Fitting (VF) places arbitrary rational poles; passivity must be
  enforced post-hoc and can fail for fractional-order (Cole-Cole/CPE) data.
- URN's bases ARE passive relaxations, so the model is causal/passive BY
  CONSTRUCTION, and each term is physically interpretable (a tau, an exponent).
- On real data (NASA 18650 battery EIS, TDK MnZn power ferrites PC47/50/95/200)
  URN matches or beats VF -- average ~22.8% lower NRMSE, with the largest gains
  (39-66%) exactly where Cole-Cole / fractional dynamics dominate; VF wins only
  on near-ideal single-Debye behaviour.

Pipeline:
  freq response Z(omega)  --train_urn-->  sparse relaxation model
       --get_active_components-->  {tau, alpha, beta, weight} per mechanism
       --generate_spice_netlist-->  equivalent circuit (== auxiliary-ODE set)
"""

URN_METHOD = r"""
# URN method

Model:  Z(omega) = Z_inf + sum_k w_k * basis_k(omega; tau_k, ...)   (series)
                          + parallel/admittance branch (Y-space)
The current SA/RM direction avoids frequency-dependent attention and instead
uses continued-fraction residual peeling to obtain frequency selectivity from
the circuit topology.

Relaxation basis library (relaxation_basis_library.py):
  Debye:             1 / (1 + j w tau)
  Cole-Cole:         1 / (1 + (j w tau)^alpha)            0<alpha<=1
  Cole-Davidson:     1 / (1 + j w tau)^beta
  Havriliak-Negami:  1 / (1 + (j w tau)^alpha)^beta       (general fractional)
  CPE / Warburg:     fractional / diffusion (alpha=1/2) elements
  Gerischer, RLC, skin-effect (Dowell RL ladder)
  + parallel (admittance-space) Debye / Cole-Cole / CPE / Warburg
Magnetic variants (permeability): magnetic_debye, magnetic_cole_cole,
  two_relaxation_permeability.

Training (train_urn): Adam + cosine LR, multi-restart, relative-error loss
  + sparsity penalty (drives unused bases to ~0 weight so the active set is
  the discovered mechanism count).  KAN-style AdaptiveURN refines tau on a
  coarse->fine grid where local error is high.

log-spaced tau parameters are stored as log_tau (nn.Parameter); active
components are read out by get_active_components(threshold) (weight above
threshold), returning tau (in seconds), alpha/beta, weight_magnitude, branch.

SA/RM-2026 paper variant (YAdmittanceURN): the research-meeting draft by
Sato/Sugahara uses a pure attention-free Y-domain dictionary:
Debye/magnetic-Debye/Cole-Cole/magnetic-Cole-Cole/inductive-CPE/
capacitive-CPE/series-RLC = 22 bases, S-domain Huber loss with z0=median(|Z|),
nonnegative gates, and output-ablation importance I_i.  This is exposed as train_y_admittance_urn,
refit_y_admittance_active_bases, and s_domain_rmse so the draft can evolve
without replacing the original NASA/TDK train_urn path.

Cauer-ladder review variant (CauerLadderURN): instead of increasing the number
of parallel bases, a small passive continued fraction is fitted:
  Z_k = R_k + s L_k + 1 / (G_k + s C_k + 1 / Z_{k+1}).
A 6-section ladder has 24 positive parameters, close to the 22-basis draft
dictionary, but represents pole-zero/anti-resonance behavior through
series/parallel nesting.  train_cauer_ladder_alternating updates R,L with a
direct impedance-domain residual and G,C with a direct admittance-domain
residual, rolling back unstable blocks and lowering the learning rate.  This is
attention-free and therefore a better candidate for direct time-domain circuit
synthesis.  use_peeling_initialization=True is experimental; naive peeling can
over-subtract and should not yet be treated as the default result path.
For resonance/anti-resonance data, set use_rational_initialization=True and
use_least_squares_polish=True: Radia first builds a small pole-zero rational
teacher, distills it into positive Cauer parameters by nonlinear least squares,
then optionally uses train_cauer_ladder_tail_then_polish to freeze outer
sections while fitting the inner ladder and finally unfreeze all sections.

CLN peeling variant (paired-basis): train_cln_peeling_urn fits, per stage,
  R_n = Z_2n + (Z_2n+1 || R_n+1)
where one 22-basis composite fit supplies the physical basis shapes and its
coefficients are split continuously between the even series branch and the odd
shunt branch (soft split a_2n,k = a_k p_k, a_2n+1,k = a_k (1-p_k)) while a
fresh 22-basis lookahead model represents R_n+1.  Accepted pairs are frozen
(no global polish; past stages are never re-trained) and the measured tail is
peeled by the exact inverse map R_n+1 = [1/(R_n - Z_2n) - 1/Z_2n+1]^-1.
Evaluation policy: report the termination="lookahead" S-domain RMSE (learned
physical tail) plus audit_passivity() -- a dense-grid min Re(Z)/min Re(Y)
audit of every frozen branch and the full lookahead ladder, with 1-decade
extrapolation.  termination="stored" re-inserts the exactly peeled measurement
residue: identity reconstruction only, never a fit-accuracy number.
Stage-wise trust region (2026-07-27): the exact peel amplifies measurement
error where R_n - Z_2n cancels, and yields sign-unstable spikes where the
peeled tail admittance nearly vanishes (parallel-resonance bands where the
tail barely loads the ladder).  On the SA/RM PCB coil the 1.36-1.49 MHz
self-resonance band gave a -59.7 kOhm negative-real spike in R_1 (median 56
Ohm), rejecting stage 2 with min_parallel_real_normalized = -1063 regardless
of the series branch.  Each stage therefore stores per-frequency trust
weights (config: denominator_margin_relative, tail_admittance_margin_relative)
computed from the relative series cancellation and the peeled-tail admittance
magnitude; the next stage fits with those weights and is accepted on trusted
points only (min_parallel_real_trusted, min_tail_admittance_real_trusted,
seed/parent_s_rmse_trusted, trusted_fraction >= min_trusted_fraction), and
inherited weights multiply stage by stage.  The stored exact tail is never
modified.
"""

URN_API = r"""
# URN API (canonical impl: radia.urn)

from radia.urn import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist)

config = URNConfig(n_debye=3, n_cole_cole=2, n_warburg=1,
                   sparsity_weight=0.01, n_epochs=2000, n_restarts=3)
model  = train_urn(freqs_Hz, Z_complex, config, verbose=False)
mech   = model.get_active_components()       # {type: [{tau, alpha, weight_magnitude,...}]}
spice  = generate_spice_netlist(model, "Z")  # SPICE netlist string (RC/RL ladders)

URNConfig fields: n_debye/n_cole_cole/n_cole_davidson/n_havriliak_negami/n_cpe/
  n_warburg/n_gerischer/n_rlc/n_skin_effect (series), *_parallel (admittance),
  sparsity_weight, lr, n_epochs(=6000), n_restarts(=10), omega_ref/Z_ref(auto).
  For a responsive tool call, lower n_epochs/n_restarts.

Y-domain SA/RM variant:
from radia.urn import (
    YAdmittanceURNConfig, refit_y_admittance_active_bases,
    s_domain_rmse, train_y_admittance_urn)
cfg = YAdmittanceURNConfig.paper_22_basis()
model = train_y_admittance_urn(freqs_Hz, Z_complex, cfg, verbose=False)
active = model.active_bases()   # output-ablation ranking
Zfit = model.predict(freqs_Hz)
rmse_s = s_domain_rmse(Zfit, Z_complex)
realizable = refit_y_admittance_active_bases(freqs_Hz, Z_complex, active, cfg)

Cauer continued-fraction review variant:
from radia.urn import (
    CauerLadderURNConfig, fit_rational_pole_zero,
    train_cauer_ladder_alternating,
    train_cauer_ladder_progressive)
cfg = CauerLadderURNConfig.twenty_two_parameter_candidate()
model = train_cauer_ladder_alternating(freqs_Hz, Z_complex, cfg, verbose=False)
Zfit = model.predict(freqs_Hz)
sections = model.parameter_summary()  # positive R,L,G,C per Cauer section

Pole-zero assisted Cauer path:
from radia.urn import train_cauer_ladder_tail_then_polish
cfg = CauerLadderURNConfig.twenty_two_parameter_candidate(
    use_rational_initialization=True,
    use_least_squares_polish=True,
    frozen_outer_sections=2)
teacher = fit_rational_pole_zero(freqs_Hz, Z_complex, order=6)
model = train_cauer_ladder_tail_then_polish(freqs_Hz, Z_complex, cfg)

CLN peeling path (paired even/odd split + stage-wise trust region):
from radia.urn import CLNPeelingConfig, train_cln_peeling_urn
cfg = CLNPeelingConfig(
    n_stages=2,
    denominator_margin_relative=3.0e-2,     # series-cancellation trust margin
    tail_admittance_margin_relative=5.0e-3,  # peeled-tail activity trust margin
    min_trusted_fraction=0.8)
model = train_cln_peeling_urn(freqs_Hz, Z_complex, cfg)
Zla = model.predict_terminated(freqs_Hz, termination="lookahead")  # REPORT THIS
audit = model.audit_passivity()  # dense grid + 1-decade extrapolation
# model.predict(freqs_Hz) / termination="stored" = identity check only;
# lookahead/constant/open/short accept arbitrary frequency grids.
# stage.metrics: *_trusted gates, trusted_fraction; stage.tail_trust_weight

This MCP module wraps it:
  run_urn_fit(freqs, Z, n_debye=.., n_cole_cole=.., n_warburg=.., n_epochs=..,
              n_restarts=.., sparsity_weight=.., spice=True) -> dict
      {"nrmse", "mechanisms", "spice_netlist", "n_active"}
  urn_fit_from_csv(path, freq_col, real_col, imag_col, delimiter, skip_rows,
                   ...same fit args..., spice_out="") -> formatted report str
"""

URN_TIMEDOMAIN = r"""
# URN -> time domain (the point for FETD / Newmark-beta)

A URN fit is a sum of passive relaxations, so it has an EXACT local-in-time
realisation.  Each Debye term  w/(1 + j w tau)  <->  an auxiliary variable phi
with    tau dphi/dt + phi = w * u ,   contributing phi to the response -- one
first-order ODE per pole, NO history convolution.  Fractional terms
(Cole-Cole/CPE) are realised as a short RC/RL ladder (Charef / Valsa / Dowell),
i.e. a few extra auxiliary ODEs that approximate the fractional pole over the
band.  generate_spice_netlist emits exactly these ladders; the same ladders are
the ADE/auxiliary-ODE set for an FETD / Newmark-beta solver.

Because the realisation is passive, the resulting time-stepping is STABLE
(no growing modes) -- this is the structural advantage over plugging a raw
complex (eps, mu) or a non-passive rational fit into the time domain, which is
generally non-causal/unstable.

Time-domain verification in the repo: demo_spice_timedomain.py,
run_ltspice_verification.py (actual LTspice), verify_timedomain_stability.py
(URN vs VF stability).
"""

URN_CQ = r"""
# URN -> convolution quadrature (CQ)

CQ is the clean bridge when the solver already has a frequency-domain acoustic,
Maxwell, impedance-boundary, or material operator and we want a causal transient
response without first hand-writing every auxiliary ODE.

Contract:
  1. Fit or identify a passive URN relaxation model from samples of H(j omega).
  2. Expose the fitted model as a Laplace evaluator H(s), not just as a table.
  3. Feed H(delta(zeta)/dt) to Lubich CQ (BDF1/BDF2 are the small teaching
     defaults).
  4. Apply the resulting weights by causal discrete convolution to the boundary
     or material input history.

Why this matters:
  - A plain FFT/IFFT demo treats the sampled time signal as periodic.  For a
    hammer/step-like transient this can create pre-hit wrap-around unless the
    padding/windowing convention is explicit.
  - CQ evaluates the same fitted frequency-domain object through the Laplace
    variable and gives a causal time-domain operator for the chosen time step.
  - Because the URN ladder is non-negative/passive, the CQ kernel inherits a
    stable physical realization; this is the useful path for FEM/BEM acoustics
    and time-domain Maxwell material kernels.

Education artifact:
  docs/universal_relaxation_network/cq_urn_bridge.ipynb
Checked evidence:
  validation_test/universal_relaxation_network/cq_urn_bridge_results.json

The result-bearing notebook fits a two-pole passive relaxation ladder on a
candidate tau grid, builds BDF2 CQ weights, and compares the causal CQ response
with a deliberately naive periodic IFFT contrast.  The point is not that IFFT is
forbidden; the point is that CQ is the native transient operator once URN has
identified a passive H(s).
"""

URN_APPLICATION = r"""
# Using URN for open-boundary / dispersive-material time-domain models

Target a frequency response that you want to march in time:
  - an exterior open-boundary DtN symbol  G_n(omega)  on a truncation circle/
    sphere (per multipole order n), OR
  - a dispersive material response eps(omega) / mu(omega) of an absorbing layer,
    OR a surface impedance / reflection coefficient.

Recipe:
  1. Sample the response over the band of interest: (freq, Re, Im).
  2. run_urn_fit / urn_fit_from_csv -> passive relaxation model + NRMSE.
  3. Each mechanism -> one (or a few) auxiliary ODE(s) on the boundary node
     (or in the absorbing layer); fold into the FETD / Newmark-beta system
     (the boundary admittance / material update each step).
  4. The passivity guarantee => stable time stepping; the sparsity => few
     auxiliary variables.

This converts a per-frequency-designed absorbing BC (or lossy layer) into a
single broadband, causal, local-in-time operator -- the time-domain
representation needed for a transient (Newmark-beta) open-boundary PoC.
"""

_TOPICS = {
    "overview": URN_OVERVIEW,
    "method": URN_METHOD,
    "api": URN_API,
    "timedomain": URN_TIMEDOMAIN,
    "cq": URN_CQ,
    "application": URN_APPLICATION,
}


def get_urn_documentation(topic: str = "all") -> str:
    """Return URN knowledge text.  topic in
    {all, overview, method, api, timedomain, cq, application}."""
    t = (topic or "all").strip().lower()
    if t == "all":
        return "\n".join(_TOPICS[k] for k in
                         ["overview", "method", "api", "timedomain", "cq", "application"])
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic '{topic}'. Options: all, "
            + ", ".join(_TOPICS) + ".")


def run_urn_fit(freqs, Z, n_debye=3, n_cole_cole=2, n_warburg=1,
                n_cole_davidson=0, sparsity_weight=0.01,
                n_epochs=2000, n_restarts=3, spice=True):
    """Fit a complex frequency response Z(freqs) with a URN and return a dict
    {nrmse, n_active, mechanisms, spice_netlist}.

    freqs : array of frequencies in Hz.   Z : complex array (same length).
    Lower n_epochs / n_restarts for a faster (rougher) fit."""
    import numpy as np
    import torch
    from radia.urn import (
        URNConfig, train_urn, generate_spice_netlist)

    freqs = np.asarray(freqs, dtype=float).ravel()
    Z = np.asarray(Z, dtype=complex).ravel()
    if freqs.shape != Z.shape:
        raise ValueError("freqs and Z must have the same length")

    cfg = URNConfig(n_debye=int(n_debye), n_cole_cole=int(n_cole_cole),
                    n_warburg=int(n_warburg), n_cole_davidson=int(n_cole_davidson),
                    sparsity_weight=float(sparsity_weight),
                    n_epochs=int(n_epochs), n_restarts=int(n_restarts))
    model = train_urn(freqs, Z, cfg, verbose=False)

    omega = torch.tensor(2 * np.pi * freqs, dtype=torch.float64)
    Zp = model(omega).detach().numpy()
    denom = np.sqrt(np.mean(np.abs(Z) ** 2))
    nrmse = float(np.sqrt(np.mean(np.abs(Zp - Z) ** 2)) / (denom + 1e-30))

    mech = model.get_active_components()
    n_active = sum(len(v) for v in mech.values())
    out = {"nrmse": nrmse, "n_active": int(n_active), "mechanisms": mech}
    if spice:
        try:
            out["spice_netlist"] = generate_spice_netlist(model, "Z")
        except Exception as e:  # noqa: BLE001
            out["spice_netlist"] = f"(SPICE synthesis failed: {e})"
    return out


def _format_report(res: dict) -> str:
    lines = [f"URN fit: NRMSE = {res['nrmse']:.4e}   active mechanisms = {res['n_active']}",
             "discovered relaxation mechanisms (tau in seconds):"]
    for mtype, comps in res["mechanisms"].items():
        if not comps:
            continue
        for c in comps:
            parts = [f"tau={c.get('tau', float('nan')):.3e}"]
            for key in ("alpha", "beta"):
                if key in c and c[key] is not None:
                    parts.append(f"{key}={c[key]:.3f}")
            if "weight_magnitude" in c:
                parts.append(f"|w|={c['weight_magnitude']:.3e}")
            if c.get("branch"):
                parts.append(str(c["branch"]))
            lines.append(f"  {mtype}: " + "  ".join(parts))
    if "spice_netlist" in res:
        lines.append("\nSPICE netlist (== auxiliary-ODE ladder for FETD):")
        lines.append(res["spice_netlist"])
    return "\n".join(lines)


def urn_fit_from_csv(data_csv, freq_col=0, real_col=1, imag_col=2,
                     delimiter=",", skip_rows=0, n_debye=3, n_cole_cole=2,
                     n_warburg=1, n_cole_davidson=0, sparsity_weight=0.01,
                     n_epochs=2000, n_restarts=3, spice_out="") -> str:
    """Read a (freq, Re, Im) CSV, fit a URN, return a text report.  If spice_out
    is a path, the SPICE netlist is also written there."""
    import numpy as np
    if not os.path.isfile(data_csv):
        return f"CSV not found: {data_csv}"
    data = np.loadtxt(data_csv, delimiter=delimiter, skiprows=int(skip_rows))
    freqs = data[:, int(freq_col)]
    Z = data[:, int(real_col)] + 1j * data[:, int(imag_col)]
    res = run_urn_fit(freqs, Z, n_debye=n_debye, n_cole_cole=n_cole_cole,
                      n_warburg=n_warburg, n_cole_davidson=n_cole_davidson,
                      sparsity_weight=sparsity_weight, n_epochs=n_epochs,
                      n_restarts=n_restarts, spice=True)
    if spice_out and "spice_netlist" in res:
        with open(spice_out, "w") as f:
            f.write(res["spice_netlist"])
        res = dict(res)
        res["spice_netlist"] = (res["spice_netlist"]
                                + f"\n* (also written to {spice_out})")
    return _format_report(res)
