r"""Universal Relaxation Network (URN) -- causal/passive rational fitting of a

SHOWCASE NOTEBOOK: docs/universal_relaxation_network/urn_showcase.ipynb -- 4 verified paper figures (URN-vs-VF ~22.8%, NASA/TDK fits, attention ablation).
frequency response, with direct time-domain (relaxation-network / SPICE / ADE)
synthesis.

URN is a KAN-inspired network that decomposes a measured or computed frequency
response Z(omega) (impedance, dispersive permittivity/permeability, or an
open-boundary DtN symbol G_n(omega)) into a sparse sum of PHYSICAL relaxation
mechanisms -- Debye, Cole-Cole, Cole-Davidson, Havriliak-Negami, CPE, Warburg,
Gerischer, RLC, skin-effect -- in both a series and a parallel (admittance)
branch, with KAN-style adaptive relaxation-time (tau) refinement and a
frequency-dependent attention gate.  Because every basis is a passive,
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
Synthesis," IEEE Access, 2026.  Canonical implementation + paper + benchmarks:
examples/universal_relaxation_network/.

API (this module):
  get_urn_documentation(topic) -- knowledge text (this module)
  run_urn_fit(freqs, Z, ...)   -- fit and return mechanisms + NRMSE + SPICE
  urn_fit_from_csv(path, ...)  -- read a (freq, Re, Im) CSV, fit, return a report
"""

import os
import sys

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
with a frequency-dependent ATTENTION gate a_k(omega) (small MLP, softmax) that
lets each basis dominate only in its own frequency band (sharpens the
decomposition; ~79-83% accuracy gain on real data in the ablation).

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
"""

URN_API = r"""
# URN API (canonical impl: examples/universal_relaxation_network/)

from universal_relaxation_network import (
    UniversalRelaxationNetwork, URNConfig, train_urn, generate_spice_netlist)

config = URNConfig(n_debye=3, n_cole_cole=2, n_warburg=1,
                   sparsity_weight=0.01, n_epochs=2000, n_restarts=3)
model  = train_urn(freqs_Hz, Z_complex, config, verbose=False)
mech   = model.get_active_components()       # {type: [{tau, alpha, weight_magnitude,...}]}
spice  = generate_spice_netlist(model, "Z")  # SPICE netlist string (RC/RL ladders)

URNConfig fields: n_debye/n_cole_cole/n_cole_davidson/n_havriliak_negami/n_cpe/
  n_warburg/n_gerischer/n_rlc/n_skin_effect (series), *_parallel (admittance),
  sparsity_weight, lr, n_epochs(=6000), n_restarts(=10), omega_ref/Z_ref(auto),
  use_attention(=True).  For a responsive tool call, lower n_epochs/n_restarts.

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
    "application": URN_APPLICATION,
}


def get_urn_documentation(topic: str = "all") -> str:
    """Return URN knowledge text.  topic in
    {all, overview, method, api, timedomain, application}."""
    t = (topic or "all").strip().lower()
    if t == "all":
        return "\n".join(_TOPICS[k] for k in
                         ["overview", "method", "api", "timedomain", "application"])
    if t in _TOPICS:
        return _TOPICS[t]
    return (f"Unknown topic '{topic}'. Options: all, "
            + ", ".join(_TOPICS) + ".")


def _resolve_urn_core() -> str:
    """Locate the canonical URN implementation (examples/universal_relaxation_network)
    by walking up from this file, and put it on sys.path.  Returns the path."""
    here = os.path.dirname(os.path.abspath(__file__))
    cur = here
    for _ in range(12):
        cand = os.path.join(cur, "examples", "universal_relaxation_network")
        if os.path.isfile(os.path.join(cand, "universal_relaxation_network.py")):
            if cand not in sys.path:
                sys.path.insert(0, cand)
            return cand
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    raise ImportError(
        "URN core not found: expected examples/universal_relaxation_network/"
        "universal_relaxation_network.py above " + here)


def run_urn_fit(freqs, Z, n_debye=3, n_cole_cole=2, n_warburg=1,
                n_cole_davidson=0, sparsity_weight=0.01,
                n_epochs=2000, n_restarts=3, spice=True):
    """Fit a complex frequency response Z(freqs) with a URN and return a dict
    {nrmse, n_active, mechanisms, spice_netlist}.

    freqs : array of frequencies in Hz.   Z : complex array (same length).
    Lower n_epochs / n_restarts for a faster (rougher) fit."""
    _resolve_urn_core()
    import numpy as np
    import torch
    from universal_relaxation_network import (
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
