# Time-domain implementation

Warburg block expansion → finite-dimensional LTI realisation.

| Script | What it does |
|---|---|
| `diffusive_quadrature.py` | Construct 50-point diffusive Foster expansion $\sqrt{s}/(s+d) \approx \sum c_k/(s+\xi_k)$; saves `warburg_diffusive_fit.npz` |
| `time_domain_verify_warburg.py` | Compare 50-pt Foster step response vs analytic Dawson $(2/\sqrt{\pi d})\,\mathrm{dawsn}(\sqrt{dt})$ ── max err $1.8\times 10^{-3}$ |
| `time_domain_verify_cylinder.py` | Bessel cylinder step response: exact (analytic) vs hierarchical Cauer LTI |
| `extract_hierarchical_cauer.py` | Production extractor: frequency-sweep Y(s) → hierarchical Cauer parameters (bulk + Warburg) for SPICE / state-space |
