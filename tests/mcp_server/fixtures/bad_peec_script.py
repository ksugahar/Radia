"""Intentionally bad PEEC/BEM script for lint testing.

Triggers PEEC and BEM-specific rules.
Used by --selftest when examples/ directory is not available.

Expected findings: 3+
"""
import numpy as np
from scipy.special import jv

# Rule: bessel-jv-not-iv (CRITICAL)
# SIBC circular conductor
ka = 1 + 1j
Z_s = jv(0, ka) / jv(1, ka)

# Rule: peec-low-nseg (MODERATE)
# PEEC coil model
n_seg = 16

# Rule: peec-p-over-jw (HIGH)
# Loop-Star PEEC formulation
omega = 2 * np.pi * 1e6
Z_SS = self.P / (1j * omega)

# Rule: efie-v-minus-sign (HIGH)
# ShieldBEMSIBC solver
Z_LL = Zs * M_LL - 1j * omega * mu_0 * V_LL

# Rule: classical-efie-breakdown (MODERATE)
from ngbem import HelmholtzSL
P = 1 / kappa ** 2 * V_SL
