"""radia-levitation: Mixed Galerkin CLN-SIBC framework for magnetic
levitation and eddy-current brake analysis.

Public API:
    radia.levitation.mixed_galerkin -- alpha(s) from any .vol mesh
    radia.levitation.ecb            -- eddy-current brake force computation
    radia.levitation.simulink       -- MATLAB / Simulink LTI export
"""
__version__ = "0.1.0"

from . import mixed_galerkin
from . import ecb
from . import simulink

__all__ = ["mixed_galerkin", "ecb", "simulink", "__version__"]
