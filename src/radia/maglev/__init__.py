"""radia-maglev: Mixed Galerkin CLN-SIBC framework for magnetic
levitation and eddy-current brake analysis.

Public API:
    radia.maglev.mixed_galerkin -- alpha(s) from any .vol mesh
    radia.maglev.ecb            -- eddy-current brake force computation
    radia.maglev.simulink       -- MATLAB / Simulink LTI export
"""
__version__ = "0.1.0"

from . import mixed_galerkin
from . import ecb
from . import simulink

__all__ = ["mixed_galerkin", "ecb", "simulink", "__version__"]
