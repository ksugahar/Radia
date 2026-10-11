"""radia-maglev: Foster + SIBC framework for magnetic levitation and
eddy-current brake analysis.

Public API:
    radia.maglev.mixed_galerkin -- alpha(s) from any .vol mesh
    radia.maglev.ecb            -- eddy-current brake force computation
    radia.maglev.state_space    -- CPU state-space realization
"""
__version__ = "0.1.0"

from . import mixed_galerkin
from . import ecb
from . import state_space
from .position_force import PositionForceCurve
from .position_foster import MovingHCurlFosterFamily

__all__ = [
    "mixed_galerkin",
    "ecb",
    "state_space",
    "MovingHCurlFosterFamily",
    "PositionForceCurve",
    "__version__",
]
