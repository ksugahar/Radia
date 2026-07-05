"""radia_mcp.magnetic_materials: comprehensive magnetic material knowledge.

Topics covered:
- Hysteresis models (13+ models with decision tree)
- Iron loss models (Steinmetz, Bertotti, MSE, iGSE, Carstensen)
- Silicon steel database (JIS grades + commercial constants)
- Permanent magnet datasheets (NdFeB, SmCo, Ferrite, AlNiCo)
- Demagnetization factor (Osborn 1945 + numerical)
- Radia implementation status (MatPlay, MatEnergy, MatSatIsoTab)

Distilled from public-safe curated corpus (151 files, 2.7 GB).
"""

from .hysteresis_models_knowledge import get_hysteresis_models_knowledge
from .iron_loss_knowledge import get_iron_loss_knowledge
from .silicon_steel_knowledge import get_silicon_steel_knowledge
from .permanent_magnet_knowledge import get_permanent_magnet_knowledge
from .demagnetization_knowledge import get_demagnetization_knowledge
from .radia_status_knowledge import get_radia_status_knowledge

__all__ = [
    "get_hysteresis_models_knowledge",
    "get_iron_loss_knowledge",
    "get_silicon_steel_knowledge",
    "get_permanent_magnet_knowledge",
    "get_demagnetization_knowledge",
    "get_radia_status_knowledge",
]
