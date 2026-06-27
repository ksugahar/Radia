"""Mixed Galerkin CLN-SIBC: bulk Foster + polyhedral surface envelope -> alpha(s)."""

from .alpha import (
    bulk_foster_via_eigen,
    bulk_foster_matrix_via_eigen,
    wedge_function,
    c1_polyhedral,
    K_SIBC_total,
    Y_mixed,
    alpha_from_Y,
    surface_moment_matrix,
    K_SIBC_matrix,
    Y_matrix_mixed,
    alpha_matrix_from_Y,
    measure_total_area,
    measure_total_area_and_edges,
)
from .cad_edges import (
    cad_topology_edges,
    cad_topology_faces,
    cad_topology_total_area,
    cad_topology_c1,
    edge_moment_matrix,
)
from .vector_bulk import bulk_foster_vector_via_eigen
from .rom_fit import (
    FosterROM,
    passive_foster_fit,
    diagonal_tensor_state_space,
)
from .references import (
    K_SIBC_cylinder,
    K_SIBC_sphere,
    Y_DC_cylinder,
    Y_DC_sphere,
    Y_exact_cylinder,
    Y_exact_sphere,
)

__all__ = [
    "bulk_foster_via_eigen",
    "bulk_foster_matrix_via_eigen",
    "wedge_function",
    "c1_polyhedral",
    "K_SIBC_total",
    "Y_mixed",
    "alpha_from_Y",
    "surface_moment_matrix",
    "K_SIBC_matrix",
    "Y_matrix_mixed",
    "alpha_matrix_from_Y",
    "measure_total_area",
    "measure_total_area_and_edges",
    "cad_topology_edges",
    "cad_topology_faces",
    "cad_topology_total_area",
    "cad_topology_c1",
    "edge_moment_matrix",
    "bulk_foster_vector_via_eigen",
    "FosterROM",
    "passive_foster_fit",
    "diagonal_tensor_state_space",
    "K_SIBC_cylinder",
    "K_SIBC_sphere",
    "Y_DC_cylinder",
    "Y_DC_sphere",
    "Y_exact_cylinder",
    "Y_exact_sphere",
]
