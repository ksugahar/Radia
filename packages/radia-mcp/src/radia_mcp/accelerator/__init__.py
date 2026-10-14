"""radia_mcp.accelerator — Accelerator magnet design with Radia.

Specialized knowledge for the design and validation of:
- Conventional accelerator magnets (dipole, quadrupole, sextupole)
- Superconducting magnets (cyclotron, synchrotron)
- End-pole chamfer analytical design (Delferriere et al)
- Multipole content analysis (integrated, harmonic)
- Rotating coil measurement / field reconstruction

Pairs with:
  - radia_mcp.electromagnet (general accelerator magnet workflow)
  - radia_mcp.differential_forms (theoretical underpinnings)
  - Radia core (the actual field solver, BNL ESRF heritage)

Distilled from:
  - public-safe curated corpus (18 PDFs)
  - Delferriere-de Menezes-Duperrier (CEA Saclay) on end pole design
  - Pradhan et al, Kolkata SC Cyclotron (RADIA validation)
  - SOLEIL, ESRF, RIBF design documents
"""
