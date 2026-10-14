# Gmsh post-processing validation

This directory owns the checked numerical evidence behind the executed public
notebooks in `docs/gmsh_post/`. The notebooks keep their display outputs but do
not write result JSON or maintain notebook checksum ledgers.

- `em_fieldlines_results.json` records equal-flux and streamline checks.
- `em_post_gallery_results.json` records the common post-processing gallery.
- `em_particle_orbits_results.json` records energy conservation, quadrupole
  focusing, and edge-focusing checks.

Run the lightweight evidence gate with:

```powershell
python -m pytest validation_test/gmsh_post/test_gmsh_post_evidence.py -q
```
