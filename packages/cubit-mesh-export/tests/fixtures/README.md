# Named curve sideset fixture

`bbnd_brick.py` independently authors a cube of side length 1 m, centred at the
origin, and meshes it with tetrahedra of nominal size 0.5 m. It labels the
vertical edge at x = y = -0.5 m using the curve sideset `test_edge`.
`bbnd_brick.vol` is the native export at order 2; no post-load naming or Python
mesh rewriting is used. `bbnd_brick.contract.json` owns its versioned labels.

Regenerate on a licensed Cubit 2026.8 machine with the built candidate DLL in
an isolated plugin directory:

```powershell
$env:CUBIT_PLUGIN_DIR = '<candidate-plugin-directory>'
$env:BBND_OUTPUT = '<output-directory>/bbnd_brick.vol'
cubitx.exe -nographics -batch -nojournal -noinitfile bbnd_brick.py
```

The runtime ACIS library directory must be on PATH (the package's
`runtime_environment` helper supplies it). `BBND_ORDER=1` checks affine export;
`BBND_CASE=conflict` or `unrepresented` exercises expected export failures.

The Cubit-free test loads the committed fixture directly in NGSolve, checks
the selected edge's position and total length, and runs check-vol against the
contract. A missing fixture fails rather than skipping. Fixtures and numerical
values here come solely from this elementary cube.
