# Third-party notices

See the component texts in LICENSE.

## External runtimes (not bundled)

| Dependency | License | Use |
| --- | --- | --- |
| NGSolve / Netgen | LGPL-2.1 | Dynamically linked, separately installed pip dependencies for finite elements and meshing. |
| NVIDIA CUDA runtime pip packages | NVIDIA EULA | Optional dependency runtimes (for example through PyTorch); never bundled. MATLAB mexcuda uses MATLAB's supplied toolkit. |
| Intel MKL | Intel Simplified Software License | Separate pip dependency for dense BLAS/LAPACK; never bundled. |
| Intel compiler runtime | Intel runtime license supplied with the dependency | Separately installed pip runtime for native dependencies; never bundled. |

Other Python dependencies are installed separately by pip, not copied into these
packages. Their installed distributions provide their own license texts.
Optional extras select the workflows below; development extras provide checks
and documentation tools. This inventory does not relicense dependencies.

## Declared Python dependencies

| Distribution | Dependency | Selected use / extra |
| --- | --- | --- |
| radia | h5py | runtime |
| radia | numpy | runtime |
| radia | scipy | runtime |
| radia | threadpoolctl | runtime |
| radia | mkl | runtime |
| radia | intel-cmplr-lib-rt | runtime |
| radia | ngsolve | runtime |
| radia | netgen-mesher | runtime |
| radia | optuna | optuna |
| radia | cubit-mesh-export | cubit |
| radia | scikit-image | topopt-cad |
| radia | trimesh | topopt-cad |
| radia | netCDF4 | topopt-cad |
| radia | fast-simplification | topopt-cad |
| radia | build123d | build123d |
| radia | xtrack | beam |
| radia | matplotlib | viz |
| radia | anywidget | docs |
| radia | jupyter | docs |
| radia | matplotlib | docs |
| radia | nbconvert | docs |
| radia | nbformat | docs |
| radia | torch | urn |
| radia | pytest | test |
| radia | pytest-cov | test |
| radia | setuptools | test |
| radia | pytest | dev |
| radia | pytest-cov | dev |
| radia | setuptools | dev |
| radia | matplotlib | dev |
| cubit-mesh-export | mcp | runtime |
| cubit-mesh-export | netgen-mesher | runtime |
| cubit-mesh-export | ngsolve | runtime |
| cubit-mesh-export | numpy | runtime |
| cubit-mesh-export | gmsh | mesh-quality |
| cubit-mesh-export | netCDF4 | sculpt |
| cubit-mesh-export | youtube-transcript-api | youtube |
| radia-mcp | mcp | runtime |
| radia-mcp | tomlkit | maintenance |
| radia-mcp | build123d | build123d |
| radia-mcp | cubit-mesh-export | cubit |
| radia-mcp | cadquery | cadquery |
| radia-mcp | gmsh | gmsh |
| radia-mcp | numpy | gmsh |
| radia-mcp | scipy | gmsh |
| radia-mcp | matplotlib | gmsh |
| radia-mcp | numpy | numeric |
| radia-mcp | scipy | numeric |
| radia-mcp | cubit-mesh-export | shape-regeneration |
| radia-mcp | radia | shape-regeneration |
| radia-mcp | scipy | shape-regeneration |
| radia-mcp | gmsh | shape-regeneration |
| radia-mcp | youtube-transcript-api | youtube |
| radia-mcp | radia | radia |
| radia-mcp | scipy | radia |
| radia-mcp | chromadb | rag |
| radia-mcp | sentence-transformers | rag |
| radia-mcp | pymupdf | rag |
| radia-mcp | easyocr | ocr |
| radia-mcp | google-cloud-vision | ocr |
| radia-mcp | radia | urn |
| radia-mcp | eqnedit64 | equation |
| radia-mcp | markdown | md2html |
| radia-mcp | pymdown-extensions | md2html |
| radia-mcp | pymupdf | document |
| radia-mcp | python-pptx | document |
| radia-mcp | lxml | document |
| radia-mcp | Pillow | document |
| radia-mcp | opencv-python-headless | document |
| radia-mcp | matplotlib | document |
| radia-mcp | numpy | document |
| radia-mcp | requests | document |
| radia-mcp | edge-tts | document |
| radia-mcp | pywin32 | document |
| radia-mcp | radia-mcp | full |
| radia-optuna | optuna | upstream |
| radia-optuna | scipy | upstream |
| radia-optuna | torch | upstream |
| radia-optuna | scikit-learn | upstream |
| radia-optuna | cmaes | upstream |
| radia-optuna | build | test |
| radia-optuna | pytest | test |
| radia-optuna | wheel | test |
| radia-simulink | numpy | exchange |
| radia-simulink | h5py | exchange |
| radia-simulink | radia | coupling |
| radia-simulink | schemdraw | ltspice |
| radia-simulink | pytest | test |
| radia-simulink | build | test |
| radia-simulink | numpy | test |
| radia-simulink | h5py | test |
| radia-simulink | schemdraw | test |
| radia-simulink | mcp | test |
