# Human Netgen .vol viewer

Run `python -B -m radia_mcp.mesh_viewer "FILE.vol"` in an unlocked interactive
Windows session. `--check` prints the actual Python, Netgen and native
module identities without opening the GUI. Missing files, other extensions and
mismatched runtimes fail explicitly. This first delivery accepts `.vol` only.

The launcher uses the installed Netgen GUI and passes the entire filename as
one `tk.call("Ng_LoadMesh", path)` argument. Japanese text, spaces, braces,
brackets, dollar signs and semicolons remain literal path data. It does not
patch Netgen, use Tcl eval for filenames, mesh, solve, write files, or modify
Windows associations. Only its owned window is destroyed on close.

The `open_mesh_viewer(filepath)` MCP tool uses this same module in a detached
process and returns `pid`, literal `argv`, desktop and runtime identity.
Preflight is asynchronous, and the server never waits for the GUI event loop.
`status="launched"` does not by itself certify a rendered image.

An operator can select an existing compatible dependency interpreter with
`RADIA_MESH_VIEWER_PYTHON`. The child's Python path selects the same reviewed
package as the server; it does not install, copy or repoint a shared runtime.
The default server interpreter must itself meet the Netgen version requirement;
there is no implicit alternate interpreter or solver fallback.

A real GUI proof uses a checked `.vol`, records the original file hash and label
contract, confirms visible mesh counts and a Japanese/space/symbol filename,
then verifies the input hash is unchanged. Keep machine/user evidence private.
Gmsh/MSH support and all-user associations are separate follow-up work.

## Declared dependency contract

The viewer feature uses the existing `radia-mcp[gmsh]` extra; no separate viewer
extra or ad-hoc Gmsh installation is introduced. Install through the declared
requirement and official development index:

```text
pip install --pre --extra-index-url https://gmsh.info/python-packages-dev "radia-mcp[gmsh]"
```

Netgen comes from the Radia/NGSolve dependency and must match the sole exact
`netgen-mesher` pin in the installed Radia package metadata. The `.vol`
route does not import Gmsh. MSH/Gmsh support and any all-user installation still
require their separate reviewed follow-up and explicit installation decision.
