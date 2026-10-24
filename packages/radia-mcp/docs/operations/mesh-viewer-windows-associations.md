# Windows mesh associations: review and installation plan

This is a design, not an installer. No registry or runtime changes are made by
the mesh-viewer launcher. Installation requires the reviewed wheel, an approved
dependency-resolution report and the user's explicit installation decision.

## Runtime choices

The existing system Python can host the reviewed `radia-mcp[gmsh]` requirement
only after confirming which development/release environments inherit its site
packages. Installing there may change their imports. A system installation that
would affect a shared release runtime is excluded from this plan; use the
protected-prefix option instead.

A protected non-venv prefix offers separation without a new interpreter or
virtual environment. Its resolver and bootstrap procedure must be reviewed
before application; do not assume that adding `--prefix` makes an installed
dependency visible to a later pip command.

Resolve and review the complete declared dependency closure before applying;
do not independently install Gmsh or use `--no-deps`. Record the exact wheels,
versions, hashes and provenance. The canonical runtime recipe is the
[two-step Gmsh installation](gmsh-development-runtime.md):

```text
pip install --pre --upgrade -i https://gmsh.info/python-packages-dev gmsh
pip install "radia-mcp[gmsh]"
```

These commands describe an approved environment, not authorization to change
system Python. A prefix-specific deployment plan must bind the qualified Gmsh
wheel and the reviewed `radia-mcp[gmsh]` dependency closure without changing the
shared release runtime. This document does not yet prescribe that installation.

The existing system Radia/Netgen dependency can be used read-only if its installed
metadata pin, native identity and accessibility for other users are verified.
VOL post conversion additionally requires a Radia distribution containing the
production MSH 4.2 exporter; an older distribution fails explicitly. Source-only
proofs do not certify the installed distribution. Any additional distribution
installation must appear in the reviewed installation plan.

A reviewed bootstrap beneath the protected prefix prepends only its own
site-packages, selects its protected Gmsh provenance record, and invokes the same
`radia_mcp.mesh_viewer` module with literal arguments. It contains no second
viewer implementation. Invoke it using the existing Python with `-I -B`;
do not change global PATH, `.pth` files, editable sources or shared runtimes.
Verify both the Python module and Gmsh DLL resolve inside the intended prefix,
and check directory ownership, writable ancestors and reparse targets before
machine-wide registration.

## Exact association scope

Register only `.vol`, `.msh` and two dedicated ProgIDs:
`Radia.MeshViewer.vol` and `Radia.MeshViewer.msh`. Each command uses the protected
Python/bootstrap paths with a quoted literal `"%1"`. `.vol` opens Netgen; `.msh`
opens Gmsh after MSH 4.2 preflight. Double-click never performs conversion or
chooses an implicit output path.

Only the approved defaults, dedicated `OpenWithProgids` values and dedicated
ProgID display-name/`shell\open\command` entries may change. Preserve unrelated
values and subkeys even beneath the two extensions. Do not register DDE/COM,
delete other ProgIDs, or load offline user hives.

Preserve Windows `UserChoice`. No hash bypass, deletion or permission changes
are allowed. HKLM registration makes a handler available to all users; it does
not prove that it overrides each user's existing choice. Users with an override
select Open with / Always themselves in their own interactive account.

## Evidence before machine-wide application

1. Qualify the exact reviewed runtime and wheel hashes. Prove current-user real
   Netgen and Gmsh opens with checked inputs and Japanese/space/symbol paths.
2. Capture registry value names, types and data in both registry views for the
   association surface and approved targets, plus the current user's FileExts.
   Record inaccessible scope explicitly; do not claim unobserved users passed.
3. Produce a dry-run value-level plan. Apply only after review and the explicit
   installation decision. Retain before/after raw differences, including any
   background changes; never suppress them to claim an unchanged registry.
4. Require every association mutation to fall within the approved value-level
   allowlist. Report unrelated changes for independent adjudication. Verify the
   command and actual runtime again from the intended user's desktop.

Machine-wide availability, current-user execution and other users' effective
defaults are separate claims with separate evidence.
