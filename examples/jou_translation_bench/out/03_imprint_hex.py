"""
Translation of 03_imprint_hex.jou to build123d.

Verdict: SKIP (do NOT translate).

Rationale per lab policy (src/radia/mcp_server/build123d/build123d_knowledge.py
topic `lab_policy` and `cubit_rosetta`, mirrored in cubit_scripting_knowledge
topic `build123d_crossref`):

    "If more than ~20% of a .jou file is in the 'no equivalent' table
     above, don't translate it - keep the .jou, run it through Cubit,
     export hex to Radia/ELF. build123d is the wrong tool for that
     design."

This file has the following commands that have NO build123d equivalent
and are explicitly listed as "keep in Cubit":

  .jou L17  imprint all          -- build123d has no shared-topology
                                    multi-body imprint operation. OCCT
                                    Compound does not enforce conformal
                                    faces across bodies.
  .jou L18  merge all             -- same; build123d has no merge-faces
                                    primitive for conformal meshing.
  .jou L21  volume all scheme sweep  -- hex-sweep scheme is Cubit-only;
                                        build123d does not mesh, and
                                        Netgen (the downstream mesher)
                                        produces tets, not hex.
  .jou L22  volume all size 2     -- meshing directive; in the build123d
                                    pipeline this lives in Netgen as
                                    OCCGeometry(...).GenerateMesh(maxh=2).
  .jou L23  mesh volume all       -- meshing is out of build123d's scope.
  .jou L26-L29  block 1 ... name "iron_core" / block 2 ... "winding"
                                  -- block tagging for Exodus export has
                                    no direct build123d equivalent. The
                                    documented substitute is
                                    part.label = "..." + Compound +
                                    run_pipeline_multi (Gmsh physical
                                    groups), but that is a downstream
                                    pipeline choice, not a translation.

The critical point is that this design *requires* hex (shared-topology,
conformal hex mesh via sweep scheme on an imprinted/merged multi-body
assembly), which is exactly the case the lab policy names as Cubit's
permanent role: "Radia/ELF requires hex - only Cubit can supply."

Keeping this file as `.jou` and running it through Cubit is the correct
action. build123d is the wrong tool.

No runnable Python geometry is produced on purpose. The pure-CAD prefix
(two bricks + subtract-keep + inner brick) *could* be expressed in
build123d as three separate solids, but that would silently discard the
imprint/merge/mesh/block semantics that are the actual point of the
script, so we refuse the translation entirely rather than produce a
misleading partial port.
"""


def build():  # noqa: D401 - intentionally not implemented
    raise NotImplementedError(
        "03_imprint_hex.jou is not translatable to build123d. "
        "It uses imprint/merge + hex-sweep + block tagging, which have "
        "no build123d equivalent per the lab policy. Keep the .jou and "
        "run it through Cubit for the hex path."
    )


if __name__ == "__main__":
    # Make the refusal explicit when someone runs the script.
    print("SKIP: 03_imprint_hex.jou - kept as Cubit .jou (hex path).")
    print("See module docstring for the full rationale.")
