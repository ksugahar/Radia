jou to build123d Translation Benchmark Report

Benchmark run: 2026-04-19
Evaluator: Claude (Opus 4.7, 1M context) restricted to lab mcp-server knowledge
  (src/radia/mcp_server/build123d/build123d_knowledge.py,
   src/radia/mcp_server/cubit/cubit_scripting_knowledge.py).
Runtime: build123d 0.10.0 on Python 3.12, Windows.

== Summary verdict ==

The mcp-server knowledge base was SUFFICIENT to correctly classify and
translate all three fixtures. The cubit_rosetta (build123d side) and
build123d_crossref (cubit side) topics are essentially identical
verb-by-verb tables; between those two plus the lab_policy topic on
each side, an LLM can decide TRANSLATE vs SKIP with no ambiguity for
these fixtures. The two TRANSLATE fixtures produce valid, labeled
build123d Parts with correct bounding boxes and geometrically plausible
volumes on the first attempt; the SKIP fixture is refused with a
specific line-by-line citation of the "no build123d equivalent"
table from the knowledge base. Overall readiness: READY for the
simple primitives+boolean+sweep/revolve/fillet regime that the rosetta
covers, with small knowledge gaps (listed below) that only appeared
in verification rather than in the translation itself.

== Per-fixture scoring ==

--- 01_pure_cad.jou ---

- Verdict: TRANSLATE (expected: TRANSLATE) -- correct.
- Python validity: PASS. python 01_pure_cad.py runs to completion
  with no exceptions. part.is_valid == True.
- Geometric correctness:
    - volume = 87433.6294 mm^3. Analytic check: plate = 100*50*20 =
      100000; two cylinder holes of r=10, h=20 -> 2 * pi * 100 * 20 =
      12566.37; 100000 - 12566.37 = 87433.63. EXACT MATCH.
    - bbox.min = (-50, -25, -10), bbox.max = (50, 25, 10),
      bbox.size = (100, 50, 20). Matches create brick x 100 y 50 z 20
      centered at origin.
    - part.label = "bracket" preserved.
- Idiomaticity: Algebra mode, one expression per .jou line, header
  comment cites each source line. Matches EXAMPLES_INTRO Ex 2.

--- 02_sweep_revolve.jou ---

- Verdict: TRANSLATE (expected: TRANSLATE) -- correct.
- Python validity: PASS. part.is_valid == True.
- Geometric correctness:
    - volume = 6175.3147 mm^3. Pre-fillet analytic: ring with inner
      radius 17.5, outer 22.5, height 10 = pi * 200 * 10 = 6283.19.
      Four r=1 fillets remove ~108 mm^3, physically consistent.
    - bbox.size = (45, 45, 10): outer diameter 45, height 10. EXACT.
    - Label "ring_bobbin" preserved.
- Idiomaticity: Plane.XZ * Rectangle(5, 10) for the revolve profile
  aligns with EXAMPLES_INTRO Ex 23 and the rosetta revolve entry.
  Fillet via fillet(ring.edges(), radius=1). Comments explain the
  zplane vs Plane.XZ remap.

--- 03_imprint_hex.jou ---

- Verdict: SKIP (expected: SKIP) -- correct.
- Python validity: PASS as stub. Prints refusal, exits 0.
  build() raises NotImplementedError with explanatory message.
- Geometric correctness: N/A by design.
- Idiomaticity: Refusal cites specific .jou lines (L17-L18
  imprint/merge, L21 sweep scheme, L22-L23 size/mesh, L26-L29
  block/name) and maps each to the "no build123d equivalent"
  bullet in cubit_rosetta / build123d_crossref.

== Knowledge gaps observed ==

1. part.is_valid attribute vs method. TOPOLOGY shows attribute form
   (correct for v0.10), but the knowledge base does not warn that
   older training data often writes is_valid(). A one-line note in
   CAE_GUIDELINES would harden this footgun.

2. Plane.XZ * Rectangle(...) placement idiom not directly in the
   knowledge base. The rosetta revolve entry does not show how to
   place the 2D sketch in a plane that contains the revolution axis.
   EXAMPLES_INTRO Ex 23 uses Polyline + make_face + split. A worked
   example mapping Cubit zplane rectangle + move + revolve z to
   build123d Plane.XZ * Rectangle + revolve Axis.Z would close this.

3. reset has no mapping entry. A no-op row in the rosetta would
   remove a common student question.

4. Compound-vs-Part for labeled assemblies. The rosetta maps Cubit
   group to Compound but does not note that Compound loses the
   fused-solid property; users often want v1 + v2 + v3 for meshing.

5. Fillet for "all edges of a volume". The rosetta fillet row does
   not explicitly show fillet(part.edges(), radius=R) for the common
   "all edges" case.

None of these gaps prevented correct translation. Gaps 1 and 2 caused
small friction during verification only. The knowledge base did not
mislead; it was just silent on these specific points.

== Overall readiness ==

Readiness: GOOD for the CAD-subset translation task. The lab policy
tables are unambiguous about TRANSLATE vs SKIP and give the correct
refusal for fixture 03 without prompting the LLM to invent
imprint/merge equivalents. The verb mapping tables are sufficient for
fixtures 01 and 02 and produced correct geometry on the first attempt.

Recommended investment, in priority order:

1. Worked "revolve profile on Plane.XZ" example in cubit_rosetta or
   EXAMPLES_INTRO Ex 23. Highest-frequency conceptual remap. (Gap #2)
2. One-line note in CAE_GUIDELINES or TOPOLOGY that is_valid is an
   attribute in current build123d. (Gap #1)
3. "All edges of a volume" row in the fillet/chamfer section of the
   rosetta. (Gap #5)
4. A jou_translation_recipe appendix topic chaining the decision
   tree: scan for imprint/merge/block/mesh -> STOP and keep .jou
   if present; else walk each CAD verb via the rosetta; verify
   is_valid, volume, bounding_box; set part.label.

No blocking deficiencies. Current mcp-server state is already good
enough to hand off routine .jou translation tasks to an LLM with
the knowledge base loaded.
