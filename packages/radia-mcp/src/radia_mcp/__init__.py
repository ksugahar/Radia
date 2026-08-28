"""radia-mcp: MCP servers for Radia CAE ecosystem.

Layout:
  radia_mcp.cubit         — standalone Cubit MCP server (Plan A: GUI + QTimer + file drop)
  radia_mcp.build123d     — standalone build123d MCP server (STEP → Cubit workflow)
  radia_mcp.gmsh          — standalone gmsh MSH v4.1 inspect/validate +
                            post-display launch artifacts + policy lint
  radia_mcp.radia_ngsolve — Radia + NGSolve general FEM/BEM/Kelvin/PEEC/MSH knowledge
  radia_mcp.force         — Shared electromagnetic-force method selection,
                            solver-independent sample integration, and validation;
                            common to Motor and MagLev
  radia_mcp.ih            — Induction Heating workflow (workpiece SIBC, ESIM,
                            Karl iteration, screening physics) — promoted from
                            legacy private source tree on 2026-04-24
  radia_mcp.peec          — LAB PEEC workflow (Loop-Star, FastHenry, PEECBuilder,
                            PEECCircuitSolver, Bessel/Dowell/ESIM, PRIMA MOR,
                            SPICE extraction) — promoted from
                            legacy private source tree on 2026-04-24
  radia_mcp.electromagnet — Accelerator electromagnet (CoilBuilder, Hantila polarization,
                            B-input Play/Energy hysteresis, IMA sign selection,
                            multipole harmonics) — promoted from
                            legacy private source tree on 2026-04-24
  radia_mcp.differential_forms — Differential geometry for computational EM:
                            tangent spaces, k-forms, wedge product, exterior
                            derivative, Hodge star, Whitney elements, de Rham
                            complex, Maxwell's house, tree-cotree gauging,
                            FEEC (Arnold-Falk-Winther 2006). Distilled from
                            Bossavit 1998 + FEEC 2006 + Kameari 2011 + 微分形式
                            lecture notes.
  radia_mcp.mathematica   — Wolfram Mathematica subprocess bridge: evaluate /
                            status / simplify / to_tex / check_identity /
                            vector_calc / unit_convert / solve / integrate /
                            differentiate via wolframscript. Promoted from
                            legacy private source tree
                            on 2026-05-20. Pairs with differential_forms for
                            symbolic verification of Maxwell identities.
  radia_mcp.topology_optimization — Shape and topology optimization for nonlinear
                            magnetostatics. Gangl-Sturm 2015 Lagrangian method
                            for nonlinear PDE shape derivatives; Sokolowski-
                            Zochowski topological derivative; IPM motor cogging-
                            torque minimization case study. Distilled from
                            arXiv:1501.04752 + Gangl PhD thesis Part I/II.
  radia_mcp.mor           — Model Order Reduction, centered on **Cauer Ladder
                            Network (CLN)** — a LAB SPECIALTY where Sugahara is
                            co-author on the canonical papers (Kameari-Ebrahimi-
                            Sugahara-Shindo-Matsuo 2018 et seq). General MOR +
                            CLN basic + multi-expansion + nonlinear + applications.
  radia_mcp.pinn          — Physics-Informed Neural Networks (Raissi 2019) and
                            Gaussian Processes (Raissi 2017, Pförtner 2023) for
                            Maxwell's equations. Inverse problems, multi-fidelity,
                            high-dim PDE. Pairs with mathematica for symbolic
                            verification of kernel transformations.
  radia_mcp.accelerator   — Accelerator magnet design with Radia: analytical
                            end-pole chamfer (Delferriere SOLEIL), Kolkata
                            superconducting cyclotron Radia+TOSCA validation
                            (Pradhan 2007), rotating-coil multipole measurement
                            and 3D field reconstruction.
  radia_mcp.presentation  — Research-talk slide lint + PPTX tools (promoted
                            2026-06-02 from LAB-private mcp-server-document;
                            2026-07-17 served by mcp-server-paper-writing —
                            the standalone server was retired).
  radia_mcp.grant_writing — Grant proposal lint, recommendation-letter
                            templates, and KDDI social-implementation checks
                            (promoted 2026-06-27 from the document stack).
  radia_mcp.poster        — Conference poster generation + lint (promoted
                            2026-06-02 from mcp-server-document).
  radia_mcp.pdf           — PDF manipulation: merge/split/metadata/watermark/
                            compress (promoted 2026-06-02 from mcp-server-document).
  radia_mcp.doc_convert   — Document format conversion: PPTX<->PDF, PDF->JPG,
                            slide extraction (promoted 2026-06-02 from
                            mcp-server-document).
  radia_mcp.bibliography  — BibTeX / citation tooling: DOI/arXiv->bibtex,
                            CrossRef search, parse/dedupe/lint (promoted
                            2026-06-02 from mcp-server-document).
  radia_mcp.document_meta — Cross-cutting helpers: deadline countdown, version
                            diff, LaTeX templates, all-domain lint orchestration
                            (promoted + REDESIGNED 2026-06-02 from
                            the document stack -- lint_all is now a registry
                            over radia-mcp's own document lints).
  radia_mcp.research_project — Project dashboard / scan: consistency, deadline
                            gantt, cross-document-type health (promoted
                            2026-06-02 from the document stack; grant health
                            now uses radia_mcp.grant_writing).

The radia_ngsolve / ih / peec / electromagnet servers reference radia
from inside example snippets in their knowledge modules, but importing
them does not require radia to be installed (knowledge is plain text +
Python). The `radia` extra (`pip install radia-mcp[radia]`) brings in
the runtime dependency for tools that actually exec radia code.
"""

__version__ = "1.4.47"
