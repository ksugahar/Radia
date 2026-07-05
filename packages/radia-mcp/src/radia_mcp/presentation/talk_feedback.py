"""
talk_feedback.py -- learned field-note catalog of REAL conference talks
(CEFC / Compumag / Intermag / IGTE / IEEJ ...) given by the lab.

Purpose
=======
The single highest-signal "training data" for presenting well at a given
venue is what that venue's audience actually ASKED, and what the session
chair / audience actually reacted to.  This is the presentation analog of
``meta/bug_patterns.py``: make the recurring Q&A and lessons explicit and
queryable, so the next talk is prepared with real history instead of
guesswork.

Two-phase entry lifecycle
=========================
- ``status="anticipated"`` -- written BEFORE the talk: the questions we
  predict (with planned answers) + the rules we are presenting against.
- ``status="recorded"``    -- updated AFTER the talk: the questions that
  were *actually* asked, the chair/audience feedback, what worked, what to
  improve.  THIS is the real venue signal; capture it within a day while
  memory is fresh.

How to query (MCP tools, registered via presentation/tools.py)
==============================================================
    presentation_qa_from_history(topic="validation")   # ★ Q&A prep
    presentation_qa_from_history(topic="FEM-comparison")
    presentation_talk_feedback_lookup(venue="CEFC")
    presentation_talk_feedback_lookup(topic="H-matrix", status="recorded")
    presentation_talk_feedback_stats()

How to add an entry after a talk
================================
Append a dict to FIELD_NOTES (or flip an "anticipated" entry to
"recorded" and fill in the real Q&A).  Keep it SHORT -- this ships in the
radia-mcp wheel; long narrative belongs in memory/ or docs/.  Fields:

  id        : stable slug (venue-year-topic)
  venue     : "CEFC 2026" / "Compumag 2025" ...
  date      : YYYY-MM-DD
  location  : city, country
  talk      : title
  presenter : name
  topics    : keyword tags (method / validation / FEM-comparison / scaling ...)
  status    : "anticipated" | "recorded"
  questions : list of {q, a, lesson}
  chair_audience_feedback : list[str]  (post-talk)
  what_worked / what_to_improve : list[str]  (post-talk)
  rules     : list[str]  -- distilled, feed these into the next talk's prep
"""
from __future__ import annotations

from collections import Counter

FIELD_NOTES: list[dict] = [
    # =====================================================
    # CEFC 2026 -- H-matrix accelerated HDiv-VIM (seed, anticipated)
    # =====================================================
    {
        "id": "cefc2026-hdiv-hacapk",
        "venue": "CEFC 2026",
        "date": "2026-06-09",
        "location": "Thessaloniki, Greece",
        "talk": "H-matrix Accelerated HDiv-VIM for Open-Boundary "
                "Magnetostatics",
        "presenter": "K. Sugahara",
        "topics": ["HDiv-VIM", "H-matrix", "HACApK", "ACA",
                   "accelerator-magnet", "validation", "FEM-comparison",
                   "scaling", "nonlinear", "open-boundary"],
        "status": "anticipated",
        "questions": [
            {
                "q": "How do accuracy and computational cost compare with FEM?",
                "a": "Field agrees with the FEM reference within 0.3 mT. HDiv-VIM "
                     "needs no air mesh and gives an exact open boundary; the "
                     "price is a dense matrix, which the H-matrix compresses "
                     "(204 GB -> 8.6 GB, ~24x).",
                "lesson": "This audience benchmarks EVERYTHING against FEM -- "
                          "put the FEM comparison early and quantitatively "
                          "(outline slide + key-result slide).",
            },
            {
                "q": "Does the ACA / H-matrix approximation degrade accuracy "
                     "or the iteration count?",
                "a": "No measurable degradation: 0.3 mT vs FEM at eps=1e-4; "
                     "~100 Picard iterations, comparable to a dense solve.",
                "lesson": "Keep an eps-vs-(accuracy, compression, iters) "
                          "backup slide.",
            },
            {
                "q": "How does memory / time scale with N? Is it sub-cubic?",
                "a": "166k DOF in 8.6 GB (4.2%); the trend suggests ~1e6 DOF "
                     "is reachable on a desktop. Honest caveat: the H-LU "
                     "factor cost dominates beyond ~15^3, so the demonstrated, "
                     "field-correct regime is mu_r<=1e4.",
                "lesson": "Do NOT overclaim scaling. State the demonstrated "
                          "point (desktop-PC scale) and have the honest "
                          "factor-cost caveat ready if pressed.",
            },
            {
                "q": "Did you validate against a TEAM benchmark or measurement, "
                     "not only FEM?",
                "a": "Here the reference is a converged FEM solution "
                     "(Bz within 0.3 mT). TEAM-problem / measurement "
                     "validation is planned.",
                "lesson": "ANTICIPATE 'TEAM / measurement?' at every CEM talk. "
                          "Have an honest one-liner; ideally add one TEAM "
                          "problem before the next talk.",
            },
            {
                "q": "Why HDiv-VIM + HACApK rather than BEM with an FMM "
                     "(e.g. ngsolve.bem)?",
                "a": "Compact, near-field-heavy magnet geometries are ~87% "
                     "near pairs, where HACApK ACA+ wins; FMM-BEM suits "
                     "smooth-surface, far-field-dominated problems. Different "
                     "geometry classes.",
                "lesson": "Be ready to defend the method choice vs FMM/BEM -- "
                          "the CEM audience knows both intimately.",
            },
            {
                "q": "Is this 3D, and does it handle nonlinear iron?",
                "a": "Yes -- 3D hexahedral surface-charge elements; nonlinear "
                     "via a Fixed-Point -> Newton hybrid (the 166k case "
                     "converged in 100 Picard iterations).",
                "lesson": "State 3D + nonlinear explicitly and early; both are "
                          "assumed must-haves at CEFC/Compumag.",
            },
            {
                "q": "Is the code available? HACApK in C or Fortran?",
                "a": "The HACApK C version is integrated into Radia "
                     "(pip install radia); a Fortran integration is planned.",
                "lesson": "Open-source / availability questions are common -- "
                          "have the one-line answer.",
            },
        ],
        "chair_audience_feedback": [],   # fill post-talk
        "what_worked": [],               # fill post-talk
        "what_to_improve": [],           # fill post-talk
        "rules": [
            "CEM (CEFC/Compumag) audiences benchmark every new method against "
            "FEM -- lead with a quantitative FEM comparison; keep "
            "cost/accuracy-vs-FEM on a backup slide.",
            "Always anticipate 'did you validate against TEAM / measurement?' "
            "-- prepare an honest one-liner; ideally add one TEAM problem.",
            "Have backup slides ready for: eps-vs-accuracy, memory/time "
            "scaling, method-choice vs FMM/BEM, the nonlinear scheme.",
        ],
    },
    # =====================================================
    # Historical orals (reconstructed from the lab's own decks;
    # actual asked-Q&A not transcribed -- status="historical")
    # =====================================================
    {
        "id": "cefc2020-retired-moment-path",
        "venue": "CEFC 2020",
        "date": "2020-11-03",
        "location": "Pisa, Italy (online)",
        "talk": "Retired Surface-Charge Inspired Magnetic Solver Study",
        "presenter": "H. Yano",
        "topics": ["HDiv-VIM", "accelerator-magnet", "convergence",
                   "FEM-comparison", "validation", "CST", "OPERA-3D",
                   "mesh-dependence"],
        "status": "historical",
        "questions": [
            {
                "q": "How do accuracy / convergence compare with the "
                     "commercial standards (OPERA-3D)?",
                "a": "The whole talk is built on it: a convergence study on a "
                     "T-shaped magnet (N=1..10) vs a stored commercial reference, "
                     "plus C-shaped and quadrupole mesh-dependence -- the proposed "
                     "method converges far faster than plain retired moment path and matches that "
                     "reference.",
                "lesson": "Position vs the SPECIFIC tools the audience uses "
                          "(OPERA-3D / Radia) and SHOW a convergence study "
                          "against them, over several example geometries.",
            },
            {
                "q": "Computational cost vs OPERA-3D?",
                "a": "OPERA-3D 72 min (8-layer IABC, 1/8 model) vs 2 min.",
                "lesson": "Lead the cost slide with one dramatic number.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Three verification examples (T / C / quadrupole) + a 10-point "
            "convergence study = very strong validation.",
            "Speaker notes written with SSML <break> tags to hit time.",
        ],
        "what_to_improve": [],
        "rules": [
            "Open by positioning vs the commercial de-facto standard "
            "(OPERA-3D TOSCA / CST), then show the meshing gap concretely.",
            "Validate with a convergence study vs a commercial reference "
            "across multiple example geometries.",
        ],
    },
    {
        "id": "compumag2023-cln-team28",
        "venue": "Compumag 2023",
        "date": "2023-05-01",
        "location": "Kyoto, Japan",
        "talk": "Cauer Ladder Network Representation with Constant Basis "
                "Functions for Eddy Current Problems Involving Conductor "
                "Movement",
        "presenter": "N. Tanimoto",
        "topics": ["CLN", "MOR", "eddy-current", "TEAM28", "motion",
                   "FEM-comparison", "BEM", "speedup", "Simulink"],
        "status": "historical",
        "questions": [
            {
                "q": "Why CLN instead of FEM or BEM for a moving-conductor "
                     "eddy-current problem?",
                "a": "Opened the talk with it: FEM needs an air mesh "
                     "(complicated for motion), BEM is time-consuming; the "
                     "reduced CLN keeps R/L constant w.r.t. levitation height, "
                     "so Simulink coupling is easy.",
                "lesson": "Open a CEM talk by positioning vs FEM (air mesh) and "
                          "BEM (cost) explicitly -- the audience expects it.",
            },
            {
                "q": "How much faster, and validated how?",
                "a": "Transient analysis in 5 s vs ~40 min for commercial EM "
                     "software; validated on TEAM Workshop Problem 28 (maglev).",
                "lesson": "TEAM 28 is the lab's go-to maglev benchmark; pair it "
                          "with one dramatic speedup number.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Tight 4-slide main body (intro / difference / result) + backup "
            "formulation slides -- good for a short Compumag oral.",
        ],
        "what_to_improve": [],
        "rules": [
            "For a short oral, keep ~4 main slides and push formulation "
            "detail to backup.",
        ],
    },
    {
        "id": "cefc2024-fp-cln",
        "venue": "CEFC 2024",
        "date": "2024-06-02",
        "location": "Jeju, Korea",
        "talk": "Fixed-Point Cauer Ladder Network Method for Eddy-current "
                "Problems",
        "presenter": "K. Sugahara",
        "topics": ["CLN", "MOR", "eddy-current", "nonlinear", "fixed-point",
                   "TEAM28", "analytic-solution", "validation"],
        "status": "historical",
        "questions": [
            {
                "q": "How does CLN differ from other model-order-reduction "
                     "methods (PVL, POD)?",
                "a": "Dedicated 'Model Order Reduction in EM' slide: unlike "
                     "PVL/POD, CLN derives the network directly as an "
                     "approximate solution of Maxwell's equations.",
                "lesson": "Have a slide positioning your method vs the other "
                          "named methods in the field (PVL/POD/...).",
            },
            {
                "q": "How is the nonlinear / fixed-point scheme validated "
                     "without FEM?",
                "a": "Examples A/B/C (1D sheet, 1D rod, 2D rod) have analytic "
                     "CLN solutions -- 'you do not have to use FEM to "
                     "construct basis functions'; plus TEAM 28 for the "
                     "extended CLN.",
                "lesson": "Analytic-solution validation is very strong at "
                          "CEFC; build from the simplest 1D analytic case up.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Credited the originator (CLN by Kameari 2016) up front.",
            "Examples built in increasing difficulty (1D -> 1D -> 2D).",
        ],
        "what_to_improve": [],
        "rules": [
            "Credit the method's originator early; show the limitation of the "
            "conventional version on its own slide before the extension.",
        ],
    },
    # =====================================================
    # CEFC 2016 -- Dirichlet IABC open-boundary BC
    # (deck new01_cefc2016; title-only extract)
    # =====================================================
    {
        "id": "cefc2016-dirichlet-iabc",
        "venue": "CEFC 2016",
        "date": "2016-11-13",
        "location": "Miami, USA",
        "talk": "A Dirichlet-type Infinite/Asymptotic Boundary Condition "
                "(IABC) for Open-Boundary Magnetostatic Field Computation",
        "presenter": "K. Sugahara",
        "topics": ["IABC", "open-boundary", "Dirichlet-BC",
                   "analytic-solution", "validation", "PML",
                   "Kelvin-transformation", "demagnetization-factor",
                   "magnetostatic"],
        "status": "historical",
        "questions": [
            {
                "q": "Why a Dirichlet-type IABC rather than PML, ballooning, "
                     "infinite elements, or a far truncated air box?",
                "a": "The IABC imposes the asymptotic far-field decay directly "
                     "as a Dirichlet condition at a finite boundary, so there "
                     "is no PML layer to tune, no infinite-element book-"
                     "keeping, and no large air box -- yet the open boundary "
                     "is reproduced exactly on the canonical test shapes.",
                "lesson": "Name every competing open-boundary technique the "
                          "audience already uses (PML / Kelvin / ballooning / "
                          "truncation) and state in one line what your BC "
                          "removes from each.",
            },
            {
                "q": "How is the boundary condition validated -- only on a "
                     "sphere?",
                "a": "On a LADDER of geometries that each have a known closed-"
                     "form internal/external field and demagnetization factor: "
                     "sphere, cylinder, ellipsoid, elliptic cylinder, plus a "
                     "periodic case -- ordered by increasing anisotropy so the "
                     "method is stress-tested, not shown on one easy shape.",
                "lesson": "Validate a new open-boundary/BC method on an "
                          "ordered ladder of analytic-solution geometries "
                          "(isotropic -> uniaxial -> biaxial -> periodic); "
                          "each rung is an exact, falsifiable check.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Method-first framing: slide 1 states the one-line identity "
            "('Dirichlet IABC') before any geometry.",
            "One geometry per two-slide block (setup then result) with a "
            "divider between -- every claim is independently auditable in "
            "the oral and maps back to its exact analytic reference.",
        ],
        "what_to_improve": [
            "Title-only extract carries NO headline number -- a per-geometry "
            "error figure (e.g. max field error vs analytic) on each result "
            "slide would let the audience size the accuracy at a glance.",
        ],
        "rules": [
            "Lead with the method NAME on slide 1; let the canonical shapes "
            "do the persuading after.",
            "Order the analytic-solution validation ladder by increasing "
            "difficulty -- the ordering itself tells the reviewer you "
            "stress-tested anisotropy.",
            "For a BC paper, position explicitly against PML / Kelvin / "
            "ballooning / truncation -- that is the comparison set this "
            "audience carries in their heads.",
        ],
    },
    # =====================================================
    # CEFC 2022 -- Multiport CLN of a skewed cage IM
    # (deck new03_cefc2022)
    # =====================================================
    {
        "id": "cefc2022-multiport-cln-skew-im",
        "venue": "CEFC 2022",
        "date": "2022-10-24",
        "location": "Denver, USA",
        "talk": "Multiport Cauer Ladder Network Model-Order Reduction of a "
                "Cage Induction Motor with Skewed Rotor Slots",
        "presenter": "K. Sugahara",
        "topics": ["CLN", "MOR", "multiport", "induction-motor", "skew",
                   "multi-slice", "space-harmonics", "FEM-comparison",
                   "speedup", "control-plant-model"],
        "status": "historical",
        "questions": [
            {
                "q": "Why a multiport CLN rather than the original single-port "
                     "(scalar) CLN, or a general POD reduced-order plant model?",
                "a": "A skewed rotor carries several space harmonics across "
                     "the slices, so a scalar CLN cannot represent it; the "
                     "multiport CLN carries multiple SHs and -- unlike a "
                     "fitted POD model -- is a passive lumped Cauer circuit "
                     "derived as an approximate Maxwell solution.",
                "lesson": "Scope the contribution to the ONE genuinely new "
                          "piece (multiport carrying space harmonics + inter-"
                          "slice bar-current continuity) and credit the "
                          "originating CLN (Cauer / Kameari lineage).",
            },
            {
                "q": "How is the reduced model validated, and how big is the "
                     "problem (DOF)?",
                "a": "Three CLN strategies (Method 1 parallel per-slice, "
                     "Method 2 direct 3-D stepwise, Method 3 continuity-aware) "
                     "are cross-validated against each other and the 2-D "
                     "multi-slice FEM: phase-current waveforms overlaid slice-"
                     "by-slice, and the Cauer circuit reproduces the FEM B and "
                     "eddy-current J from the CLN basis vectors. Problem size "
                     "is stated as Ns slices x M space harmonics = MNs DOF.",
                "lesson": "Earn trust before claiming speed -- show the "
                          "reduced model reproduces the FEM field FIRST, then "
                          "the 'computational time' slide reads 'fast AND "
                          "faithful'. State the DOF (MNs) so reviewers can "
                          "size it.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Controlled comparison of 3 numbered methods on a fixed colour "
            "legend (Method 1 green / 2 blue / 3 red) reused on every results "
            "slide -- accuracy AND cost differences readable at a glance.",
            "Cost slide as the climax (dedicated 'Comparison of computational "
            "time') after the field-fidelity slides.",
            "Duplicate 'Table of contents' slide splitting formulation half "
            "from results half -- a navigation reset before results.",
        ],
        "what_to_improve": [
            "The dramatic speedup (lab CLN precedent: hours-per-frequency FEM "
            "-> minutes/real-time) is implied by the dedicated cost slide but "
            "the exact figure did not survive the extract -- put one bold "
            "speedup number on the cost slide.",
        ],
        "rules": [
            "Make the talk a controlled comparison, not a single-method "
            "pitch: numbered method variants on a fixed colour legend let "
            "the structure do the persuading.",
            "Anchor a niche MOR method to a real engineering need on slide 2 "
            "(model-based motor control needs an accurate plant model).",
            "Mechanism-before-results: derive the equivalent Cauer circuit "
            "and show B/J from basis vectors before any validation, so the "
            "reduced model reads as physical, not curve-fit.",
        ],
    },
    # =====================================================
    # Compumag 2025 -- Warm-Starting CMA-ES topology opt (IPM rotor)
    # (deck new08_compumag2025)
    # =====================================================
    {
        "id": "compumag2025-warmstart-cmaes-ipm",
        "venue": "Compumag 2025",
        "date": "2025-06-08",
        "location": "Thessaloniki, Greece",
        "talk": "Data-Driven Topology Optimization of IPM-Motor Rotors via "
                "Warm-Starting CMA-ES for Constrained Design",
        "presenter": "Yama- (山, exact name unconfirmed)",
        "topics": ["topology-optimization", "CMA-ES", "warm-start",
                   "IPM-motor", "torque-ripple", "iron-loss",
                   "constrained-optimization", "penalty-method",
                   "GA-comparison"],
        "status": "historical",
        "questions": [
            {
                "q": "How does Warm-Starting CMA-ES compare with a GA or with "
                     "plain (cold-start) CMA-ES?",
                "a": "A dedicated search-performance slide races Warm-Starting "
                     "CMA-ES against GA, and the warm-started variant is "
                     "compared against the un-warm-started CMA-ES baseline; "
                     "the warm start seeds parameters from the top-gamma "
                     "individuals so the search begins from good flux-barrier "
                     "shapes.",
                "lesson": "For an optimization talk, race the proposed "
                          "variant against the NAMED baselines (GA, plain "
                          "CMA-ES) the audience expects -- method-vs-baseline "
                          "framing is the persuasion.",
            },
            {
                "q": "Does it only work on your one example?",
                "a": "Two FE design cases are shown: Case A (torque-ripple "
                     "minimization) and Case B (iron-loss minimization under "
                     "an average-torque constraint, V-shape rotor) framed as "
                     "'additional verification' of generality, each with an "
                     "initial-solution-vs-optimized-shape comparison.",
                "lesson": "Validate by GENERALITY, not just one case: Case A "
                          "then a second 'additional verification' case "
                          "answers the reviewer's 'only your example?' before "
                          "it is asked.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Recurring identical agenda slide reused as a section divider + "
            "progress indicator before each of the 3 acts -- keeps a "
            "methodology-heavy talk legible.",
            "Side-by-side 'initial solution (warm start) vs optimized shape' "
            "comparison slides -- the persuasive artifact for an opt talk.",
        ],
        "what_to_improve": [
            "No headline number survived: the gains are shown as Case-A/"
            "Case-B shape + objective comparisons but with no single dramatic "
            "figure (iterations-to-converge or objective gap) -- add one for "
            "a CEFC/Compumag oral.",
        ],
        "rules": [
            "A single agenda slide reused as a progress divider keeps a "
            "3-act methodology talk legible -- the audience always knows "
            "which act they are in.",
            "Show the init-vs-optimized shape pair raced against a named "
            "baseline (GA / plain CMA-ES), and quantify the win with ONE "
            "headline number.",
            "Frame a second case as 'additional verification' to pre-empt "
            "the generality question.",
        ],
    },
    # =====================================================
    # Compumag 2025 -- B-H loop fitting for energy-based hysteresis
    # (deck new09_compumag2025; data-conditioning sub-component)
    # =====================================================
    {
        "id": "compumag2025-bh-loop-fitting",
        "venue": "Compumag 2025",
        "date": "2025-06-08",
        "location": "Thessaloniki, Greece",
        "talk": "Hysteresis-Curve Fitting from Measured B-H Loops for an "
                "Energy-Based Stop/Play Model (data-conditioning pipeline)",
        "presenter": "Sugahara lab (Kindai); data credited to Toyo University",
        "topics": ["hysteresis", "B-H-loop", "curve-fitting", "measurement",
                   "energy-based-model", "Stop-model", "Play-model",
                   "synthetic-data", "validation"],
        "status": "historical",
        "questions": [
            {
                "q": "Where does the B-H data come from, and does it cover the "
                     "full operating range?",
                "a": "Measured minor/major loops from Toyo University, "
                     "credited on every data slide, swept one-slide-per-"
                     "amplitude-band across six 0.25 T windows from 0.05 to "
                     "1.50 T -- so the full envelope is shown to be data-"
                     "backed, not cherry-picked, before any fit.",
                "lesson": "Spend a real block PROVING you cover the "
                          "measurement space before modelling it; credit the "
                          "borrowed measurement data on the slide where it "
                          "appears, not just in acknowledgements.",
            },
            {
                "q": "How is the fit validated, and why fit instead of a "
                     "Preisach / Jiles-Atherton analytic model?",
                "a": "The reference IS the measurement: the rational-function "
                     "Bmax-Hmax envelope + a normalized 17th-order polynomial "
                     "must reproduce the measured loop family, and the "
                     "regenerated synthetic loops (constant dHmax/dBmax "
                     "increments) are checked for self-consistency. It is the "
                     "data-driven alternative that supplies the energy-based "
                     "Stop/Play model's input curves.",
                "lesson": "When measurement IS the ground truth, judge the "
                          "method by curve-fit fidelity to the measured loops "
                          "and self-consistency of the regenerated data; keep "
                          "the dense fit machinery in a backup deck.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Clean measure -> fit -> generate arc: raw loops first, then "
            "envelope/polynomial extraction, then synthetic-loop generation.",
            "One-slide-per-amplitude-band sweep (6 bands) is a credibility "
            "move -- the audience sees the whole 0.05-1.50 T envelope.",
        ],
        "what_to_improve": [
            "Working/backup deck with no dramatic cost number; for the oral, "
            "present only the measure->fit->generate storyline and bank the "
            "key fit-fidelity number, not every coefficient.",
        ],
        "rules": [
            "Prove measurement-space coverage before modelling (one slide per "
            "amplitude band).",
            "Credit the originator of borrowed measurement data on its own "
            "slide.",
            "Keep dense data-conditioning math (rational envelope, 17th-order "
            "polynomial, synthetic generation) in a separate backup deck.",
        ],
    },
    # =====================================================
    # CEFC 2026 (sibling new10) -- IEM-FEM weak coupling for maglev
    # (Yano; remote/pre-recorded delivery)
    # =====================================================
    {
        "id": "cefc2026-iem-fem-maglev-coupling",
        "venue": "CEFC 2026",
        "date": "2026-06-09",
        "location": "Thessaloniki, Greece",
        "talk": "Bi-directional IEM-FEM Weak-Coupling Framework (Radia HDiv-VIM + "
                "NGSolve FEM) for Electromagnetic Analysis of Maglev Systems",
        "presenter": "Takaaki Yano",
        "topics": ["IEM-FEM", "weak-coupling", "maglev", "eddy-current",
                   "Radia", "NGSolve", "A-phi", "T-Omega",
                   "cross-validation", "Kelvin-transformation",
                   "reduced-potential", "MCP-server", "LLM",
                   "open-boundary", "FEM-comparison"],
        "status": "historical",
        "questions": [
            {
                "q": "You have no TEAM benchmark or commercial-code "
                     "comparison -- how do you know the result is correct?",
                "a": "The SAME rotating+translating PM-cube eddy-current "
                     "problem (1 mm cube, Br = 0.2 T, X-translation -6 to "
                     "4 mm with rotation) is solved two independent ways -- "
                     "the A-phi and the T-Omega potential formulations -- and "
                     "their time histories agree to 0.21% mean relative error "
                     "in Joule loss and 0.11% in Lorentz force.",
                "lesson": "When you lack a TEAM/commercial reference, solve "
                          "the SAME problem two independent ways (A-phi vs "
                          "T-Omega) and report the small mutual error -- "
                          "formulation cross-validation is a credible self-"
                          "contained correctness story.",
            },
            {
                "q": "Isn't this just Radia? What did you actually add, and "
                     "why couple it to FEM at all?",
                "a": "Radia is credited explicitly as ESRF (Grenoble) and we "
                     "state we forked + extended it. Pure FEM must mesh the "
                     "air domain and re-mesh on magnet motion with open-"
                     "boundary truncation error; the integral route gives the exact open-"
                     "boundary external field analytically, FEM solves only "
                     "the conductor reaction field (reduced potential), and "
                     "the two exchange fields via NGSolve CoefficientFunction "
                     "/ GridFunction -- so magnet motion needs only a field "
                     "update, no re-mesh.",
                "lesson": "Always name and credit the upstream tool (Radia = "
                          "ESRF) and state precisely what you added ('forked + "
                          "extended'); front-load a 1-2 slide primer on the "
                          "enabling abstraction (CF/GF field exchange) so the "
                          "coupling is legible to a mixed FEM/IEM audience.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Two-formulation cross-validation used as the headline "
            "correctness argument (0.21% Joule / 0.11% force) instead of a "
            "single reference solution.",
            "Tool-primer slide (CoefficientFunction vs GridFunction) before "
            "results so the coupling mechanism is legible.",
            "Cost-vs-DoF + Japanese background kept as titled BACKUP slides "
            "after the Thank-you, for Q&A rather than padding the main arc.",
            "SSML/scripted speaker notes for a remote / non-native-English "
            "delivery (the presenter was absent due to flight disruption).",
        ],
        "what_to_improve": [
            "The speed argument lives only on a 'Computational Cost vs DoF' "
            "backup slide with no explicit minute-level number -- promote one "
            "concrete speed figure into the main arc.",
            "No TEAM-problem / measurement anchor yet -- add one to harden "
            "the maglev validation beyond formulation agreement.",
        ],
        "rules": [
            "Formulation cross-validation (A-phi vs T-Omega, small mutual "
            "error) is a self-contained correctness story when no TEAM/"
            "commercial reference is available.",
            "Front-load a primer on the enabling abstraction before results "
            "for a mixed-audience coupling talk.",
            "Credit the upstream tool by name and say exactly what you "
            "extended; keep cost-vs-DoF + extra background as titled backup "
            "slides for Q&A.",
        ],
    },
    # =====================================================
    # CEFC 2026 (sibling new11) -- self-built field-measurement bench
    # (Doshisha; experimental progress talk in the CEFC2026 bucket)
    # =====================================================
    {
        "id": "cefc2026-field-measurement-bench",
        "venue": "CEFC 2026",
        "date": "2026-06-09",
        "location": "Thessaloniki, Greece",
        "talk": "A Self-Built Magnetic-Field Measurement Rig: Tesla-Meter on "
                "a Custom Frame, Calibrated and Compared Against Simulation",
        "presenter": "Doshisha University member (name not stated)",
        "topics": ["measurement", "experimental", "Tesla-meter",
                   "calibration", "parallelism", "dial-gauge",
                   "measurement-vs-simulation", "coil-design",
                   "field-gradient"],
        "status": "historical",
        "questions": [
            {
                "q": "How accurate is the rig -- can we trust the field "
                     "numbers it produces?",
                "a": "Apparatus is validated before the physics: the named "
                     "commercial sensor (F.W.Bell Model 8030 digital Tesla-"
                     "meter) gets its own instrument + specifications slides, "
                     "and stage parallelism is calibrated with a dial gauge "
                     "shown on a dedicated parallelism-results slide.",
                "lesson": "For an experimental EM talk, validate the "
                          "apparatus before the physics -- name and spec the "
                          "sensor and show the calibration tool + its result; "
                          "a measurement claim is only as credible as the "
                          "stated instrument and its accuracy.",
            },
            {
                "q": "What is the measured field compared against, and what "
                     "is the rig FOR?",
                "a": "Measured field is compared to an in-house simulation "
                     "model, and the validated model is then used downstream "
                     "to derive a coil from the simulated field gradient -- "
                     "closing the loop measurement -> simulation -> design.",
                "lesson": "Close the loop: don't stop at 'measurement matches "
                          "model'; show what the validated model is FOR. A "
                          "talk that ends on a forward use lands harder than "
                          "one that ends on an agreement plot.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Agenda/roadmap slide up front listing the four topics in order "
            "('tell them what you'll tell them').",
            "Methods-transparency for an experimental talk: the dial gauge "
            "gets its own slide BEFORE the parallelism result.",
            "Builds to a measurement-vs-simulation payoff, then forwards the "
            "simulated gradient into a coil design.",
        ],
        "what_to_improve": [
            "No quantitative figure survived (no field magnitude, no "
            "parallelism tolerance, no % error) -- an experimental talk needs "
            "the instrument's accuracy and the measurement-vs-model error "
            "stated on the slide.",
            "Terse Japanese noun-phrase titles + a blank divider read as an "
            "internal progress deck, not a polished CEFC podium deck -- the "
            "filing under the CEFC2026 bucket may be a seminar deck.",
        ],
        "rules": [
            "Validate the apparatus before the physics (calibration tool + "
            "its result visible, not buried).",
            "Name and spec the measurement hardware explicitly (sensor model "
            "+ a specifications slide).",
            "Close the loop measurement -> simulation -> design; end on a "
            "forward use, not an agreement plot.",
        ],
    },
    # =====================================================
    # IGTE 2026 -- Energy-based Stop hysteresis model
    # (deck new12_igte2026; interim status deck. Included as the
    # strongest problem->trap->fix exemplar + an honest missing-
    # validation gap worth cataloging.)
    # =====================================================
    {
        "id": "igte2026-energy-based-stop-hysteresis",
        "venue": "IGTE 2026",
        "date": "2026-09-16",
        "location": "Graz, Austria",
        "talk": "An Energy-Based Stop Model for Magnetic Hysteresis "
                "(Play -> Stop -> Egger Variational Extension)",
        "presenter": "Sugahara-lab hysteresis-modeling line (name not in "
                     "extract)",
        "topics": ["hysteresis", "Stop-model", "Play-model",
                   "energy-based-model", "Egger-2025", "shape-functions",
                   "ridge-penalty", "hysteron-placement", "validation"],
        "status": "historical",
        "questions": [
            {
                "q": "Why move from the Play model to an energy-based Stop "
                     "model -- isn't Play already established?",
                "a": "The first ~5 slides motivate it: the Play model is hard "
                     "to make energy-based, its shape functions favour the "
                     "Stop formulation, and the Egger 2025 variational "
                     "structure gives a convex energy with gamma-controlled "
                     "operating modes -- so Stop is the natural carrier for an "
                     "energy-based extension.",
                "lesson": "Lead with the WHY before the formulation: spend the "
                          "opening on 'why Play -> energy-based is happening' "
                          "so the audience accepts the formulation as "
                          "motivated, not arbitrary. Name the originator "
                          "(Egger 2025) in the slide title.",
            },
            {
                "q": "How is the fitted model validated, and what are its "
                     "limits?",
                "a": "Self-consistency / fit-quality only: measured B-H loops "
                     "are reconstructed from the fitted shape functions "
                     "(showcased at N=25 hysterons), hysteron count + "
                     "placement are swept, and three candidate formulations "
                     "are compared numerically; a non-monotone-shape-function "
                     "trap is cured by a ridge penalty. Explicit limits + "
                     "future-work slides precede the summary.",
                "lesson": "The problem->trap->fix triplet (non-monotone shape-"
                          "function trap immediately followed by the ridge-"
                          "penalty fix) convinces an expert audience far more "
                          "than only the polished final method.",
            },
        ],
        "chair_audience_feedback": [],
        "what_worked": [
            "Motivation-before-math opening (~5 slides on why the extension "
            "is needed) before the Stop formulation.",
            "problem -> trap -> fix narrative on one tight slide group.",
            "Honest explicit limits + future-work slides before the summary "
            "-- the credibility posture IGTE/CEFC reviewers reward.",
        ],
        "what_to_improve": [
            "Validation is self-consistency / loop-reconstruction ONLY -- no "
            "TEAM problem, no commercial-code cross-check, no measurement-vs-"
            "FEM field comparison. For the actual oral, add a hard external "
            "validation and one headline number (the most quantitative "
            "anchor in the deck is just N=25 hysterons).",
        ],
        "rules": [
            "Motivation before math: open on the physical/method WHY so the "
            "formulation reads as motivated.",
            "Use the problem->trap->fix triplet on one slide group -- show "
            "the failure you hit and the cure.",
            "Credit the originator in the slide title (Egger 2025) and keep "
            "explicit limits/future-work slides; an interim status deck still "
            "needs a hard external validation + one headline number before "
            "it is oral-ready.",
        ],
    },
]


# ------------------------------------------------------------------
# helpers / MCP tools (registered by name via presentation/tools.py)
# ------------------------------------------------------------------
def _topic_hit(note: dict, t: str) -> bool:
    if any(t in x.lower() for x in note.get("topics", [])):
        return True
    if t in note.get("talk", "").lower() or t in note.get("id", "").lower():
        return True
    for q in note.get("questions", []):
        if t in (q.get("q", "") + " " + q.get("a", "")).lower():
            return True
    return False


def presentation_talk_feedback_lookup(venue: str = "",
                                      topic: str = "",
                                      status: str = "") -> dict:
    """Query the learned conference-talk field-note catalog.

    USE THIS WHEN PREPARING A TALK for a venue the lab has presented at
    before (CEFC / Compumag / Intermag / IGTE ...).  Each entry records
    the questions asked (or anticipated) and the lessons distilled, so a
    new talk is prepared against real history.

    Args:
        venue: substring filter, e.g. "CEFC", "Compumag 2025". Empty = any.
        topic: substring filter against tags + title + id + Q&A text,
            e.g. "validation", "FEM-comparison", "scaling". Empty = any.
        status: "anticipated" (pre-talk predictions) | "recorded"
            (actually-asked at the talk) | "historical" (reconstructed
            from a past deck; Q&A not transcribed). Empty = any.

    Returns ``{"matched": N, "notes": [...]}`` with full field-note dicts.
    """
    out = FIELD_NOTES
    if venue:
        v = venue.lower()
        out = [n for n in out if v in n.get("venue", "").lower()]
    if status:
        s = status.lower()
        out = [n for n in out if n.get("status", "").lower() == s]
    if topic:
        t = topic.lower()
        out = [n for n in out if _topic_hit(n, t)]
    return {
        "matched": len(out),
        "notes": out,
        "hint": ("Read each note's 'rules' (how to present) and 'questions' "
                 "(what gets asked). status='recorded' = real venue signal; "
                 "'anticipated' = predicted, still to be confirmed."),
    }


def presentation_qa_from_history(topic: str = "") -> dict:
    """★ Q&A REHEARSAL: surface the real (or anticipated) questions asked
    at past talks, to rehearse before a new one.

    The single most useful pre-talk action: run this for your talk's
    topic and make sure you can answer every question in <=30 s.

    Args:
        topic: substring filter against tags / title / Q&A text
            (e.g. "validation", "FEM", "scaling", "nonlinear", "BEM").
            Empty = every recorded/anticipated question.

    Returns ``{"matched": N, "qa": [{venue, status, q, a, lesson}, ...]}``.
    """
    t = (topic or "").lower()
    qa: list[dict] = []
    for n in FIELD_NOTES:
        topical = (not t) or any(t in x.lower() for x in n.get("topics", [])) \
            or t in n.get("talk", "").lower()
        for q in n.get("questions", []):
            qtext = (q.get("q", "") + " " + q.get("a", "")).lower()
            if topical or (t and t in qtext):
                qa.append({
                    "venue": n.get("venue", ""),
                    "status": n.get("status", ""),
                    "q": q.get("q", ""),
                    "a": q.get("a", ""),
                    "lesson": q.get("lesson", ""),
                })
    return {
        "matched": len(qa),
        "qa": qa,
        "hint": ("Rehearse a <=30 s answer to each. 'recorded' = a question "
                 "actually asked at a past talk (highest signal); "
                 "'anticipated' = predicted for an upcoming talk."),
    }


def presentation_talk_feedback_stats() -> dict:
    """Counts of the conference-talk field-note catalog (by venue / status
    / topic). Quick view of how much real venue signal has accumulated."""
    by_venue = Counter(n.get("venue", "?") for n in FIELD_NOTES)
    by_status = Counter(n.get("status", "?") for n in FIELD_NOTES)
    by_topic = Counter(tp for n in FIELD_NOTES for tp in n.get("topics", []))
    n_q = sum(len(n.get("questions", [])) for n in FIELD_NOTES)
    n_recorded_q = sum(len(n.get("questions", [])) for n in FIELD_NOTES
                       if n.get("status") == "recorded")
    return {
        "n_talks": len(FIELD_NOTES),
        "n_questions": n_q,
        "n_questions_recorded": n_recorded_q,
        "by_venue": dict(by_venue),
        "by_status": dict(by_status),
        "top_topics": dict(by_topic.most_common(15)),
        "hint": ("Grow this after every talk: flip the 'anticipated' entry "
                 "to 'recorded' and fill the real Q&A + chair feedback "
                 "within a day."),
    }
