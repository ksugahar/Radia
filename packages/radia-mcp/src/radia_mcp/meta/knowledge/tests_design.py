"""Read-only authoring guide for AI-written electromagnetic analysis code.

Policy belongs to CLAUDE.md. This module explains how to apply it; it neither
approves reference values nor executes cases, mutations or lifecycle changes.
"""

SECTIONS = {
    "overview": ("1. Ownership and three purposes", """
AI writes the implementation; humans own acceptance criteria, reference values
and numerical tolerances. An agent may propose a change with evidence, never
approve its own oracle or silently update a baseline to match current output.

tests prevents future changes from breaking a capability. Its numerical oracle
is a saved output of an identified earlier implementation: a regression lock,
not independent evidence of physical correctness. validation_test checks the
formulation against an analytic solution, independent method or measurement.
validation_records records a reached milestone and the conditions under which
it becomes stale. A saved notebook explains the capability and cites its record.

A surrogate test is a light regression test of an adopted validation case:
execute convergence level 1 (the coarsest/lowest-order level) through the same
production entry and compare against its saved earlier computed output.
Use the full term "surrogate test", not "surrogate" alone: it is not a surrogate
model, response surface, ROM or machine-learning approximation. Its historical
oracle checks change, not agreement with the analytic solution at coarse scale.
The adoption record must show that level 1 belonged to the accepted convergence
sequence. About 30 seconds is a design budget, not a timing assertion.
Assumes must cite: "surrogate test of adoption record <ID> (convergence level 1)".
Replace <ID> only with the actual approved record identifier. Until approval
and level mapping exist, mark the test as a candidate, never a certified claim.
Existing different-scale locks may continue; record the mismatch and align at
the next explicitly reviewed baseline regeneration, without relabeling bytes.

Begin with independent validation, retain the executable case, and create a
small historical lock through the same production entry. Re-certify validation
before intentional baseline changes. Link feature, route, validation, record,
lock and introduction rather than copying results into competing registries.
"""),
    "feature": ("2. Identify the feature and its production entry", """
Write a feature row: public production entry, owner, formulation/backend and
important branches, validation case, light lock and source contracts. Select
one sensitive end-to-end lock per feature rather than enumerating variants.
Keep additional cases only when they have a distinct purpose or failure signal.

Calculation time is scale: mesh, order, elements or frequency samples change the
size of one numerical case. Use a validation scale with an independent oracle
and a light scale that calls the same production entry. Light setup may be
written independently in tests, including a frozen mesh and small runner.
Never duplicate solver logic. Sharing a parameterized case runner is useful,
not mandatory. Record when validation and light scales are intentionally equal.

About 30 seconds is a light-model design budget, not an elapsed-time assertion.
Store measured runtime_s and conditions; let CI report durations. Keep the heavy
validation runnable and do not remove its existing CI coverage before its
replacement validation lane actually runs.
"""),
    "quantity": ("3. Choose a sensitive computed quantity", """
Choose a response that depends on the code being protected: a computed field
sample, integral or coupled response influenced by the declared feature. Explain
the causal path from operator, source and boundary treatment to that quantity.
Bad locks include the applied field itself, copied input arrays and a component
forced to zero by symmetry. They may remain identical while the solver is broken.

Accidental numerical agreement in electromagnetic responses is very unlikely,
so one sensitive quantity can protect a substantial path. A combined test that
requires A AND B AND C supports all three only where they execute and
materially influence its verdict. Matching a number is not a proof of every
component: small contributions can be hidden and rare cancellations can occur.
Demonstrate sensitivity to each claimed Guard with production mutations.
"""),
    "lock": ("4. Create a historical lock", """
1. Retain a passing independent validation receipt with exact source, conditions,
   gates and artifact hashes; do not relabel an old failure as a pass.
2. Run the light case through the declared production entry. Save a machine-
   readable numerical result before coverage export, plotting or other analysis.
3. Repeat at the same source SHA under the relevant admitted thread settings.
   Measure spread, including absolute behavior near zero; propose each quantity's
   atol/rtol with a reproducibility rationale for human review.
4. Store source SHA, conditions, quantity names/shapes/units/tolerances and
   runtime_s metadata. Use tools/regression_lock.py for comparison and proposal-
   only regeneration; do not implement a second comparison framework.
5. Execute coverage comparison and the mutation/restoration checks. Submit those
   results with the baseline proposal to an independent reviewer.

Validation thresholds do not become light-scale assertions. A light model may
have discretization error; its comparison is to its own approved historical
output. Reproducibility bounds justify a proposed lock tolerance, not a change
to the human-owned physics gate. Follow the helper's finite-value, shape and
condition checks. Split complex quantities into real/imaginary parts when needed.
"""),
    "declarations": ("5. Declare Guards and Assumes", """
Every new, moved or newly registered test module states:

Guards: <which future change or failure this test detects>.
Assumes: <entry/backend, scale, conditions, baseline source and applicability>.

For a numerical lock, Assumes also explains how the selected output depends on
the guarded production path. Name only features with sensitivity evidence. A
declaration is a scoped claim, not a list of everything imported by the case.
For a source contract, identify the API, packaging, failure behavior or fixture
compatibility being guarded instead of inventing a numerical oracle.

Keep the feature map, baseline metadata and docstring consistent. If a mutation
survives, remove the unsupported Guard or strengthen the quantity/conditions and
repeat the evidence; never retain a broad claim because the normal case passes.
"""),
    "mutation": ("6. Demonstrate production sensitivity", """
At initial lock creation and every baseline regeneration, select at least one
small artificial corruption per declared Guard in the actual production code.
Examples: perturb a constitutive coefficient by 1e-3, reverse a sign, swap a
dispatch branch, reduce integration order, or omit one interface condition.
Choose a change that exercises that Guard, not an unrelated crash.

Use an owned checkout with no other job reading the mutable source. Save the
normal source identity, exact patch and patch SHA256. Never commit the corruption
or modify tests/mocks to manufacture a failure. Run the actual lock; record the
failed assertion and its actual/expected values. Restore exact original bytes
in finally-style cleanup and run the normal lock again: PASS is required.

Report feature x mutation with patch hash, killed/survived, failure signal,
normal-restoration PASS and evidence identity. Mark interrupted or unrelated-error
runs inconclusive, not killed. A survived corruption is evidence the feature is
not guarded. Add a sensitive quantity/condition or propose a separate contract;
do not loosen tolerance. Mutations are authoring/regeneration evidence, not
mandatory work on every CI run. Confirm a clean tree before normal precheck.
"""),
    "coverage": ("7. Use coverage to find holes", """
Capture validation and light production coverage with matching source/artifacts,
entry, formulation, backend, element/order and branch settings. Use the manual
tools/coverage_difference.py comparator and its declared evidence schema rather
than a second coverage verdict. Keep the original numerical gates and receipts.

Inspect heavy-only lines/branches and the route matrix. For an uncovered path,
choose conditions that exercise it, add a distinct contract for behavior that
does not affect the response, or propose retirement with consumer evidence.
Do not count tests or coverage percentages as correctness. Passing coverage
does not establish numerical sensitivity; mutation checks answer that question.

Missing capture or subprocess sentinels remain UNKNOWN. Named native/JIT entry
hits show entry use only, not internal coverage. Do not invent C++ branch claims
from Python coverage, instrument release kernels as a prerequisite, automatically
demote cases, or turn this manual diagnostic into a new light-CI dependency.
"""),
    "contracts": ("8. Protect behavior that values cannot show", """
Validation is event-driven, not periodically scheduled. Development-only unit
checks can be written freely on a work branch, but are classified one by one
before a PR: stable public promises invisible to end-to-end values belong in
tests with Guards; independent-reference comparisons belong in validation_test;
checks of internal shape alone are disposable scaffolding. Do not create a
permanent unit_tests layer.

Retain small source contracts for abnormal inputs, unsupported states, API and
packaging boundaries, policy, external fixtures and silent-fallback rejection.
A correct number on the normal branch cannot certify an unexecuted error path
or that a backend was explicitly selected. Check the real library's public
behavior and fail loudly; a handcrafted fake is not evidence of library support.

Keep production logic out of tests. Do not mirror the implementation to compute
an expected answer. Historical locks intentionally use past output; independent
physics oracles remain in validation_test. Coverage overlap alone is not a
retirement reason. Retire only complete
redundancy: identical purpose, quantity, model and failure signal. Different
quantities, models, references or routes can be independent detectors worth
keeping when their independence is justified.
"""),
    "review": ("9. Check each change independently", """
Execute at least one independent, failure-capable check for every change and
report exact source, artifact identities and what actually ran. Independence
may come from a different reviewer or an analytic/independent method or
measurement; rereading the same summary or reusing the same derivation is weak.

The 90% to 99% example is intuition for independent errors, not a guarantee or
measured reliability claim. Correlated implementation and checker mistakes can
survive both. Scale scrutiny to impact and follow the stronger requirements for
releases, deletion and human-owned criteria changes. Review oracle independence,
mutation sensitivity and coverage limits, not merely the number of green tests.
"""),
    "regeneration": ("10. Regenerate and retire deliberately", """
Before an intentional baseline change, re-certify the independent validation,
repeat same-SHA reproducibility and production mutation/restoration evidence,
then generate a proposal with reason and numerical diff. The helper must not
overwrite an approved baseline. Human review owns tolerance/reference changes;
a clean source tree and provenance identify the proposal, not its approval.

For retirement, provide exact-SHA consumer-scanner evidence covering imports,
workflows, manifests, docs, knowledge and other packages, then obtain independent
review. No consumers is a necessary check, not deletion permission by itself.
For a complete duplicate, show identical purpose, quantity, model and failure
signal, and name the retained detector.
Preserve any analytic/reference oracle elsewhere. Keep deletion in a small
separate change, not mixed with moves or numerical modifications. History is
provenance; never silently retire an active route or its only gate.
"""),
    "checklist": ("11. Execution preflight and evidence checklist", """
Before execution:
- Resolve dependencies from package metadata and the declared CI profile; do not
  hand-pick an environment until imports happen to pass. Record the interpreter,
  versions and actual import paths in an isolated job environment.
- Check attributes/functions against the real installed library, using its
  exposed signature and a small valid call when appropriate. Handwritten fake
  APIs do not prove the real binding supports an operation.
- Dry-run child-interpreter startup in the actual worker environment, including
  its import paths, identity and exit status, before scheduling numerical work.
- Use native artifacts only with known generating source SHA and verified hash;
  document any same-source reuse with the production/build-input comparison.
- Construct child environments with one function returning a fresh mapping.
  Record explicit threads/backend/paths and test that environment; do not use
  global text replacement to edit commands or environment settings.
- Persist numerical results and identities before optional post-processing
  (including coverage export). Keep partial/failed evidence if later steps fail.

Submission checklist:
feature/entry and owner; Guards/Assumes; validation receipt with unchanged gates;
light scale and sensitive output; baseline/proposed tolerances with same-SHA
repeatability; runtime_s as metadata only; coverage report and UNKNOWN limits;
feature x production mutation/patch SHA256/failure/restored PASS table; exact
source/artifact hashes; proposal reason/diff when regenerating; independent
review and required precheck. Preserve originals. A complete checklist is not
itself a numerical PASS or permission to publish private comparisons.
"""),
    "adoption": ("12. Agree adoption criteria before implementation", """
Before implementation, propose the feature's production entry and a criterion:
independent reference type/source (analytic solution, independent route,
measurement or benchmark), sensitive quantity and units, exact comparison
operator/threshold, scale/conditions and failure definition. A previous output
from the same solver is a regression oracle, not an independent adoption oracle.
For an analytic reference, use a bib_key from the single canonical bibliography
(bibliography_canonical_path). Its entry needs a DOI or checksum-valid ISBN-13;
a URL alone is not sufficient. Keep the hashed symbolic_verification receipt.
Do not invent a DOI for a textbook or relabel the reference to bypass this rule.
Identifier checks do not verify the reference's mathematical adequacy.
Obtain actual user approval; an agent or reviewer cannot approve on the user's
behalf. Preserve the dated consent identifier. Public approval metadata contains
only date and identifier; private consent locations remain outside public data.
Do not change an existing gate, including strict versus inclusive comparison,
without explicit approval. Never backdate approval to a historical execution.

validation_test executes the approved criterion. Save values, PASS/FAIL, exact
source and artifact hashes before post-processing. Put the criterion, approval
and result in the validation record, then link its ID from the lock's Assumes:.
A passing execution without approved criteria is evidence, not adopted status.
A record containing approval text does not authenticate consent: independent
review must check the actual consent record and its scope. Check the recorded
threshold against validation code mechanically when possible; otherwise mark
reviewer confirmation required, never silently call it machine-verified.

radia-mcp and docs may describe a feature as validated only when citing its
applicable validation record. Existence alone is insufficient: check PASS,
user approval, exact declared scale/backend/reference and current/stale status.
A failed, stale, pre-rewrite or unapproved result cannot establish a current
verified claim. Unknown/missing links stay unverified. Do not generalize one
stencil to a whole interface, one linear case to nonlinear physics, or CPU to
GPU. The documentation/catalog projects the record; it is not another source
of adoption state. Runtime is metadata unless separately approved as a criterion.
"""),
    "layers": ("14. Protection, adoption and champion evidence", """
Three layers, with different execution purposes:

 tests: convergence level1 -> saved-self surrogate test -> every PR
 validation_test: approved convergence series -> adoption -> trigger/release
 validation_records: beyond adoption -> finest/highest-p champion once,
                     or speed/scale record -> schema-only checks

A validation record can preserve the adoption decision; a champion references
that decision by ID and hash. An accuracy-champion continues the sequence to
finer meshes/higher p, retaining the complete table through the selected finest
level, each level's source SHA/artifact hash/error/runtime, declared conditions
and limits. Speed and scale records preserve their accuracy and comparability
conditions; a measured timing alone does not prove a faster solver in general.

Use supersedes to name the previous champion and explain improves/corrects/
replaces. Keep old evidence, even stale or superseded. Covered-code changes
trigger a stale mark and reconsideration, not deletion. Do not schedule champion
reruns: one-off heavy work beyond 30 minutes is not periodically repeated or made
a PR gate. Routine recertification runs the adoption series in validation_test.
A surrogate test is the lowest adopted level, not the finest champion and not
a new comparison with the analytic reference. Runtime is metadata, never a time
assertion. No champion, record or registry fabricates user approval.
"""),
    "inference": ("13. Judgment inference rules", """
Record each judgment as premises -> scoped conclusion, known failure modes and
recheck action. These are defeasible engineering inferences, not logical proofs
of universal correctness. Preserve the reasoning, not merely a green label.

1. CI green + qualified shared precheck + independent review -> eligible for
merge within the reviewed scope. A no-op step can be green without checking its
claim (reported PowerPoint step and two policy audits). Recheck commands,
selected files, actual execution receipts and an injected failure; status alone
is insufficient. GitHub execution remains the assigned executor's action.

2. Production mutation killed -> this lock is sensitive to that particular
fault on the executed path. A mutation in a public wrapper survived because the
solve followed another path. Recheck actual entry/coverage, mutation patch and
restored PASS; revise Guards or add a sensitive quantity/contract. Do not loosen
tolerance, or infer coverage of every fault from one killed mutation.

3. Convergence to an analytic reference at the expected rate -> evidence that
discretization is consistent under the declared regularity, geometry, norm and
series conditions. A copied reference formula can be wrong: check it with
symbolic identities and a primary source. Kelvin P2->P3 also changed geometry
order, so that historical pair cannot attribute improvement to p alone. Recheck
reference, fixed-axis series, geometry and observation integration separately.

4. Independent solvers agree -> increased confidence within their common tested
scope. Shared mesh, integration points or constitutive inputs can carry the same
error into both. Recheck the independence inventory and shared components with
analytic identities or another genuinely independent reference. Agreement is
not proof that both are correct.

5. A surrogate test agrees with past self -> the adoption-era numerical behavior
is preserved within measured tolerances and assumptions. Many individually
small changes can accumulate; dependency upgrades can invalidate assumptions.
Recheck cumulative deviation against the approved baseline, dependency identity
and trigger a fresh validation before intentional baseline regeneration. A
surrogate test does not independently recertify the physical formulation.

6. A matched Kelvin-on/off comparison can support the open-boundary treatment
on its declared case. APPROVAL-20261023-KELVIN-SPHERE-01 approves a test plan
including integral-energy error <=1% with Kelvin and >=5% for its truncated
control. Approval of that criterion is not a measured result or an adoption
verdict; do not report the thresholds as achieved errors. Recheck
matched mesh/source/material/stencil/residuals, boundary placement and reference;
obtain user-approved criteria before adopting the conclusion. Do not substitute
these illustrative numbers for any existing gate.

For each adoption record store a scoped inference and assumptions with stable
IDs (record ID, assumption ID), evidence IDs, invalidation conditions and actions.
Include Assumes, coverage, independence and reference verification explicitly
where relevant. When an assumption fails (e.g. survived mutation, wrong analytic
formula), traverse the IDs to all dependent judgments and reconsider them. A
matching JSON/hash authenticates neither consent nor scientific truth.
For adopted features, aim for two genuinely independent detectors per
production mutation, e.g. a surrogate test plus validation. Record each
detector failure and restored PASS. Repeats and shared coverage do not prove
independence. A coarse surrogate test may differ from an analytic field by1%
or more: do not enlarge it or add an analytic accuracy gate for that reason.
It detects change, supported by mutation and adoption links. Keep
private receipts private; public data carries identifiers, dates and hashes.
"""),

}


def get_tests_design_guide(section: str = "overview") -> dict:
    """Return one authoring section; unknown selectors fail without side effects."""
    if not isinstance(section, str) or section not in SECTIONS:
        raise ValueError(f"Unknown tests-design section {section!r}; expected one of {list(SECTIONS)}")
    title, text = SECTIONS[section]
    return {
        "section": section,
        "title": title,
        "content": text.strip(),
        "sections": [{"id": key, "title": value[0]} for key, value in SECTIONS.items()],
        "policy_source": "CLAUDE.md: Tests / Tests For AI-Authored Code",
        "related_tools": ["radia_mcp_golden_gate", "bug_patterns_lookup"],
        "helpers": ["tools/regression_lock.py", "tools/coverage_difference.py"],
        "read_only": True,
    }
