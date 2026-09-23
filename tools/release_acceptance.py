"""The one place the release acceptance hosts are named.

Release-quad requires LAB, 100号機, mdx1 and mdx2 for the same release commit
(AGENTS.md / CLAUDE.md, "CI Execution, Validation Evidence, and Notebook
Policy"). hibino is a computation host and is not an acceptance target.

Two tools consume this: `release_quad.py`, which produces the evidence, and
`verify_radia_promotion.py`, which refuses to publish without it. They used
to each carry their own list, and drifted -- the promotion gate still asked
for the two-host quad it was written against, so a release with two of the
four hosts' evidence would have passed the public gate. Importing the tuple
from here is what makes that drift impossible rather than merely unlikely.

Deliberately dependency-free so the promotion verifier, which runs in a
GitHub-hosted workflow, can import it without pulling in the release
tooling's SSH, MATLAB and filesystem assumptions.
"""

# Evidence directory names under
# validation_test/<suite>/results/candidate_<sha9>/<host>/, in gate order.
RELEASE_ACCEPTANCE_HOSTS = ("lab", "100", "mdx1", "mdx2")

# Evidence directories use deployment roles; JUnit records the OS hostname.
# 100 is an SSH/deployment alias for INTEL11, not a second physical machine.
RELEASE_ACCEPTANCE_HOSTNAMES = {
    "lab": frozenset({"lab"}),
    "100": frozenset({"100", "intel11"}),
    "mdx1": frozenset({"mdx1"}),
    "mdx2": frozenset({"mdx2"}),
}

# Human names, for messages that name a machine rather than a directory.
RELEASE_ACCEPTANCE_HOST_LABELS = {
    "lab": "LAB",
    "100": "100号機",
    "mdx1": "mdx1",
    "mdx2": "mdx2",
}

assert tuple(RELEASE_ACCEPTANCE_HOST_LABELS) == RELEASE_ACCEPTANCE_HOSTS
