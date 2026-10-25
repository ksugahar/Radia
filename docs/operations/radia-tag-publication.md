# Radia exact-tag publication

Dispatch Release on the new release tag ref (`v<version>`), never on main.
Supply the successful native tag run ID, the accepted wheel SHA-256, and the
main-reachable acceptance-data commit. The dispatch tag must resolve to that
native run's source commit. No wheel or native binary is rebuilt here.

The read-only promotion audit retains both host acceptance gates and native
producer CI checks. It freezes wheel-identical native binaries and binds the
publishing repository, tag/source SHA, run ID and attempt into its receipt.
Immutable accepted artifacts are downloaded by artifact ID after approval.

Only publish-pypi uses the pypi environment. After user approval it rechecks the
live remote tag, publishing attempt, input identities and accepted bytes before
publication. Existing versions fail. The native binaries upload needs successful
publish-pypi and repeats the identity and byte checks; it requires no second
approval. Checkout credentials are not persisted. No deployment-approval API is
part of this route.

A publish-only rerun changes run_attempt and rejects the old receipt. Recover by
rerunning all jobs and obtaining fresh approval. If PyPI succeeded but the later
native upload failed, stop and ask repository management to repair the partial
publication; do not republish or bypass the checks. Tags predating this workflow
are unsupported. Tag creation, publication and deployment require their separate
user approvals and are performed by repository management.
