# Shared MCP runtime policy

This is the operational contract referenced by AGENTS.md and CLAUDE.md.
It governs development, testing, release deployment and reconnect work on LAB
and 100号機. The maintenance CLI implements cooperative host-wide exclusion and
activation receipts; it does not intercept arbitrary installers or provide
per-client reconnect control. See [maintenance](../maintenance.md) for its scope.

## 1. Separate development from active service

Use a task-owned worktree and virtual environment with an explicitly selected
interpreter for installs and tests. A process-local PYTHONPATH override is
acceptable for a disposable probe when dependencies are already available;
check the resolved module path and never persist that override in a shared
client or user environment. Neither probe verifies an existing live client.
Do not use global/user-site pip install, uninstall, manual .pth edits, symlink
retargeting or shared MCP configuration edits to make ordinary tests pass.

For long-lived service, retain editable support but pin the source to an
approved commit in a dedicated snapshot. A snapshot name containing a SHA is
not proof: verify the actual Git commit and working-tree state. Complete builds
before activation and record native artifact hashes where relevant. Do not
edit code, advance Git HEAD, replace binaries or delete that source while any
consumer still depends on it, including consumers on another host via the NAS.

Prefer a dedicated, explicitly configured runtime interpreter for each service
generation so later shared site-packages changes cannot affect lazy imports.
Migration to such interpreters is a separately scoped deployment, not an
automatic action authorized by this policy. Until migrated, freeze the shared
interpreter as well as its active source while its consumers remain in flight.

## 2. One deployment owner, not one owner per chat

Before any shared mutation, identify the host, resolved interpreter, affected
distributions, users, clients and server names. Normalize mapped-drive and UNC
aliases when identifying the same resource. A LAB install is not evidence for
100号機, and an administrator's client is not evidence for every user's client.

Select one owner task for that complete deployment scope and record the handoff
in a shared maintenance record at a path all participating operators can read.
Declare that path before starting. Keep it outside Git and outside any snapshot
that may be retired; C:\temp is permitted for this operational record under the
lab temp policy, provided it is accessible and retained through the handoff.
An inaccessible record is not coordination. A second task must hand off to the
owner, not perform another install, restoration or reconnect concurrently.

Use an existing exclusive deployment lock when provided by the workflow. A Git
index lock or a release receipt's file-write lock is not a deployment lock.
If no deployment lock exists, explicit owner assignment and acknowledgement
from affected operators is required before mutation; if this cannot be
established, defer. Do not claim a prose record provides atomic exclusion.
Never steal ownership because a PID is idle, a lock is old, or a task is quiet.
After an interruption, reconcile the actual state before accepting a handoff.

Required record fields:

- Change ID, authorization/scope, owner task and named handoff recipient.
- Host, interpreter, affected users/clients/servers and distributions.
- Previous and proposed source roots, actual commits, package versions and
  native artifact identities where applicable.
- Reason, expected tool additions/removals/schema changes and test evidence.
- Observation and operation timestamps, including timezone.
- Staging/activation outcome, per-client verification results and pending
  targets, rollback source and conditions, and final owner handoff.

Store no credentials or private client tokens in this record. Record failures
as well as successes; do not erase earlier evidence when another task takes over.

## 3. Stage, coordinate, activate, verify

1. Read the record and live status before deciding that an install has drifted.
   The approved deployment target, not whichever checkout is newest, defines
   expected state. A comparison with origin/main is an advisory freshness check;
   being older can be an intentional pin, and being newer is not authorization.
2. Build/test the candidate in isolation. A main merge or a package publication
   does not switch any live service. A radia-mcp-only change does not repoint
   radia or cubit-mesh-export, and their releases do not silently repoint MCP.
   A coordinated cross-package deployment must name and validate all of them.
3. Check in-flight work and state ownership for the complete activation scope.
   Keep the old environment available. If the change mutates an interpreter
   still used by busy clients, do not change its registration first and promise
   to reconnect them later: defer or stage a separate environment instead.
4. The designated owner performs the approved activation once and records its
   result. Reconnect only authorized idle targets using supported client
   controls and the mcp-reconnect skill. Do not kill all Python/client processes
   or reload an unrelated CAD/MATLAB stateful service to force convergence.
5. Through each original client, inspect live discovery, loaded-source evidence
   and a harmless changed tool. A fresh subprocess, module import or successful
   tools/list alone is not full live acceptance. Record busy, unreachable and
   manually pending clients individually, without retry loops.
6. Rollback is another coordinated activation, not an automatic global repair.
   Do not remove the old snapshot until every known dependent has retired and
   unknown consumers have been resolved. Retaining it is safer than guessing.

Release completion is not permission to restore an older dirty canonical tree.
Neither a build failure nor a metadata mismatch alone authorizes restoration.
Use the approved source and ownership contract even if a legacy script offers
a broader restore command. Conflicting legacy automation must not be run as-is.

## 4. Evidence and reporting

Keep these observations distinct:

| Observation | What it establishes | What it does not establish |
| --- | --- | --- |
| .pth / direct_url / installed metadata | Registration for an interpreter | Code already loaded in an existing process |
| Fresh-process import and source path | What that new process resolves | Another client's loaded code or schema |
| Live loaded-module path and registration-time identity | Evidence about that target process | Every domain module or every other client |
| Live tool discovery and harmless changed call | That target exposes and executes the checked tool | Global adoption or whole-package identity |

A status response may mix new installed metadata with an old module path. Do
not choose one field and label the whole response current. Matching package
versions do not prove matching commits; a current on-disk hash is not a hash
of already loaded Python code. Check registration-time identities for affected
modules where exposed, plus changed behavior; missing provenance stays unknown.

Report each target with its host/user/client/server identity and timestamp as
`verified`, `deferred-busy`, `manual-action-required`, `failed` or `unverified`.
Use `unverified` with a `mixed-generation` reason for contradictory evidence.
A deferred result expires as an observation: re-read before repeating it in a
later handoff or final answer. One verified target never promotes others.

Always distinguish source integration, installation/activation, live code and
client schema. "This task did not deploy" describes an action boundary, not
the current global state. Do not repeat "live reflection remains deferred"
unless that named target has been freshly checked. For code-only work say
"no runtime changes made by this task; live state not checked" instead.

## 5. Enforcement follow-up

The policy and its synchronized documentation regression apply immediately.
They do not make existing sessions restart or older checkouts read new rules.
Provide this policy and the current deployment record at the next owner handoff.
Do not patch another task's dirty checkout just to distribute instructions.

The maintenance CLI now provides exclusive host-wide ownership, compare-before-
change checks and atomic receipts for config edits and radia-mcp-only editable
activation. Tests cover competing processes, interrupted work and same-version
source changes. Idle status is operator-attested and client results remain
unverified; this is not unattended deployment acceptance.
Remaining work includes integrating legacy release entry points, cross-host
coordination, automatic busy-client discovery and per-client generation evidence.
No new watcher, service or client restart is installed by this implementation.
