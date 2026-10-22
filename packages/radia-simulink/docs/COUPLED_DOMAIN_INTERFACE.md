# Coupled domain host contract

This interface owns a single coherent accepted state: Eddy solution, heat samples,
temperature coefficients, diagnostics, revision, step count and time. Language
adapters retain immutable host snapshots and own checked handles, generation,
sample identity and serialization through the common runtime policy.

## Current execution boundary

The config is the existing `linear_temperature` / `independent_cpu` contract.
Structural validity does not imply an executable physical provider. The canonical
production factory currently validates config and raises `unsupported_physics`.
Its diagnostic names the requested backend and this document. No initialized
domain or solved output is published after failed construction.

`radia_ih_coupled_domain` is a **production-source link boundary check in the test
tree**, not a shipped native solver target. It links only `ih_coupled_domain.cpp`
and `ih_coupled_backend.cpp`. The rejection test links this target. Transaction
scenarios define synthetic physics in a separate executable and must never be
linked into production. Existing MEX build inputs are unchanged.

No GPU backend, resolved nonlinear cycle, physical solver implementation, release
scope change or CPU/GPU numerical qualification is introduced by this contract.

## Atomic advance and reset

1. Validate config, controls, counters and the simultaneous accepted/trial memory
   budget before work. A future provider must account for all physical storage
   and workspace, not just the published vectors.
2. Stage all fallible operations on independent trial storage. Accepted Eddy,
   heat, temperature and diagnostics must not alias mutable trial workspace.
3. Finish solve/projection operations, acceptance gates, synchronization and host
   snapshot capture before publication. The host owner checks complete shapes,
   finite values, reports, counters and capacities before commit.
4. Commit uses only noexcept ownership swaps of physical state and the complete
   immutable host snapshot. Retired state remains in the trial and is destroyed
   outside this switch. Reset uses the same staging rule and configured initial
   temperature with zero heat/Eddy state.
5. A precommit failure preserves both identity and contents of the last coherent
   snapshot. A device-failure classification poisons future computation but
   preserves host observation. Close is idempotent; previously retained snapshots
   remain valid after close. No implicit replacement backend is permitted.

## Host evidence and limits

The CTest scenario executable injects failures after trial modification at Eddy,
heat projection, thermal solve, power gates, synchronization, host capture and
reset. It also checks capacity budgets, malformed snapshots, poisoning, owner
sample idempotency and reset generation. These are source/lifecycle checks with
synthetic data, not physical validation or proof of rollback in an unimplemented
physical provider.

The separately linked production-unavailable executable checks the real factory
and repeated failed domain construction. Neither suite establishes numerical
accuracy, independent FEM agreement, GPU operation or installed-artifact readiness.
