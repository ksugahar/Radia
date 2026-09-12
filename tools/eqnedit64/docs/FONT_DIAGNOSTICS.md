# Font-host incident diagnostics

Release remains on hold. The crash in isolated CI run `34657957211` matches
the earlier Server 2022 signature: fontdrvhost 10.0.20348.5256, c0000005,
offset 366a2. Its event timestamp lies in the three hidden executable checks,
but delayed reporting means the responsible command/API is not established.
Passing 32 earlier lifecycle iterations does not exonerate later operations.

Set `EQNEDIT64_FONT_TRACE_DIR` to an existing, dedicated artifact directory
**before the first native process** in disposable CI. This enables append-only
`font-<PID>.jsonl` records. Collect them even when tests fail. Each line contains
UTC millisecond time, process/thread/session IDs, a fixed stage name and numeric
result. No source text, clipboard data, or user document paths are recorded.
Unset the variable for normal use. No font registration or removal policy is
changed by this instrumentation, and it is not a fix for the incident.

Stages: `app.enter`, `app.return`, `cache.begin/end`, `register.begin/end`,
`measure.begin/end`, `metrics.destroy.begin/end`, `drawfonts.destroy.begin/end`.
Registration results are the AddFontResourceExW count; measurement results are
boolean. Measurement begin reports the retry index. Application return precedes
static cache destruction; it is not proof of process termination. Pair PID/time
with the external command interval, font-host PID history and WER events.
Missing records, write failure or absent terminal records are not success.
The trace preserves thread last-error state; opt-in disk flushes perturb timing,
so non-reproduction with tracing does not establish safety.

`test_font_trace` tests disabled behavior, append records and last-error
preservation without loading a font or linking the renderer. The actual
instrumented executable must not be run on shared LAB/100 sessions for this
investigation. The gate/evidence workflow is maintained separately in PR #205.

Microsoft documents matching flags for Add/RemoveFontResourceExW. There is no
evidence yet that adding an explicit removal fixes this crash; do not introduce
one as an unverified remedy or claim private registration isolates the host.
Reference: https://learn.microsoft.com/en-us/windows/win32/api/wingdi/nf-wingdi-removefontresourceexw

## API and dump documentation check

Microsoft's AddFontResourceExW contract explicitly provides no extended error
information on failure. Do not report GetLastError as the registration failure
reason. Count returned-zero calls, retain timestamps, and correlate external
host events/dumps. FR_PRIVATE resources are removed by Windows at process
termination; private visibility does not imply a private font-driver host.
https://learn.microsoft.com/en-us/windows/win32/api/wingdi/nf-wingdi-addfontresourceexw

WER LocalDumps supports per-executable configuration. Its output directory must
grant write access to the crashing process identity, not just the CI account.
Configure only the disposable CI machine, retain dumps privately and briefly,
and regard an empty dump directory as missing evidence rather than no crash.
https://learn.microsoft.com/en-us/windows/win32/wer/collecting-user-mode-dumps

Exact searches for fontdrvhost with offset 366a2, version 20348.5256, Latin
Modern, and AddFontResourceEx did not identify a matching published fix. This
negative search result does not establish that no upstream defect exists.

## Independent trace reading: run 34658678863

Source: PR #205 at `775a03984`; executable SHA-256
`E883E56313593A6072CE295B7077936211D49BB3BC1221994DFAB287F9984DD5`.
Evidence is retained under `C:/temp/eqnedit-font-acceptance-34658678863`.
All times below are UTC on 2026-09-11.

| Time | Observation |
| --- | --- |
| 23:38:29.298 | Model process 6572 registers successfully. |
| 23:38:32.597 | Model metrics-cache destructor begins and ends. |
| 23:38:32.972 | Self-test process 3596 enters; first registration returns zero in the same millisecond. |
| 23:38:33.3632636 | WER event records fontdrvhost PID 2920 crash (c0000005/366a2). |
| 23:38:33.475 | Ninth failed registration returns zero. |
| 23:38:33.538 | Registration succeeds; only now does the first successful measurement run. |
| 23:38:33.757 | Self-test returns and both cache destructors complete. |

The model destructor precedes self-test entry by 375 ms. The recorded self-test
measurement and its cache destruction occur after the incident event; they
cannot explain it as synchronous preceding operations. Registration itself is
not exonerated: the first call could trigger failure, or encounter a host already
failing after model-process teardown. WER time is reporting time, not a precise
death timestamp. Nine failed calls do not prove nine crashes.

Next controlled comparison uses the identical binary on independent fresh CI
machines: self-test alone; model suite followed by idle observation with no next
EXE; model suite immediately followed by self-test. Capture host process
liveness from before the first operation through the idle tail, along with
command intervals and JSONL. Keep diagnostic results separate from release
acceptance, and do not weaken the existing failure gate or add retry-until-green.

## Offline dump analysis: run 34661533840 (2026-09-12)

Two independent disposable windows-2022 jobs using the same binary payload
captured full dumps. Both model and self-test returned zero despite the host
crash. Dumps and debugger logs remain ACL-restricted local evidence, not public
repository artifacts.

Microsoft CDB with the public PDB matching the dump's CodeView GUID
`DE2401FD-C131-8258-9A23-A072EE4125B1`, age 1, resolves both crashes to
`fontdrvhost!CheckMemoryUsage+0x26` (module offset `0x366a2`). Both are reads
of address `0x8`, with RDX zero, in fontdrvhost 10.0.20348.5256.

The matching call chain is:

```
DrvLoadFontFile -> LoadFont -> AllocateFontCollection
               -> ATMallocExt -> CheckMemoryUsage
```

Both dumps have identical relevant global state: the first two DWORDs at
`atm` are `0x0053e134, 0`; `lruListHead` is null. Disassembly enters the
list-processing path when the first DWORD exceeds `0x500000` or the second
exceeds `0x0fa00000`, then dereferences `[lruListHead+8]` without a null check.
The first condition holds in both dumps. This establishes the immediate fault
mechanism, not how this inconsistent management state arose. It is not evidence
of system-wide out-of-memory, nor proof that malformed font data is responsible.

The symbol server image emits a checksum mismatch warning; do not omit it.
CodeView GUID/age match, and `!lmi` reports the image read from dump memory and
public PDB symbols loaded successfully. Preserve original dumps for independent
confirmation rather than treating the downloaded image as the original binary.

Next isolate which registration/measurement/process-lifetime sequence leaves
this state, using disposable workers and retaining external host monitoring.
Do not patch the Windows binary, weaken the release gate, or claim that an
explicit RemoveFontResourceEx call or retry delay is a fix without controlled
evidence. Product release remains on hold.

## Connection to laboratory evidence

The laboratory's historical CLAUDE.md section "Do Not Break the Interactive
Session's Fonts (2026-09-03)" records 202 incidents since August 23 at offset
0x366a2. The reviewer reports September 5 symbolization of 203 LAB crashes at
the same CheckMemoryUsage location; that later count is reviewer-provided, not
a new recount of the original event logs in this investigation.

Independent non-CI evidence is recorded in
`C:/temp/eqnedit64-palette-claude-review/REVIEW.md`, section "fontdrvhost crashed
in 100's session 0": September 11 at 14:21:11 JST, ja-JP session 0, host
10.0.20348.5256, c0000005/366a2, coincident with self-test startup. The review
reports no competing Eqnedit64/python/node process in that session and only
this event in the preceding 30 days. Both application checks returned zero.
This makes disposable CI or en-US locale unnecessary conditions for the
observed signature; it does not establish identical preceding host state.

## Count-controlled diagnostic (not a remedy)

The identical `atm` value in two independent dumps motivates a history-dependent
accounting hypothesis. A retained counter after list teardown is one possible
explanation, not an established leak. Deterministic counts can also arise from
deterministic timing; variable counts do not alone prove a race. The threshold
0x500000 is 5 MiB (5,242,880), not decimal 5 MB.

`font_lifecycle_probe` is an EXCLUDE_FROM_ALL console target with no editor or
renderer linkage, no retries and no cache writes. Build it explicitly. Its CLI:

```
font_lifecycle_probe FONT exit|remove|hold none|measure HOLD_SECONDS
```

HOLD_SECONDS is 1..60. Only hold mode sleeps. Each process registers once;
remove mode explicitly removes once with matching flags; exit mode relies on
Windows teardown. Optional measure mode selects Latin Modern Math, verifies
the resolved face and measures one parenthesis before deleting its GDI objects.
Hold mode keeps the registration alive until its bounded timeout. JSONL includes
PID, precise FILETIME and operation boundaries/results, but not user contents.
Both GITHUB_ACTIONS=true and EQNEDIT64_ISOLATED_TEST_SESSION=1 are required;
these flags are a guard, not proof of actual machine isolation.

Use one immutable probe/font payload on fresh hosted workers for each arm and
replicate. Record hashes, OS/font-host version, baseline PID/start time, child
intervals, events and first failing k. Start with exit/remove crossed with
none/measure, three workers per arm, at most 64 sequential processes each.
Observe a fixed idle interval before/after, stop at the first host death/change
or child error, and report NO_REPRODUCTION_WITHIN_BOUND separately from success.
The hold arm needs overlapping children to distinguish retained registrations
from process teardown; cap concurrency explicitly and terminate only owned
children after observations. Do not label an ordinary sequential hold loop as
this control. Preserve timing as a factor rather than replacing order tests.

The bundled font license is GUST Font License (LPPL), not OFL; see
`assets/GUST-FONT-LICENSE.txt`. A CFF-to-glyf conversion is a separate diagnostic
candidate, not an approved replacement. Preserve provenance and applicable
license/renaming conditions, validate changed tables and outline approximation,
and test rendering separately. Failure of CFF but not glyf within a finite
bound supports a route-specific hypothesis; it does not prove ATM is a necessary
condition for every incident. No Windows binary patching is contemplated.

### Count-controlled result: 34663488070

All 12 fresh-worker trials (exit/remove crossed with none/one-glyph measurement,
three replicates each) completed 64 successful registrations: 768 launches in
total. No child, measurement, removal, host-exit/PID-change, WER or monitor
query error was recorded. All report NO_REPRODUCTION_WITHIN_BOUND. Evidence:
`C:/temp/eqnedit-font-probe-34663488070/summary.json` and adjacent artifacts.
Probe SHA256: BFE93F6EDF8F1A4D397974844C3F1E9D491441046E0F6F7CD5D84FE96CA22146.
Font SHA256: 6075562B771F8B82F0C179E363389684F2DD09DE30038269E2628E504BD7BE0F.
Host StartTime was not recorded in this round; retain that limitation.

This does not support adopting explicit removal as a remedy: neither control
nor removal reproduced. Registration count alone under these conditions is
insufficient to trigger the known crash. It does not exclude accounting driven
by a richer preceding drawing history.

Next use `run_model_tests.py --diagnostic-prefix N`, followed immediately by
the fixed reproducing self-test EXE, on fresh workers. Prefix zero still imports
the module and loads/measures the embedded font; prefix ten is the full known
reproducing workload. Prefixes 0,3,5,6,8,10 with two replicates each localize the
required history without changing the default acceptance suite. Reuse the same
EXE/pyd payload across arms and record runner and payload source separately.
If the full-prefix control does not reproduce, do not interpret shorter prefixes
as exoneration. Query failures are INCONCLUSIVE; record host StartTime as well
as PID. Product release stays on hold.
