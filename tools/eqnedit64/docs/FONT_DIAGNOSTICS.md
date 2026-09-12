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
