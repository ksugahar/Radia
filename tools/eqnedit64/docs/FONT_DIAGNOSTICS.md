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
