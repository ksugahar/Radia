# Read-only session font diagnostics

`doc_convert_session_font_check(faces="Arial, Times New Roman")` inspects the
requested faces through GDI in the server process/session. It reports the
realized face, substitution, font-data and glyph-query results. Empty requests,
substitutions, query failures and cleanup failures must not pass.

The implementation declares pointer-sized Win32 handles and restores/deletes
only the logical font it created, then releases its own DC. It never changes
font registrations, unloads session-shared font references, stops processes or
broadcasts notifications. If restoring the old selection fails, it reports the
cleanup failure instead of deleting a potentially selected font.

The CLI is check-only:

```powershell
python -m radia_mcp.doc_convert.plans.T17_session_fonts --faces Arial --crash-window 60
```

The optional bounded event-log query reads XML. Command errors, invalid output,
timeouts and truncated coverage remain unknown and cause a nonzero CLI result.
Application-log events are not session-specific proof, and an observed crash
does not establish the cause of a glyph-query failure.
When `--crash-window` is explicitly requested, both an observed crash and an
unknown query cause nonzero exit. `face_ok` reports only the face checks;
aggregate `ok` additionally requires a successful query with no matching crash
in that window. This does not certify accident-free operation outside that
window or attribute an observed event to the checked face or current session.

The shared WIP's `--repair`, `--force`, `--notify`, registration/unregistration,
automatic-repair and `msg *` behavior were deliberately not recovered. Existing
scheduled repair commands must not be repointed to this module: unsupported
options fail closed. This change does not edit any existing scheduled task.

Acceptance for this recovery is mocked ABI, ownership, failure and MCP registry
tests only. No LAB/100 native font API, fontdrvhost workload, font registration,
Office export or session-repair validation was performed. A successful mocked
test does not certify an actual font or another client's document export.
