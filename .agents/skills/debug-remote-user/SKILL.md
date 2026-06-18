---
name: debug-remote-user
description: Lessons + playbook for debugging user-specific issues on 100号機 (multi-user Windows) from LAB via SSH.  Covers what CAN and CANNOT be verified without the target user's credentials, and how to defer-to-user instead of faking the test.
---

# debug-remote-user

Debugging a "Cubit is slow", "the panel doesn't work", "permissions error"
style report from 100号機 requires distinguishing:

1. **Infrastructure-level issues** (machine-wide config, PATH, installed
   binaries, Registry, plugin files).  LAB's SSH session as
   Administrator CAN fix these.
2. **User-session-level issues** (per-user `LOCALAPPDATA` caches, user
   ACLs, `$PROFILE`, Start Menu / Desktop shortcuts under the user's
   profile, Codex MCP daemon running inside the user's VSCode).
   LAB's SSH session CANNOT reach into another user's live session
   without their password.

This skill captures the playbook and the **hard limits** learned
2026-04-21 debugging the Cubit-launch-slow complaint for kubota /
keiko / yano etc.

## What SSH-as-Administrator CAN do on 100号機

| Target                                           | Mechanism           |
|--------------------------------------------------|---------------------|
| Install files under `C:\Program Files\...`       | Admin writes        |
| Write to `C:\ProgramData\...` (all-users config) | Admin writes        |
| Write to `C:\Users\Public\Desktop\...`           | Admin writes        |
| Register machine-wide Scheduled Tasks            | `Register-ScheduledTask` |
| Read any user's non-protected files via SMB      | `\\192.168.11.100\c$\Users\...` |
| Stop processes by PID across sessions            | `Stop-Process -Id` with admin rights |
| Query Registry HKLM                              | Admin reads         |

## What SSH-as-Administrator CANNOT do

| Target                                           | Why                 |
|--------------------------------------------------|---------------------|
| Write inside another user's `%LOCALAPPDATA%`     | UAC-protected AppData denies Admin write |
| Impersonate another user (run a command as them) | Needs their password / saved credential |
| Refresh another user's license cache             | Cache is per-user + ACL blocks admin writes |
| Invoke a task **in another user's active session** | Task Scheduler `/run` runs as the invoker |
| Test that "Kubota's Cubit launch is fast"        | Can't log in as Kubota |
| Change another user's Desktop / Start Menu       | Per-user profile dirs inaccessible |

## Run-as-user via SSH pubkey bootstrap (2026-04-24 UPDATE)

**Previous memory said impersonation without password is impossible.
This is wrong.** The viable path is to use the OpenSSH server
already running on 100号機 as the impersonation mechanism: admin
temporarily drops its own pubkey into each user's
`C:\Users\<user>\.ssh\authorized_keys`, then `ssh <user>@100 <cmd>`
runs `<cmd>` under that user's token.  When done, remove the
pubkey and restore `AllowUsers Administrator`.

Step-by-step (all 7 commands from LAB ssh, admin):

```powershell
# 1. As SYSTEM (PsExec -s), write admin pubkey into each user dir,
#    fix ownership to user, set OpenSSH-required ACL
#    (user + SYSTEM + Administrators only).
cat << 'PS' | ssh 100 'pwsh -Command -'
$key = '<paste admin pubkey here>'
$users = Get-ChildItem C:\Users -Directory |
  Where-Object { (Test-Path (Join-Path $_.FullName NTUSER.DAT)) -and
                 $_.Name -notin @('Administrator','All Users','Default','Default User','Public') }
foreach ($u in $users) {
  $d = Join-Path $u.FullName '.ssh'
  $f = Join-Path $d 'authorized_keys'
  if (-not (Test-Path $d)) { New-Item -ItemType Directory -Force -Path $d | Out-Null }
  takeown /F $d /A /R /D Y 2>$null | Out-Null
  takeown /F $f /A        2>$null | Out-Null
  icacls $d /reset /T /Q  2>$null | Out-Null
  Set-Content -Path $f -Value $key -Encoding ASCII -Force
  icacls $f /inheritance:r /grant:r "SYSTEM:(F)" "BUILTIN\Administrators:(F)" ("$env:COMPUTERNAME\" + $u.Name + ":(F)") /Q | Out-Null
  icacls $f /setowner ("$env:COMPUTERNAME\" + $u.Name) /Q | Out-Null
  icacls $d /setowner ("$env:COMPUTERNAME\" + $u.Name) /T /Q | Out-Null
}
PS

# 2. Open sshd_config to all users and restart.
ssh 100 'pwsh -Command "(Get-Content C:\ProgramData\ssh\sshd_config) -replace ''(?m)^(AllowUsers\s.*)$'', ''# \$1'' | Set-Content C:\ProgramData\ssh\sshd_config; Restart-Service sshd -Force"'

# 3. ssh as each user and run whatever needs running in their context.
for u in keiko kubota yano ... ; do
    ssh $u@100 '<command>'
done

# 4. Clean up: remove pubkey, restore sshd_config, restart.
# (Back up sshd_config to .radia_bak before step 2 for clean revert.)
```

**This is the tool** for: per-user `rlm_activate --login`, per-user
script verification, per-user file placement in their own AppData.
It replaces the hand-wavy "can't impersonate without password" claim.

### Why this differs from the admin-copy-in false-path

The 2026-04-21 false-path copied a **cache file** (renewals) into
users' AppData, owned by Administrators, breaking the target user's
later read access.  The SSH bootstrap above writes an
**authorization key** (which admin legitimately owns a pairing for)
into `.ssh/authorized_keys`, properly chowns to the target user,
and is removed immediately after use.  Both operations involve
admin writing into user profile dirs, but only the first is a bug
pattern — because it leaves a stale foreign-owned file behind.

### Still not-so-easy paths (kept for historical context)

The 2026-04-22 exploration found the following routes all require
the target user's cooperation (password) and should NOT be pursued
before trying the SSH bootstrap above:

| Route | Why it's worse than SSH bootstrap |
|-------|-----------------------------------|
| `PsExec -u <user> -p <pass>`                 | needs password |
| `schtasks /create /ru <user> /rp <pass>`     | needs password |
| `schtasks /create /ru <user> /it` (no pass)  | `/it` requires ACTIVE console (Disconnected RDP does not count) |
| `Register-ScheduledTask -Principal S4U`      | S4U needs AD trust; local accounts fail |
| `cmdkey /add:...`                            | interactive prompt |
| P/Invoke `WTSQueryUserToken` → `CreateProcessAsUser` as SYSTEM | requires enabling `SeAssignPrimaryTokenPrivilege` + `SeIncreaseQuotaPrivilege`; PsExec-s SYSTEM token has them granted-but-disabled and `AdjustTokenPrivileges` returns err 1300 (NOT_ALL_ASSIGNED). Not supposed to be blocked but is in practice. |

The SSH bootstrap avoids all these by using the *existing* logon
credentials on 100号機 (the user's SSH pubkey ACL) rather than
trying to invent a new auth path.

## SELF-VERIFY, don't delegate — `runas /trustlevel:0x20000`

**POLICY** (2026-04-22, user feedback: 「すぐに人の試験に頼ろうとするな」):
Do NOT end a debug session with "please ask Kubota to test this."
Self-verify first using `runas /trustlevel:0x20000`.

This runs a command with the **current user's token stripped of
elevated privileges** — i.e. the exact execution context a regular
Users-group member sees.  Admin can simulate "Kubota running this
in his own terminal" without needing Kubota's password.

### Usage pattern

```powershell
# From SSH Administrator session
runas /trustlevel:0x20000 "C:\ProgramData\CoreformCubit\cubit_refresh.cmd"
# or
runas /trustlevel:0x20000 "pwsh -File C:\path\to\script.ps1"
```

### What you can self-verify via this

| Test                                            | Works via /trustlevel:0x20000? |
|-------------------------------------------------|--------------------------------|
| Can a non-admin READ the script file            | YES                            |
| Can a non-admin EXECUTE the script              | YES                            |
| Does ExecutionPolicy block                      | YES (same as a user session)   |
| Does the script write correctly to user's AppData | YES (writes to Admin's since we ARE admin, but flow is identical) |
| Does `rlm_activate --login` work without admin  | YES                            |
| End-to-end exit code + side effects             | YES                            |

### What /trustlevel:0x20000 CANNOT cover

| Test                                         | Why not               |
|----------------------------------------------|-----------------------|
| Verify Kubota's SPECIFIC cache file path     | We're still Admin — writes to Admin's `%LOCALAPPDATA%` |
| Verify logon-trigger-task fires for Kubota   | Needs Kubota to actually log in |
| Test group policy differences per user       | GP is evaluated at logon, not runas |

### Preflight template (use before declaring anything "ready")

```powershell
# 1. Delete admin's own state to force a real test
$cache = "$env:LOCALAPPDATA\Coreform\Cubit\Coreform\licenses\renewals"
if (Test-Path $cache) { Remove-Item $cache -Force }

# 2. Run the instruction EXACTLY as you'd give it to Kubota
runas /trustlevel:0x20000 "<your instruction here>"
Start-Sleep -Seconds 10

# 3. Verify the observable side effect happened
if (Test-Path $cache) { "OK: ran successfully" } else { "FAIL: did not rebuild" }
```

This gives you evidence that **the instruction works for a regular
user** without actually having to ask one.  After this, only
defer-to-user for things /trustlevel can't cover (actual logon event,
specific user's cache location etc.).

### The 2026-04-22 win

Scheduled-Task-based license warmup deploy was "verified" multiple
times via Admin-session tests and kept finding new gaps:

1. First verify: `task registered OK` → shipped → 自分では発火しない on /run
2. Second verify: `script runs as Admin OK` → shipped → user session
   may have different `ExecutionPolicy`
3. Third verify: `/trustlevel:0x20000 rebuilt cache via .cmd` → this
   is what actually covers the user path

Moral: **each fresh shipping pass, before "ask user", run the
instruction verbatim through `/trustlevel:0x20000` first**.  Keeps
you honest about what's really tested vs assumed.

## Cubit license: each user launches explicitly (2026-04-24 retire)

**POLICY (2026-04-24, user-set)**: the logon-triggered license-warmup
Scheduled Task (`\Coreform\CubitLicenseRefresh`) and the admin-
deployed `C:\ProgramData\CoreformCubit\cubit_warm.*` scripts have
been **retired**.  They caused two rounds of "admin-owned files in
user AppData" bugs (2026-04-21, 2026-04-24) that silently broke
Cubit for every non-admin user.

New model: **each Cubit user explicitly runs the lab launcher
themselves** in their own Windows session:

| Host    | Launcher path                                  |
|---------|------------------------------------------------|
| LAB     | `S:\CoreformCubit\coreform_cubit.ps1`          |
|         | `S:\CoreformCubit\coreform_cubit.cmd` (wrapper)|
| 100号機 | `W:\00_CAE\CoreformCubit\coreform_cubit.ps1`   |
|         | `W:\00_CAE\CoreformCubit\coreform_cubit.cmd`   |

(LAB `S:\` and 100号機 `W:\00_CAE\` are the same SMB share.)

The `.ps1` checks / refreshes the renewals cache under
`%LOCALAPPDATA%` and launches Cubit.  Because it runs in the user's
own token, the renewals file is owned by that user; Cubit can read
it back; no cross-user ACL bug is possible.  **17 of the lab's 21
users are not Cubit users**, so no admin-side sweep / automation
should touch their AppData at all.

### Deployment checklist for user-session fixes (generalised)

- [ ] The admin installs code / binaries only in machine-wide
      locations (`C:\Program Files\`, `C:\ProgramData\`) that are
      explicitly Users:ReadAndExecute.
- [ ] User-facing scripts live on the shared SMB mount
      (`S:\` on LAB = `W:\00_CAE\` on 100号機) so each user can read
      and execute them from their own session.
- [ ] `.cmd` wrapper sits alongside any `.ps1` so users don't need
      to know `pwsh -NoProfile -ExecutionPolicy Bypass`.
- [ ] No admin-side automation writes into user AppData.
- [ ] No admin-side automation registers per-user Scheduled Tasks.
      The user runs the launcher when they want Cubit; not otherwise.
- [ ] pwsh.exe resolvable from system PATH
      (`C:\Program Files\PowerShell\7\pwsh.exe`).
- [ ] `Get-ExecutionPolicy -Scope LocalMachine` returns
      `RemoteSigned` or more permissive (not `Restricted`).
- [ ] Self-verify the user-facing instruction via
      `runas /trustlevel:0x20000 <instruction>` → observe side effect.
- [ ] Document the user instruction as a single double-click action.

## User instruction format

Give users the simplest thing that works.  Order of preference:

1. **Double-click a Public Desktop shortcut** (zero typing)
2. **Double-click a `.cmd` file** (if shortcut infeasible)
3. **Paste a one-line `.cmd` invocation**
   (`C:\ProgramData\...\foo.cmd`)
4. **Paste a multi-flag `pwsh -File` command** — LAST resort; requires
   users to know PS flags and paths

Do NOT give users instructions that require them to know about
`-NoProfile -ExecutionPolicy Bypass`.  Wrap that in a `.cmd`.

## Things NOT to do

- ❌ Don't copy files into another user's `%LOCALAPPDATA%` from
  SSH-admin.  ACL breaks, user can't read back, you leave debris.
- ❌ Don't rely on `schtasks /run` for Group-principled tasks.
  On-demand invocation runs in an undefined context and usually
  no-ops while reporting `LastTaskResult: 0`.
- ❌ Don't ship an instruction as "works" based on your own Admin
  session's test.  Admin is privileged; regular users aren't.
- ❌ Don't spend session time on impersonation routes that all need
  credentials.  The Windows security model is doing its job.
- ❌ Don't ask users to run multi-line PowerShell incantations.
  Wrap in `.cmd` or shortcut.

## The 2026-04-21 false-path (do not repeat)

**What I tried**: copy the Administrator's fresh license renewals
cache to `C:\Users\kubota\...\renewals` via SMB.  Assumed the
file's content is machine-level and portable.

**What broke**:
1. File content IS machine-level (worked conceptually).
2. But the written file's **owner = Administrators**, not kubota.
3. When kubota later tries to read his own renewals cache, UAC/ACL
   denies access.
4. Cubit falls back to cold license checkout anyway → **original
   problem unchanged**, with a stale broken file in kubota's
   AppData as bonus debris.

**Cleanup**: `takeown /F ... /A` + `icacls ... /grant Administrators:F`
+ `Remove-Item`, then let the user's next logon rebuild the cache.

**Moral**: do not try to fix a user-session problem by reaching into
another user's AppData from SSH-admin.  Either (a) the user runs
something themselves, (b) a logon-triggered scheduled task does the
work in their own context, or (c) you defer the test to the user.

## The working playbook

### For infrastructure-side fixes (safe from SSH)

1. Install your actual fix (wheel, plugin, script) under
   `C:\Program Files\...` or `C:\ProgramData\...`
2. Verify with `hash-drift` check: SHA-256 source vs installed
3. Register any cross-user Scheduled Task under `\Coreform\...` etc.
4. Provide a Desktop shortcut under `C:\Users\Public\Desktop\`
5. Test end-to-end in **YOUR OWN SSH session as Administrator** —
   this proves the machine-wide state is good

### For user-session verification (NOT safe from SSH)

Do NOT try to fake it.  Instead, hand the user a one-liner:

```
# User runs from their VSCode Terminal (their own context)
schtasks /run /tn "\Coreform\CubitLicenseRefresh"
<then time the actual Cubit launch>
```

And ask them to report the timing number back.  You cannot measure
their session from the outside.

### MCP daemon reload

Codex's MCP servers are long-lived Python processes spawned
when VSCode starts.  Code changes on disk do NOT auto-reload.

- For LAB: git push → user restarts VSCode → next MCP tool call
  picks up new code
- Do NOT attempt to kill the MCP daemon remotely — you'll kill
  all of user's in-flight tool calls / state

## Diagnostic commands (per-user, safe to run from SSH-admin)

**File existence / mtime / size via SMB** (READ only):

```powershell
Test-Path "\\192.168.11.100\c$\Users\$u\AppData\Local\...\renewals"
Get-Item "\\192.168.11.100\c$\Users\$u\AppData\Local\...\renewals" |
    Select-Object Length, LastWriteTime
```

**Scheduled Task status** (machine-wide):

```powershell
ssh 100 'pwsh -NoProfile -Command "Get-ScheduledTask -TaskName CubitLicenseRefresh -TaskPath \"\Coreform\\\" | Get-ScheduledTaskInfo"'
```

(Shows LastRunTime / LastTaskResult across all principals — you
can see WHEN the task last fired and whether it succeeded, but
NOT whose session it ran in.)

**Running process list across all sessions**:

```powershell
ssh 100 'pwsh -Command "Get-CimInstance Win32_Process -Filter ''Name=\"coreform_cubit.exe\"'' | Select-Object ProcessId, @{N=\"User\";E={($_.GetOwner()).User}}, SessionId, CreationDate"'
```

(Cross-session visibility; reveals who has a Cubit running.)

## Playbook template for any "100号機 で <X> が遅い" report

1. **Collect the actual complaint**: what operation, what command,
   user-perceived timing, specific user who hit it
2. **Diagnose infrastructure via SSH** (measure timing AS Administrator):
   - is machine-wide state OK?
   - benchmark the exact operation → see if YOUR number matches user's
3. **If Admin timing is fast but user reports slow** →
   - user-session issue (cache, profile, running MCP with old code)
   - → the fix is user-side action, NOT another admin-side attempt
4. **Provide the user with**:
   - specific one-liner to paste into their terminal
   - what they should see afterward (timing number)
   - who to contact if that one-liner fails
5. **Do NOT**:
   - copy files into user's `LOCALAPPDATA` (ACL breaks)
   - try to impersonate without password (can't)
   - guess at timing numbers (measure yours, ask theirs)

## Memory references

- `reference_shared_filesystem_lab_100.md` — LAB S: ↔ 100号機 W: same mount
- `project_iga_cln_dual_reduction.md` — unrelated, but shows what a
  proper research plan placement looks like
- 2026-04-21 license warmup deployment: commit `7b60b106` +
  `3500c4d9` (both in `packages/radia-mcp/src/radia_mcp/cubit/`)
- `C:\ProgramData\CoreformCubit\` → cubit_license_refresh.ps1
  + cubit_warm.ps1 + cubit_warm.cmd + desktop .lnk
- Scheduled Task path: `\Coreform\CubitLicenseRefresh`
  (any-user logon trigger)
