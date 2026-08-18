# Drive the real EQNEDT64 window with real key presses.
#
# LAST RESORT.  Almost everything this checks is now checked by
# tests/equation/test_key_dispatch.py, which presses the same keys through
# Equation.press() -- no window, no keyboard, no foreground.  That became
# possible once the window stopped deciding anything: it reads the modifier
# state, which only a window can, and hands the press to a function that takes
# the modifiers as arguments.
#
# Taking over somebody's keyboard is what you do when there is no seam.  There
# is one now, so reach for the headless tests first and keep this for the two
# things they genuinely cannot see:
#
#   * that the window is wired to the dispatcher at all
#   * the clipboard round trip, which needs a real clipboard
#
# It refuses to run without -IAmNotUsingThisMachine, because it was being run
# without asking the person sitting at the machine.

param(
    [string]$Exe = "$PSScriptRoot\..\..\build_eq\eqnedt64.exe",
    # It refuses to run without this.  The probe takes the keyboard and the
    # foreground for about a minute, which is intolerable if somebody is at
    # the machine -- and it was being run without asking them.
    [switch]$IAmNotUsingThisMachine
)

if (-not $IAmNotUsingThisMachine) {
    Write-Host "This probe takes over the keyboard and the foreground for about a"
    Write-Host "minute.  Do not run it while you are working at this machine."
    Write-Host ""
    Write-Host "  pwsh -File `"$PSCommandPath`" -IAmNotUsingThisMachine"
    exit 2
}

Add-Type @"
using System;
using System.Runtime.InteropServices;
using System.Text;
public class Fg {
  [DllImport("user32.dll")] public static extern IntPtr GetForegroundWindow();
  [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
  [DllImport("user32.dll")] public static extern uint GetWindowThreadProcessId(IntPtr h, IntPtr pid);
  [DllImport("user32.dll")] public static extern bool AttachThreadInput(uint a, uint b, bool f);
  [DllImport("kernel32.dll")] public static extern uint GetCurrentThreadId();
  [DllImport("user32.dll", CharSet=CharSet.Unicode)] public static extern int GetWindowTextW(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] public static extern bool BringWindowToTop(IntPtr h);
  /* Windows refuses a plain SetForegroundWindow from a background process, so
   * borrow the foreground thread's input state for the length of the call. */
  public static bool Take(IntPtr h) {
    uint me = GetCurrentThreadId();
    uint fg = GetWindowThreadProcessId(GetForegroundWindow(), IntPtr.Zero);
    AttachThreadInput(me, fg, true);
    BringWindowToTop(h); bool ok = SetForegroundWindow(h);
    AttachThreadInput(me, fg, false);
    return ok;
  }
  public static string Title() {
    StringBuilder sb = new StringBuilder(256);
    GetWindowTextW(GetForegroundWindow(), sb, 256);
    return sb.ToString();
  }
}
"@
Add-Type -AssemblyName System.Windows.Forms

if (-not (Test-Path $Exe)) { Write-Error "no editor at $Exe -- build it first"; exit 2 }

$script:failed = 0

function Probe([string]$name, [string[]]$keys, [string]$want) {
    Get-Process eqnedt64 -ErrorAction SilentlyContinue | Stop-Process -Force
    Start-Sleep -Milliseconds 300
    Start-Process $Exe
    Start-Sleep -Milliseconds 1300
    [void][Fg]::Take((Get-Process eqnedt64).MainWindowHandle)
    Start-Sleep -Milliseconds 400

    # Never send keys we cannot prove are going to the editor.
    if ([Fg]::Title() -notlike "EQNEDT64*") {
        "SKIP {0,-24} foreground is '{1}'" -f $name, [Fg]::Title()
        $script:failed++
        return
    }

    Set-Clipboard -Value "(nothing)"
    foreach ($k in $keys) {
        [System.Windows.Forms.SendKeys]::SendWait($k)
        Start-Sleep -Milliseconds 200
    }
    Start-Sleep -Milliseconds 450
    $got = Get-Clipboard

    if ($got.TrimEnd() -eq $want.TrimEnd()) {
        "ok   {0,-24} -> {1}" -f $name, $got
    } else {
        "FAIL {0,-24} -> {1}" -f $name, $got
        "     expected                 -> $want"
        $script:failed++
    }
}

# The templates, each built the way a person builds it: chord, then Tab between
# slots.  The LaTeX is the whole equation, so a stray character shows up.
Probe "Ctrl+F  fraction"    @("^f","dB","{TAB}","dt","^c")           '\frac{dB}{dt}'
Probe "Ctrl+R  root"        @("^r","2","^c")                         '\sqrt{2}'
Probe "Ctrl+H  superscript" @("x","^h","2","^c")                     'x^{2}'
Probe "Ctrl+L  subscript"   @("x","^l","i","^c")                     'x_{i}'
Probe "Ctrl+J  sub+sup"     @("B","^j","z","{TAB}","2","^c")         'B_{z}^{2}'
Probe "Ctrl+I  integral"    @("^i","0","{TAB}","T","{TAB}","f","^c") '\int _{0}^{T}f'
Probe "Ctrl+9  parentheses" @("^9","x","^c")                         '\left( x \right)'

# Two-step chords: the second key must build the template WITHOUT also being
# typed into it.
Probe "Ctrl+T,S summation"  @("^t","s","n","^c")                     '\sum _{n}^{} '
Probe "Ctrl+T,P product"    @("^t","p","k","^c")                     '\prod _{k}^{} '

# Up and down, which no arrow reached before: into a fraction, then out of the
# denominator to the numerator and back.
Probe "Up/Down in a fraction" @("^f","a","{TAB}","b","{UP}","x","^c")  '\frac{ax}{b}'
Probe "Insert is a second Tab" @("^f","a","{INSERT}","b","^c")         '\frac{a}{b}'

# Selection, and the styles that act on it.  Equation Editor names each
# style by its effect: B for bold gives Matrix-Vector, G gives Greek.
Probe "Ctrl+A  select all"  @("a","b","c","^a","^c")                 'abc'
Probe "Ctrl+Shift+B vector" @("B","+{LEFT}","^+b","^c")              '\mathbf{B}'
Probe "Ctrl+Shift+G greek"  @("a","+{LEFT}","^+g","^c")              '\alpha '

# The plain-ctrl editing keys still do their own job.
Probe "Ctrl+X   cut"        @("a","b","+{LEFT}","^x","^c")           'a'
Probe "Ctrl+Z   undo"       @("a","b","^z","^c")                     'a'

Get-Process eqnedt64 -ErrorAction SilentlyContinue | Stop-Process -Force
if ($script:failed) { "`n$script:failed failed"; exit 1 }
"`nall chords reached the editor"
