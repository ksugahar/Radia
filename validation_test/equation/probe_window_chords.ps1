# Drive the real EQNEDT64 window with real key presses.
#
# tests/equation/ drives the model directly, which is why 619 of them passed
# while every keyboard shortcut in the editor was dead.  Nothing below the
# Python binding can see a key press: the window is a WIN32 executable, and the
# path from "the user pressed Ctrl+F" to "a fraction appeared" runs through
# WM_KEYDOWN, WM_CHAR and GetKeyState, none of which a unit test touches.
#
# So this presses the keys.  Each case starts a clean editor, sends a chord,
# copies the result to the clipboard and compares the LaTeX.  It found two bugs
# that tests/equation/test_chords.py structurally cannot:
#
#   * the second key of "Ctrl+T, S" was ALSO typed as text, so summation came
#     out as \sum_{sn} -- the chord consumed the key, WM_CHAR delivered it again
#   * Ctrl+Shift+X CUT the selection instead of making it a bold vector, because
#     the built-in Ctrl+X binding never checked whether shift was down
#
# It is not a pytest: it needs a desktop session, it steals the foreground, and
# it types into whatever holds focus if that steal fails -- so it checks the
# window title first and refuses to send anything otherwise.  Run it by hand
# after touching eq_window.cpp or the shortcut table.
#
#   pwsh -File validation_test/equation/probe_window_chords.ps1
#
# Do not run it while you are using the machine; it takes the keyboard.

param(
    [string]$Exe = "$PSScriptRoot\..\..\build_eq\eqnedt64.exe"
)

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
