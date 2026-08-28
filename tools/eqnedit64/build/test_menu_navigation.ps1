# Opening and walking the menus must not freeze the program.
#
# The responsiveness test presses every command through WM_COMMAND, which
# proved they all RETURN -- but the user's report is that touching some menus
# freezes the editor, and touching happens before any command runs: opening
# the popup, drawing its items, moving the highlight.  None of that is
# exercised by dispatching the command directly.
#
# So this drives the real menu loop: it starts the editor visibly, enters the
# menu bar with the keyboard, walks every top-level menu and steps through all
# of its items and submenus, and between steps asks the window whether it is
# still answering messages.  A freeze names the menu being walked at the time.

$ErrorActionPreference = 'Stop'

Add-Type @'
using System;
using System.Runtime.InteropServices;
public class MenuNav {
  [DllImport("user32.dll")] public static extern bool PostMessageW(IntPtr h, uint m, IntPtr w, IntPtr l);
  [DllImport("user32.dll", CharSet=CharSet.Auto)]
  public static extern IntPtr SendMessageTimeout(IntPtr h, uint m, IntPtr w, IntPtr l,
                                                 uint flags, uint timeoutMs, out IntPtr result);
  [DllImport("user32.dll")] public static extern bool SetForegroundWindow(IntPtr h);
  public const uint WM_NULL = 0x0000;
  public const uint WM_KEYDOWN = 0x0100;
  public const uint WM_KEYUP = 0x0101;
  public const uint WM_SYSCOMMAND = 0x0112;
  public const uint SC_KEYMENU = 0xF100;
  public static bool Alive(IntPtr h, uint timeoutMs) {
    IntPtr r;
    return SendMessageTimeout(h, WM_NULL, IntPtr.Zero, IntPtr.Zero, 2, timeoutMs, out r) != IntPtr.Zero;
  }
  public static void Key(IntPtr h, int vk) {
    PostMessageW(h, WM_KEYDOWN, (IntPtr)vk, IntPtr.Zero);
    PostMessageW(h, WM_KEYUP, (IntPtr)vk, IntPtr.Zero);
  }
}
'@

$root = Split-Path -Parent $PSScriptRoot
$exe = Join-Path $root 'build\Eqnedit64.exe'
if (-not (Test-Path $exe)) { $exe = Join-Path $root 'dist\Eqnedit64.exe' }
if (-not (Test-Path $exe)) {
    Write-Output 'skip  Eqnedit64.exe has not been built'
    exit 0
}

$VK_ESC = 0x1B; $VK_LEFT = 0x25; $VK_UP = 0x26; $VK_RIGHT = 0x27; $VK_DOWN = 0x28

# Seven top-level menus; the deepest popup holds 30 cells.  Walking down past
# the end wraps, so a fixed step count covers every item without needing the
# exact counts.
$topLevelMenus = 7
$stepsPerMenu = 40

$proc = Start-Process -FilePath $exe -PassThru
try {
    $deadline = [DateTime]::Now.AddSeconds(15)
    while ([DateTime]::Now -lt $deadline -and $proc.MainWindowHandle -eq 0) {
        Start-Sleep -Milliseconds 200
        $proc.Refresh()
    }
    $hwnd = $proc.MainWindowHandle
    if ($hwnd -eq 0) { Write-Output 'FAIL  the editor never opened a window'; exit 1 }
    [void][MenuNav]::SetForegroundWindow($hwnd)
    Start-Sleep -Milliseconds 500

    $frozeAt = $null
    for ($menu = 0; $menu -lt $topLevelMenus -and -not $frozeAt; $menu++) {
        # Enter the menu bar, move right to this menu, open it.
        [void][MenuNav]::PostMessageW($hwnd, [MenuNav]::WM_SYSCOMMAND,
                                      [IntPtr][MenuNav]::SC_KEYMENU, [IntPtr]0)
        Start-Sleep -Milliseconds 60
        for ($i = 0; $i -lt $menu; $i++) { [MenuNav]::Key($hwnd, $VK_RIGHT); Start-Sleep -Milliseconds 25 }
        [MenuNav]::Key($hwnd, $VK_DOWN)          # open the popup
        Start-Sleep -Milliseconds 60

        for ($step = 0; $step -lt $stepsPerMenu; $step++) {
            [MenuNav]::Key($hwnd, $VK_DOWN)
            if ($step % 5 -eq 2) { [MenuNav]::Key($hwnd, $VK_RIGHT) }  # dip into submenus
            if ($step % 5 -eq 4) { [MenuNav]::Key($hwnd, $VK_LEFT) }
            Start-Sleep -Milliseconds 15
            if (-not [MenuNav]::Alive($hwnd, 3000)) {
                $frozeAt = "menu #$menu, step $step"
                break
            }
        }
        [MenuNav]::Key($hwnd, $VK_ESC); Start-Sleep -Milliseconds 30
        [MenuNav]::Key($hwnd, $VK_ESC); Start-Sleep -Milliseconds 30
        if (-not [MenuNav]::Alive($hwnd, 3000)) { $frozeAt = "closing menu #$menu"; break }
    }

    if ($frozeAt) {
        Write-Output ("FAIL  the editor stopped answering while walking " + $frozeAt)
        exit 1
    }
    Write-Output ("ok    walked {0} menus, {1} steps each; the editor kept answering" -f `
                  $topLevelMenus, $stepsPerMenu)
    exit 0
} finally {
    if (-not $proc.HasExited) { try { $proc.Kill() } catch { } }
}
