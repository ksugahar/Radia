<#
  Render one .eqn in Equation Editor 3.0 and capture what IT draws.

  The question this answers: when our reader produces a mess, is the DOCUMENT
  malformed or is the reader wrong?  EE3's own drawing is the arbiter, and it
  is the only one -- a corpus score says a document carries a defect marker,
  never whose fault that is.

  It earned its place on 2026-08-21.  Five corpus documents had been written
  off as "wrong before this work"; EE3 drew the first of them perfectly, which
  turned a closed question back into a bug and led to the fix in
  line_pass.cpp's deepestBareBigOp.

  The window is captured with PrintWindow, so nothing is stolen from whoever
  is at the machine.  EE3 is launched with the file as its argument, via an
  ASCII scratch copy -- it does not open a Japanese path.

  Usage:
    pwsh -File render_in_ee3.ps1 -Eqn <path.eqn> -Out <shot.png>
#>
param(
  [Parameter(Mandatory=$true)][string]$Eqn,
  [Parameter(Mandatory=$true)][string]$Out,
  [string]$Exe = 'C:\Program Files\Microsoft Office\root\vfs\ProgramFilesCommonX64\Microsoft Shared\EQUATION\EQNEDT32.EXE'
)
Add-Type -AssemblyName System.Drawing
Add-Type @'
using System;
using System.Text;
using System.Runtime.InteropServices;
public class R {
  public delegate bool EnumProc(IntPtr h, IntPtr l);
  [DllImport("user32.dll")] public static extern bool EnumWindows(EnumProc cb, IntPtr l);
  [DllImport("user32.dll", CharSet=CharSet.Unicode)] public static extern int GetClassNameW(IntPtr h, StringBuilder s, int n);
  [DllImport("user32.dll")] public static extern uint GetWindowThreadProcessId(IntPtr h, out uint pid);
  [DllImport("user32.dll")] public static extern bool GetClientRect(IntPtr h, out RECT r);
  [DllImport("user32.dll")] public static extern bool PrintWindow(IntPtr h, IntPtr dc, uint flags);
  [DllImport("user32.dll")] public static extern bool SetWindowPos(IntPtr h, IntPtr a, int x, int y, int w, int t, uint f);
  [StructLayout(LayoutKind.Sequential)] public struct RECT { public int L,T,R,B; }
  public static IntPtr Find(uint pid, string cls) {
    IntPtr found = IntPtr.Zero;
    EnumWindows((h,l) => {
      uint p; GetWindowThreadProcessId(h, out p);
      if (p != pid) return true;
      var c = new StringBuilder(256); GetClassNameW(h, c, 256);
      if (c.ToString() == cls) { found = h; return false; }
      return true;
    }, IntPtr.Zero);
    return found;
  }
}
'@

# EE3 opens a path it is given; copy to an ASCII scratch name first.
$tmp = 'C:\temp\_eq_render.eqn'
Copy-Item -LiteralPath $Eqn -Destination $tmp -Force

$p = Start-Process -FilePath $Exe -ArgumentList $tmp -PassThru
Start-Sleep -Seconds 5
$h = [R]::Find([uint32]$p.Id, 'EQNWINCLASS')
if ($h -eq [IntPtr]::Zero) { $p.Kill(); throw 'no EQNWINCLASS window' }
[void][R]::SetWindowPos($h, [IntPtr]::Zero, 100, 100, 900, 360, 0x14)
Start-Sleep -Milliseconds 1200

[R+RECT]$rc = New-Object R+RECT
[void][R]::GetClientRect($h, [ref]$rc)
$bmp = New-Object System.Drawing.Bitmap ($rc.R - $rc.L), ($rc.B - $rc.T)
$g = [System.Drawing.Graphics]::FromImage($bmp)
$dc = $g.GetHdc(); [void][R]::PrintWindow($h, $dc, 1); $g.ReleaseHdc($dc); $g.Dispose()
$bmp.Save($Out, [System.Drawing.Imaging.ImageFormat]::Png); $bmp.Dispose()
$p.Kill()
Write-Output "rendered $(Split-Path $Eqn -Leaf) -> $Out"
