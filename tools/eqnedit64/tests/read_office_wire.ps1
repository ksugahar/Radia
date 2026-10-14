param([Parameter(Mandatory=$true)][string]$OutputPath)
$ErrorActionPreference='Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1' -or $env:COMPUTERNAME -in @('INTEL11','LAB')) { throw 'Clipboard wire evidence requires a disposable isolated host.' }
Add-Type -TypeDefinition @"
using System; using System.Runtime.InteropServices;
public static class OfficeWire {
 [DllImport("user32.dll")] public static extern bool OpenClipboard(IntPtr hwnd);
 [DllImport("user32.dll")] public static extern bool CloseClipboard();
 [DllImport("user32.dll",CharSet=CharSet.Unicode)] public static extern uint RegisterClipboardFormat(string name);
 [DllImport("user32.dll")] public static extern bool IsClipboardFormatAvailable(uint format);
 [DllImport("user32.dll")] public static extern IntPtr GetClipboardData(uint format);
 [DllImport("kernel32.dll")] public static extern IntPtr GlobalLock(IntPtr value);
 [DllImport("kernel32.dll")] public static extern bool GlobalUnlock(IntPtr value);
 [DllImport("kernel32.dll")] public static extern UIntPtr GlobalSize(IntPtr value);
}
"@
$opened=$false
for($attempt=0;$attempt -lt 40;$attempt++) { if([OfficeWire]::OpenClipboard([IntPtr]::Zero)){$opened=$true;break}; Start-Sleep -Milliseconds 25 }
if(-not $opened){throw 'Cannot open clipboard.'}
try {
 foreach($name in @('MathML','MathML Presentation')) { if([OfficeWire]::IsClipboardFormatAvailable([OfficeWire]::RegisterClipboardFormat($name))){throw "Competing $name published in normal copy."} }
 $handle=[OfficeWire]::GetClipboardData([OfficeWire]::RegisterClipboardFormat('HTML Format'))
 if($handle -eq [IntPtr]::Zero){throw 'Missing HTML Format.'}
 $pointer=[OfficeWire]::GlobalLock($handle)
 if($pointer -eq [IntPtr]::Zero){throw 'Cannot lock HTML Format.'}
 try { $length=[OfficeWire]::GlobalSize($handle).ToUInt64(); if($length -gt 16777216){throw 'Unexpected clipboard size.'}; $bytes=[byte[]]::new([int]$length); [Runtime.InteropServices.Marshal]::Copy($pointer,$bytes,0,$bytes.Length); [IO.File]::WriteAllBytes($OutputPath,$bytes) } finally { $null=[OfficeWire]::GlobalUnlock($handle) }
} finally { $null=[OfficeWire]::CloseClipboard() }
