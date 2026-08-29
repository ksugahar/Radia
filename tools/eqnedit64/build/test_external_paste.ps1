param(
    [string]$GooglePngArtifact,
    [string]$PowerPointPngArtifact,
    [string]$PowerPointPptxArtifact,
    [string]$AppPath
)

$ErrorActionPreference = 'Stop'

if ([Threading.Thread]::CurrentThread.GetApartmentState() -ne
    [Threading.ApartmentState]::STA) {
    throw 'External paste test must run in an STA PowerShell process.'
}

Add-Type -AssemblyName System.Drawing
Add-Type -AssemblyName System.IO.Compression.FileSystem
Add-Type -TypeDefinition @'
using System;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using System.Threading;

public static class EqneditClipboardNative {
    [DllImport("ole32.dll")]
    public static extern int OleInitialize(IntPtr reserved);

    [DllImport("ole32.dll")]
    public static extern void OleUninitialize();

    [DllImport("user32.dll", CharSet = CharSet.Unicode)]
    public static extern uint RegisterClipboardFormat(string format);

    [DllImport("user32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool IsClipboardFormatAvailable(uint format);

    [DllImport("user32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool OpenClipboard(IntPtr owner);

    [DllImport("user32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool CloseClipboard();

    [DllImport("user32.dll")]
    public static extern IntPtr GetClipboardData(uint format);

    [DllImport("user32.dll", SetLastError = true)]
    public static extern uint EnumClipboardFormats(uint format);

    [DllImport("user32.dll", SetLastError = true)]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool EmptyClipboard();

    [DllImport("user32.dll", SetLastError = true)]
    public static extern IntPtr SetClipboardData(uint format, IntPtr memory);

    [DllImport("kernel32.dll")]
    public static extern IntPtr GlobalLock(IntPtr memory);

    [DllImport("kernel32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool GlobalUnlock(IntPtr memory);

    [DllImport("kernel32.dll")]
    public static extern UIntPtr GlobalSize(IntPtr memory);

    [DllImport("kernel32.dll", SetLastError = true)]
    public static extern IntPtr GlobalAlloc(uint flags, UIntPtr bytes);

    [DllImport("kernel32.dll")]
    public static extern IntPtr GlobalFree(IntPtr memory);

    [DllImport("gdi32.dll", CharSet = CharSet.Unicode)]
    public static extern IntPtr CopyEnhMetaFile(IntPtr source, string fileName);

    [DllImport("gdi32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool DeleteEnhMetaFile(IntPtr metafile);

    [DllImport("gdi32.dll", CharSet = CharSet.Unicode)]
    public static extern IntPtr CopyMetaFile(IntPtr source, string fileName);

    [DllImport("gdi32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool DeleteMetaFile(IntPtr metafile);

    [DllImport("user32.dll", SetLastError = true)]
    public static extern IntPtr CopyImage(
        IntPtr image, uint type, int width, int height, uint flags);

    [DllImport("gdi32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    public static extern bool DeleteObject(IntPtr value);

    [StructLayout(LayoutKind.Sequential)]
    private struct BitmapObject {
        public int Type;
        public int Width;
        public int Height;
        public int WidthBytes;
        public ushort Planes;
        public ushort BitsPixel;
        public IntPtr Bits;
    }

    [DllImport("gdi32.dll", EntryPoint = "GetObjectW")]
    private static extern int GetBitmapObject(
        IntPtr value, int bytes, out BitmapObject bitmap);

    [DllImport("user32.dll")]
    private static extern IntPtr GetDC(IntPtr window);

    [DllImport("user32.dll")]
    private static extern int ReleaseDC(IntPtr window, IntPtr dc);

    [DllImport("gdi32.dll")]
    private static extern IntPtr CreateCompatibleDC(IntPtr dc);

    [DllImport("gdi32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool DeleteDC(IntPtr dc);

    [DllImport("gdi32.dll")]
    private static extern IntPtr CreateCompatibleBitmap(
        IntPtr dc, int width, int height);

    [DllImport("gdi32.dll")]
    private static extern IntPtr SelectObject(IntPtr dc, IntPtr value);

    [DllImport("gdi32.dll")]
    [return: MarshalAs(UnmanagedType.Bool)]
    private static extern bool BitBlt(
        IntPtr destination, int xDestination, int yDestination,
        int width, int height, IntPtr source, int xSource, int ySource,
        uint rasterOperation);

    private static IntPtr DuplicateBitmap(IntPtr source) {
        BitmapObject description;
        if (GetBitmapObject(source, Marshal.SizeOf<BitmapObject>(),
                            out description) == 0 ||
            description.Width <= 0 || description.Height <= 0)
            return IntPtr.Zero;
        IntPtr screen = GetDC(IntPtr.Zero);
        IntPtr sourceDc = IntPtr.Zero;
        IntPtr destinationDc = IntPtr.Zero;
        IntPtr copy = IntPtr.Zero;
        IntPtr oldSource = IntPtr.Zero;
        IntPtr oldDestination = IntPtr.Zero;
        try {
            if (screen == IntPtr.Zero) return IntPtr.Zero;
            sourceDc = CreateCompatibleDC(screen);
            destinationDc = CreateCompatibleDC(screen);
            copy = CreateCompatibleBitmap(screen, description.Width,
                                          description.Height);
            if (sourceDc == IntPtr.Zero || destinationDc == IntPtr.Zero ||
                copy == IntPtr.Zero) return IntPtr.Zero;
            oldSource = SelectObject(sourceDc, source);
            oldDestination = SelectObject(destinationDc, copy);
            if (oldSource == IntPtr.Zero || oldDestination == IntPtr.Zero ||
                !BitBlt(destinationDc, 0, 0, description.Width,
                        description.Height, sourceDc, 0, 0, 0x00CC0020))
                return IntPtr.Zero;
            IntPtr result = copy;
            copy = IntPtr.Zero;
            return result;
        } finally {
            if (oldSource != IntPtr.Zero && sourceDc != IntPtr.Zero)
                SelectObject(sourceDc, oldSource);
            if (oldDestination != IntPtr.Zero && destinationDc != IntPtr.Zero)
                SelectObject(destinationDc, oldDestination);
            if (copy != IntPtr.Zero) DeleteObject(copy);
            if (sourceDc != IntPtr.Zero) DeleteDC(sourceDc);
            if (destinationDc != IntPtr.Zero) DeleteDC(destinationDc);
            if (screen != IntPtr.Zero) ReleaseDC(IntPtr.Zero, screen);
        }
    }

    [DllImport("gdi32.dll")]
    public static extern uint GetPaletteEntries(
        IntPtr palette, uint start, uint count, IntPtr entries);

    [DllImport("gdi32.dll")]
    public static extern IntPtr CreatePalette(IntPtr logicalPalette);

    [StructLayout(LayoutKind.Sequential)]
    private struct MetafilePict {
        public int MappingMode;
        public int XExtent;
        public int YExtent;
        public IntPtr Metafile;
    }

    private sealed class ClipboardItem {
        public uint Format;
        public IntPtr Handle;
        public int Kind;
    }

    public sealed class ClipboardSnapshot : IDisposable {
        private readonly List<ClipboardItem> items = new List<ClipboardItem>();
        private bool restored;

        public int FormatCount { get { return items.Count; } }

        private static void OpenWithRetry() {
            for (int attempt = 0; attempt != 40; ++attempt) {
                if (OpenClipboard(IntPtr.Zero)) return;
                Thread.Sleep(25);
            }
            throw new InvalidOperationException("OpenClipboard failed after 1 second.");
        }

        private static ClipboardItem CopyItem(uint format, IntPtr source) {
            if (source == IntPtr.Zero)
                throw new InvalidOperationException("Clipboard format " + format + " has no data handle.");

            if (format == 14 || format == 0x008E) {
                IntPtr copy = CopyEnhMetaFile(source, null);
                if (copy == IntPtr.Zero)
                    throw new InvalidOperationException("CopyEnhMetaFile failed for format " + format + ".");
                return new ClipboardItem { Format = format, Handle = copy, Kind = 1 };
            }
            if (format == 2 || format == 0x0082) {
                IntPtr copy = CopyImage(source, 0, 0, 0, 0x00002000);
                if (copy == IntPtr.Zero)
                    copy = CopyImage(source, 0, 0, 0, 0);
                if (copy == IntPtr.Zero)
                    copy = DuplicateBitmap(source);
                if (copy == IntPtr.Zero)
                    throw new InvalidOperationException("CopyImage failed for format " + format + ".");
                return new ClipboardItem { Format = format, Handle = copy, Kind = 2 };
            }
            if (format == 9) {
                uint count = GetPaletteEntries(source, 0, 0, IntPtr.Zero);
                if (count == 0)
                    throw new InvalidOperationException("GetPaletteEntries failed.");
                int entryBytes = checked((int)count * 4);
                IntPtr entries = Marshal.AllocHGlobal(entryBytes);
                IntPtr logicalPalette = Marshal.AllocHGlobal(entryBytes + 4);
                try {
                    if (GetPaletteEntries(source, 0, count, entries) != count)
                        throw new InvalidOperationException("GetPaletteEntries returned incomplete data.");
                    Marshal.WriteInt16(logicalPalette, 0, unchecked((short)0x0300));
                    Marshal.WriteInt16(logicalPalette, 2, checked((short)count));
                    byte[] bytes = new byte[entryBytes];
                    Marshal.Copy(entries, bytes, 0, entryBytes);
                    Marshal.Copy(bytes, 0, IntPtr.Add(logicalPalette, 4), entryBytes);
                    IntPtr copy = CreatePalette(logicalPalette);
                    if (copy == IntPtr.Zero)
                        throw new InvalidOperationException("CreatePalette failed.");
                    return new ClipboardItem { Format = format, Handle = copy, Kind = 3 };
                } finally {
                    Marshal.FreeHGlobal(entries);
                    Marshal.FreeHGlobal(logicalPalette);
                }
            }
            if (format == 3 || format == 0x0083) {
                IntPtr sourcePointer = GlobalLock(source);
                if (sourcePointer == IntPtr.Zero)
                    throw new InvalidOperationException("GlobalLock failed for METAFILEPICT.");
                MetafilePict value;
                try {
                    value = Marshal.PtrToStructure<MetafilePict>(sourcePointer);
                } finally {
                    GlobalUnlock(source);
                }
                value.Metafile = CopyMetaFile(value.Metafile, null);
                if (value.Metafile == IntPtr.Zero)
                    throw new InvalidOperationException("CopyMetaFile failed.");
                UIntPtr bytes = (UIntPtr)(uint)Marshal.SizeOf<MetafilePict>();
                IntPtr copy = GlobalAlloc(0x0002, bytes);
                if (copy == IntPtr.Zero) {
                    DeleteMetaFile(value.Metafile);
                    throw new InvalidOperationException("GlobalAlloc failed for METAFILEPICT.");
                }
                IntPtr copyPointer = GlobalLock(copy);
                try {
                    Marshal.StructureToPtr(value, copyPointer, false);
                } finally {
                    GlobalUnlock(copy);
                }
                return new ClipboardItem { Format = format, Handle = copy, Kind = 4 };
            }

            ulong byteCount = GlobalSize(source).ToUInt64();
            if (byteCount == 0 || byteCount > int.MaxValue)
                throw new InvalidOperationException(
                    "Clipboard format " + format + " is not a supported HGLOBAL payload.");
            IntPtr copyHandle = GlobalAlloc(0x0002, (UIntPtr)byteCount);
            if (copyHandle == IntPtr.Zero)
                throw new InvalidOperationException("GlobalAlloc failed for format " + format + ".");
            IntPtr sourcePointer2 = GlobalLock(source);
            IntPtr copyPointer2 = GlobalLock(copyHandle);
            try {
                if (sourcePointer2 == IntPtr.Zero || copyPointer2 == IntPtr.Zero)
                    throw new InvalidOperationException("GlobalLock failed for format " + format + ".");
                byte[] bytes = new byte[(int)byteCount];
                Marshal.Copy(sourcePointer2, bytes, 0, bytes.Length);
                Marshal.Copy(bytes, 0, copyPointer2, bytes.Length);
            } catch {
                GlobalFree(copyHandle);
                throw;
            } finally {
                if (copyPointer2 != IntPtr.Zero) GlobalUnlock(copyHandle);
                if (sourcePointer2 != IntPtr.Zero) GlobalUnlock(source);
            }
            return new ClipboardItem { Format = format, Handle = copyHandle, Kind = 0 };
        }

        public static ClipboardSnapshot Capture() {
            Exception lastFailure = null;
            for (int attempt = 0; attempt != 5; ++attempt) {
                var result = new ClipboardSnapshot();
                bool opened = false;
                try {
                    OpenWithRetry();
                    opened = true;
                    uint format = 0;
                    while ((format = EnumClipboardFormats(format)) != 0) {
                        result.items.Add(CopyItem(format, GetClipboardData(format)));
                    }
                    return result;
                } catch (Exception error) {
                    lastFailure = error;
                    result.Dispose();
                } finally {
                    if (opened) CloseClipboard();
                }
                Thread.Sleep(50);
            }
            throw new InvalidOperationException(
                "Clipboard snapshot was unstable for five attempts.", lastFailure);
        }

        public void Restore() {
            OpenWithRetry();
            try {
                if (!EmptyClipboard())
                    throw new InvalidOperationException("EmptyClipboard failed while restoring data.");
                foreach (ClipboardItem item in items) {
                    if (SetClipboardData(item.Format, item.Handle) == IntPtr.Zero)
                        throw new InvalidOperationException(
                            "SetClipboardData failed for format " + item.Format + ".");
                    item.Handle = IntPtr.Zero;
                }
                restored = true;
            } finally {
                CloseClipboard();
            }
        }

        public void Dispose() {
            if (restored) return;
            foreach (ClipboardItem item in items) {
                if (item.Handle == IntPtr.Zero) continue;
                if (item.Kind == 1) DeleteEnhMetaFile(item.Handle);
                else if (item.Kind == 2 || item.Kind == 3) DeleteObject(item.Handle);
                else if (item.Kind == 4) {
                    IntPtr pointer = GlobalLock(item.Handle);
                    if (pointer != IntPtr.Zero) {
                        MetafilePict value = Marshal.PtrToStructure<MetafilePict>(pointer);
                        GlobalUnlock(item.Handle);
                        if (value.Metafile != IntPtr.Zero) DeleteMetaFile(value.Metafile);
                    }
                    GlobalFree(item.Handle);
                } else GlobalFree(item.Handle);
                item.Handle = IntPtr.Zero;
            }
        }
    }
}
'@

function Release-ComObject([object]$Value) {
    if ($null -ne $Value -and [Runtime.InteropServices.Marshal]::IsComObject($Value)) {
        [void][Runtime.InteropServices.Marshal]::FinalReleaseComObject($Value)
    }
}

function Open-ClipboardWithRetry([string]$Context) {
    for ($attempt = 0; $attempt -lt 40; ++$attempt) {
        if ([EqneditClipboardNative]::OpenClipboard([IntPtr]::Zero)) { return }
        Start-Sleep -Milliseconds 25
    }
    throw "OpenClipboard failed while $Context."
}

function Read-ClipboardUtf8([uint32]$Format) {
    Open-ClipboardWithRetry 'reading the raw LaTeX format'
    try {
        $handle = [EqneditClipboardNative]::GetClipboardData($Format)
        if ($handle -eq [IntPtr]::Zero) { return $null }
        $pointer = [EqneditClipboardNative]::GlobalLock($handle)
        if ($pointer -eq [IntPtr]::Zero) { return $null }
        try {
            $size = [int][EqneditClipboardNative]::GlobalSize($handle).ToUInt64()
            $bytes = [byte[]]::new($size)
            [Runtime.InteropServices.Marshal]::Copy($pointer, $bytes, 0, $size)
            $end = [Array]::IndexOf($bytes, [byte]0)
            if ($end -lt 0) { $end = $bytes.Length }
            return [Text.Encoding]::UTF8.GetString($bytes, 0, $end)
        } finally {
            [void][EqneditClipboardNative]::GlobalUnlock($handle)
        }
    } finally {
        [void][EqneditClipboardNative]::CloseClipboard()
    }
}

function Read-ClipboardUtf16([uint32]$Format) {
    Open-ClipboardWithRetry 'reading UTF-16 clipboard data'
    try {
        $handle = [EqneditClipboardNative]::GetClipboardData($Format)
        if ($handle -eq [IntPtr]::Zero) { return $null }
        $pointer = [EqneditClipboardNative]::GlobalLock($handle)
        if ($pointer -eq [IntPtr]::Zero) { return $null }
        try {
            return [Runtime.InteropServices.Marshal]::PtrToStringUni($pointer)
        } finally {
            [void][EqneditClipboardNative]::GlobalUnlock($handle)
        }
    } finally {
        [void][EqneditClipboardNative]::CloseClipboard()
    }
}

function Read-ClipboardBytes([uint32]$Format) {
    Open-ClipboardWithRetry 'reading binary clipboard data'
    try {
        $handle = [EqneditClipboardNative]::GetClipboardData($Format)
        if ($handle -eq [IntPtr]::Zero) { return $null }
        $pointer = [EqneditClipboardNative]::GlobalLock($handle)
        if ($pointer -eq [IntPtr]::Zero) { return $null }
        try {
            $size = [int][EqneditClipboardNative]::GlobalSize($handle).ToUInt64()
            $bytes = [byte[]]::new($size)
            [Runtime.InteropServices.Marshal]::Copy($pointer, $bytes, 0, $size)
            return $bytes
        } finally {
            [void][EqneditClipboardNative]::GlobalUnlock($handle)
        }
    } finally {
        [void][EqneditClipboardNative]::CloseClipboard()
    }
}

function Read-BigEndianUInt32([byte[]]$Bytes, [int]$Offset) {
    return [uint32](
        ([uint32]$Bytes[$Offset] -shl 24) -bor
        ([uint32]$Bytes[$Offset + 1] -shl 16) -bor
        ([uint32]$Bytes[$Offset + 2] -shl 8) -bor
        [uint32]$Bytes[$Offset + 3])
}

function Get-PngContract([byte[]]$Bytes) {
    $signature = [byte[]](137, 80, 78, 71, 13, 10, 26, 10)
    if ($Bytes.Length -lt 33 -or
        [Convert]::ToHexString($Bytes[0..7]) -ne [Convert]::ToHexString($signature)) {
        throw 'Registered PNG clipboard payload has an invalid signature.'
    }
    $width = 0
    $height = 0
    $xPixelsPerMetre = 0
    $yPixelsPerMetre = 0
    $unit = 0
    $offset = 8
    while ($offset + 12 -le $Bytes.Length) {
        $length = [int](Read-BigEndianUInt32 $Bytes $offset)
        if ($length -lt 0 -or $offset + 12 + $length -gt $Bytes.Length) { break }
        $type = [Text.Encoding]::ASCII.GetString($Bytes, $offset + 4, 4)
        if ($type -eq 'IHDR' -and $length -ge 8) {
            $width = Read-BigEndianUInt32 $Bytes ($offset + 8)
            $height = Read-BigEndianUInt32 $Bytes ($offset + 12)
        } elseif ($type -eq 'pHYs' -and $length -eq 9) {
            $xPixelsPerMetre = Read-BigEndianUInt32 $Bytes ($offset + 8)
            $yPixelsPerMetre = Read-BigEndianUInt32 $Bytes ($offset + 12)
            $unit = $Bytes[$offset + 16]
        }
        $offset += 12 + $length
    }
    if ($width -le 0 -or $height -le 0 -or $unit -ne 1 -or
        [Math]::Abs([double]$xPixelsPerMetre - 11811.0) -gt 1.0 -or
        [Math]::Abs([double]$yPixelsPerMetre - 11811.0) -gt 1.0) {
        throw "PNG is not a 300 dpi physical image: ${width}x${height}, ${xPixelsPerMetre}x${yPixelsPerMetre} px/m, unit=$unit"
    }
    return [pscustomobject]@{
        Width = [int]$width
        Height = [int]$height
        XPixelsPerMetre = [int]$xPixelsPerMetre
        YPixelsPerMetre = [int]$yPixelsPerMetre
    }
}

function Read-ClipboardUnicodeText {
    Open-ClipboardWithRetry 'reading CF_UNICODETEXT'
    try {
        $handle = [EqneditClipboardNative]::GetClipboardData(13)
        if ($handle -eq [IntPtr]::Zero) { return $null }
        $pointer = [EqneditClipboardNative]::GlobalLock($handle)
        if ($pointer -eq [IntPtr]::Zero) { return $null }
        try {
            return [Runtime.InteropServices.Marshal]::PtrToStringUni($pointer)
        } finally {
            [void][EqneditClipboardNative]::GlobalUnlock($handle)
        }
    } finally {
        [void][EqneditClipboardNative]::CloseClipboard()
    }
}

function New-ClipboardHGlobal([byte[]]$Bytes) {
    $handle = [EqneditClipboardNative]::GlobalAlloc(
        0x0002, [UIntPtr][uint64]$Bytes.Length)
    if ($handle -eq [IntPtr]::Zero) {
        throw 'GlobalAlloc failed while preparing clipboard input.'
    }
    $pointer = [EqneditClipboardNative]::GlobalLock($handle)
    if ($pointer -eq [IntPtr]::Zero) {
        [void][EqneditClipboardNative]::GlobalFree($handle)
        throw 'GlobalLock failed while preparing clipboard input.'
    }
    try {
        [Runtime.InteropServices.Marshal]::Copy($Bytes, 0, $pointer, $Bytes.Length)
    } finally {
        [void][EqneditClipboardNative]::GlobalUnlock($handle)
    }
    return $handle
}

function Set-ClipboardTexInput(
    [string]$UnicodeText,
    [AllowNull()][string]$RawLatex = $null
) {
    $unicodeBytes = [Text.Encoding]::Unicode.GetBytes($UnicodeText + [char]0)
    $unicodeHandle = New-ClipboardHGlobal $unicodeBytes
    $rawHandle = [IntPtr]::Zero
    if ($null -ne $RawLatex) {
        $rawBytes = [Text.Encoding]::UTF8.GetBytes($RawLatex + [char]0)
        $rawHandle = New-ClipboardHGlobal $rawBytes
    }
    $opened = $false
    try {
        for ($attempt = 0; $attempt -lt 40 -and -not $opened; ++$attempt) {
            $opened = [EqneditClipboardNative]::OpenClipboard([IntPtr]::Zero)
            if (-not $opened) { Start-Sleep -Milliseconds 25 }
        }
        if (-not $opened) { throw 'OpenClipboard failed while setting TeX input.' }
        if (-not [EqneditClipboardNative]::EmptyClipboard()) {
            throw 'EmptyClipboard failed while setting TeX input.'
        }
        if ([EqneditClipboardNative]::SetClipboardData(13, $unicodeHandle) -eq
            [IntPtr]::Zero) {
            throw 'SetClipboardData failed for CF_UNICODETEXT input.'
        }
        $unicodeHandle = [IntPtr]::Zero
        if ($rawHandle -ne [IntPtr]::Zero) {
            $latexFormat = [EqneditClipboardNative]::RegisterClipboardFormat('LaTeX')
            if ($latexFormat -eq 0 -or
                [EqneditClipboardNative]::SetClipboardData($latexFormat, $rawHandle) -eq
                    [IntPtr]::Zero) {
                throw 'SetClipboardData failed for registered LaTeX input.'
            }
            $rawHandle = [IntPtr]::Zero
        }
    } finally {
        if ($opened) { [void][EqneditClipboardNative]::CloseClipboard() }
        if ($unicodeHandle -ne [IntPtr]::Zero) {
            [void][EqneditClipboardNative]::GlobalFree($unicodeHandle)
        }
        if ($rawHandle -ne [IntPtr]::Zero) {
            [void][EqneditClipboardNative]::GlobalFree($rawHandle)
        }
    }
}

function Get-ClipboardFormatCount {
    Open-ClipboardWithRetry 'counting formats'
    try {
        $count = 0
        $format = [uint32]0
        while (($format = [EqneditClipboardNative]::EnumClipboardFormats($format)) -ne 0) {
            ++$count
        }
        return $count
    } finally {
        [void][EqneditClipboardNative]::CloseClipboard()
    }
}

function Get-SlideXml([string]$PptxPath) {
    $archive = [IO.Compression.ZipFile]::OpenRead($PptxPath)
    try {
        $entry = $archive.GetEntry('ppt/slides/slide1.xml')
        if ($null -eq $entry) { throw 'PowerPoint output has no slide1.xml.' }
        $reader = [IO.StreamReader]::new($entry.Open(), [Text.Encoding]::UTF8)
        try { return $reader.ReadToEnd() } finally { $reader.Dispose() }
    } finally {
        $archive.Dispose()
    }
}

function Assert-ImageHasInk([string]$Path, [string]$Label = 'Image') {
    $bitmap = [Drawing.Bitmap]::FromFile($Path)
    try {
        if ($bitmap.Width -le 1 -or $bitmap.Height -le 1) {
            throw "$Label has an invalid image size: $($bitmap.Width)x$($bitmap.Height)"
        }
        $ink = $false
        $xStep = [Math]::Max(1, [int]($bitmap.Width / 256))
        $yStep = [Math]::Max(1, [int]($bitmap.Height / 256))
        for ($y = 0; $y -lt $bitmap.Height -and -not $ink; $y += $yStep) {
            for ($x = 0; $x -lt $bitmap.Width; $x += $xStep) {
                $pixel = $bitmap.GetPixel($x, $y)
                # Transparent PNG backgrounds commonly store RGB=0. They are
                # not visible ink, so alpha must contribute to the pixel too.
                if ($pixel.A -gt 0 -and
                    ($pixel.R -lt 245 -or $pixel.G -lt 245 -or $pixel.B -lt 245)) {
                    $ink = $true
                    break
                }
            }
        }
        if (-not $ink) { throw "$Label is visually blank (no non-white pixels)." }
        return "$($bitmap.Width)x$($bitmap.Height)"
    } finally {
        $bitmap.Dispose()
    }
}

function Assert-DibV5OpaqueBlackOnWhite([byte[]]$Bytes) {
    if ($null -eq $Bytes -or $Bytes.Length -lt 128) {
        throw 'DIBV5 payload is missing or too short.'
    }
    $headerSize = [BitConverter]::ToInt32($Bytes, 0)
    $width = [BitConverter]::ToInt32($Bytes, 4)
    $height = [BitConverter]::ToInt32($Bytes, 8)
    $bitCount = [BitConverter]::ToInt16($Bytes, 14)
    $compression = [BitConverter]::ToInt32($Bytes, 16)
    $imageBytes = [BitConverter]::ToInt32($Bytes, 20)
    if ($headerSize -ne 124 -or $width -le 0 -or $height -ge 0 -or
        $bitCount -ne 32 -or $compression -ne 0 -or $imageBytes -le 0 -or
        $headerSize + $imageBytes -gt $Bytes.Length) {
        throw "Invalid DIBV5 header: size=$headerSize ${width}x${height}, bpp=$bitCount, compression=$compression, bytes=$imageBytes"
    }
    $hasInk = $false
    $hasWhite = $false
    $transparent = 0
    for ($offset = $headerSize; $offset -lt $headerSize + $imageBytes; $offset += 4) {
        if ($Bytes[$offset + 3] -ne 255) { ++$transparent }
        if ($Bytes[$offset] -lt 245 -or $Bytes[$offset + 1] -lt 245 -or
            $Bytes[$offset + 2] -lt 245) {
            $hasInk = $true
        }
        if ($Bytes[$offset] -eq 255 -and $Bytes[$offset + 1] -eq 255 -and
            $Bytes[$offset + 2] -eq 255 -and $Bytes[$offset + 3] -eq 255) {
            $hasWhite = $true
        }
    }
    if ($transparent -ne 0 -or -not $hasInk -or -not $hasWhite) {
        throw "DIBV5 opacity/ink contract failed: transparent=$transparent, ink=$hasInk, white=$hasWhite"
    }
    return "${width}x$(-$height), alpha=255"
}

$projectRoot = Split-Path -Parent $PSScriptRoot
$app = if ($AppPath) {
    (Resolve-Path -LiteralPath $AppPath -ErrorAction Stop).Path
} else {
    Join-Path $projectRoot 'dist\Eqnedit64.exe'
}
$irfanView = 'C:\Program Files\IrfanView\i_view64.exe'
$runId = [Guid]::NewGuid().ToString('N')
$pptxOutput = "C:\temp\Eqnedit64-PowerPoint-paste-$runId.pptx"
$powerPointPngOutput = "C:\temp\Eqnedit64-PowerPoint-paste-$runId.png"
$pngOutput = "C:\temp\Eqnedit64-IrfanView-paste-$runId.png"
$texclipPngOutput = "C:\temp\Eqnedit64-texclip-paste-$runId.png"
$cliTexInput = "C:\temp\Eqnedit64-cli-input-$runId.tex"
# A delimiter space after a control word is valid TeX and proves that the
# public file-based CLI normalises the supplied equation rather than copying a
# private fixed fixture.
$expectedRaw = 'x+\sum_{n=1}^{m} a^3 \int_{a}^{b} \frac{f(x)}{\sqrt{y}}\, dx^3 + \overline{u} + \underline{v}'
$expectedOffice = '\[' + $expectedRaw + '\]'
[IO.File]::WriteAllText(
    $cliTexInput, $expectedOffice, [Text.UTF8Encoding]::new($false))

if (-not (Test-Path -LiteralPath $app)) { throw "Portable app is missing: $app" }
if (-not (Test-Path -LiteralPath $irfanView)) { throw "IrfanView is missing: $irfanView" }
$existingPowerPoint = @(Get-Process -Name POWERPNT -ErrorAction SilentlyContinue)
if ($existingPowerPoint.Count -ne 0) {
    throw 'PowerPoint is already running; refusing to attach to or close a user session.'
}

function Assert-PowerPointEquationRendering([string]$Path) {
    $bitmap = [Drawing.Bitmap]::FromFile($Path)
    try {
        if ($bitmap.Width -le 1 -or $bitmap.Height -le 1) {
            throw "PowerPoint rendered an invalid image size: $($bitmap.Width)x$($bitmap.Height)"
        }

        $inkCount = 0
        $minX = $bitmap.Width
        $minY = $bitmap.Height
        $maxX = -1
        $maxY = -1
        $maxHorizontalRun = 0
        for ($y = 0; $y -lt $bitmap.Height; $y++) {
            $horizontalRun = 0
            for ($x = 0; $x -lt $bitmap.Width; $x++) {
                $pixel = $bitmap.GetPixel($x, $y)
                $isInk = $pixel.A -gt 0 -and
                    ($pixel.R -lt 245 -or $pixel.G -lt 245 -or $pixel.B -lt 245)
                if ($isInk) {
                    $inkCount++
                    $horizontalRun++
                    $minX = [Math]::Min($minX, $x)
                    $minY = [Math]::Min($minY, $y)
                    $maxX = [Math]::Max($maxX, $x)
                    $maxY = [Math]::Max($maxY, $y)
                    $maxHorizontalRun = [Math]::Max($maxHorizontalRun, $horizontalRun)
                } else {
                    $horizontalRun = 0
                }
            }
        }
        if ($inkCount -eq 0) {
            throw 'PowerPoint-rendered pasted equation is visually blank.'
        }

        if ($minX -gt 24) {
            throw ("PowerPoint-rendered equation is not left-aligned: " +
                "first ink is at x=$minX px.")
        }

        $inkWidth = $maxX - $minX + 1
        $inkHeight = $maxY - $minY + 1
        # The fixed acceptance equation contains a fraction.  A replacement
        # glyph or tofu box has ink but cannot contain the long fraction rule.
        # Requiring its silhouette prevents the former "any black pixel" false
        # positive while remaining independent of font antialiasing.
        # A fixed 40 px rule detects the known fraction in this 24 pt fixture
        # without making the threshold depend on unrelated terms added to the
        # left-alignment sentinel prefix.
        $requiredFractionRun = 40
        if ($inkWidth -lt 40 -or $inkHeight -lt 20 -or
            $maxHorizontalRun -lt $requiredFractionRun) {
            throw ("PowerPoint did not render the expected fraction silhouette: " +
                "ink=${inkWidth}x${inkHeight}, longest horizontal run=" +
                "$maxHorizontalRun px; required width>=40 px, height>=20 px, " +
                "horizontal run>=$requiredFractionRun px.")
        }
        return "$($bitmap.Width)x$($bitmap.Height); ink-left=$minX; ink=${inkWidth}x${inkHeight}; fraction-run=$maxHorizontalRun"
    } finally {
        $bitmap.Dispose()
    }
}

function Write-ImageOnWhite([string]$InputPath, [string]$OutputPath) {
    $source = [Drawing.Bitmap]::FromFile($InputPath)
    $composite = [Drawing.Bitmap]::new(
        $source.Width, $source.Height,
        [Drawing.Imaging.PixelFormat]::Format32bppArgb)
    $graphics = [Drawing.Graphics]::FromImage($composite)
    try {
        $graphics.Clear([Drawing.Color]::White)
        $graphics.DrawImage($source, 0, 0, $source.Width, $source.Height)
        $composite.Save($OutputPath, [Drawing.Imaging.ImageFormat]::Png)
    } finally {
        $graphics.Dispose()
        $composite.Dispose()
        $source.Dispose()
    }
}

$oleResult = [EqneditClipboardNative]::OleInitialize([IntPtr]::Zero)
if ($oleResult -lt 0) { throw ('OleInitialize failed: 0x{0:X8}' -f $oleResult) }
$clipboardSnapshot = [EqneditClipboardNative+ClipboardSnapshot]::Capture()
$originalFormatCount = $clipboardSnapshot.FormatCount
$originalHadUnicodeText = [EqneditClipboardNative]::IsClipboardFormatAvailable(13)
$originalUnicodeText = if ($originalHadUnicodeText) { Read-ClipboardUnicodeText } else { $null }
$powerPoint = $null
$presentation = $null
$slide = $null
$shapeRange = $null
$restoreFailure = $null
$completed = $false

try {
    $publisher = Start-Process -FilePath $app `
        -ArgumentList @('--copy-tex-file', $cliTexInput) `
        -WorkingDirectory (Split-Path -Parent $app) -WindowStyle Hidden -Wait -PassThru
    if ($publisher.ExitCode -ne 0) {
        throw "Eqnedit64 clipboard publication failed: $($publisher.ExitCode)"
    }

    $latexFormat = [EqneditClipboardNative]::RegisterClipboardFormat('LaTeX')
    $mathMlFormat = [EqneditClipboardNative]::RegisterClipboardFormat('MathML')
    $mathMlPresentationFormat = [EqneditClipboardNative]::RegisterClipboardFormat(
        'MathML Presentation')
    $htmlFormat = [EqneditClipboardNative]::RegisterClipboardFormat('HTML Format')
    $requiredFormats = @(
        13, 14, 17, $latexFormat, $htmlFormat)
    foreach ($requiredFormat in $requiredFormats) {
        if (-not [EqneditClipboardNative]::IsClipboardFormatAvailable($requiredFormat)) {
            throw "Required clipboard format is missing: $requiredFormat"
        }
    }
    if ((Read-ClipboardUnicodeText) -ne $expectedOffice) {
        throw 'CF_UNICODETEXT does not contain the Office-recognisable TeX wrapper.'
    }
    $actualRaw = Read-ClipboardUtf8 $latexFormat
    if ($actualRaw -ne $expectedRaw) {
        throw ("Registered LaTeX clipboard format does not contain the raw " +
            "TeX fragment: actual=<$actualRaw> expected=<$expectedRaw>")
    }
    $officeHtml = Read-ClipboardUtf8 $htmlFormat
    $startMarker = '<!--StartFragment-->'
    $endMarker = '<!--EndFragment-->'
    $fragmentStart = $officeHtml.IndexOf($startMarker)
    $fragmentEnd = $officeHtml.IndexOf($endMarker)
    if ($fragmentStart -lt 0 -or $fragmentEnd -le $fragmentStart) {
        throw 'CF_HTML fragment markers are missing.'
    }
    $fragmentStart += $startMarker.Length
    $fragment = $officeHtml.Substring(
        $fragmentStart, $fragmentEnd - $fragmentStart)
    $mathMl = if ($fragment.EndsWith('&#160;')) {
        $fragment.Substring(0, $fragment.Length - '&#160;'.Length)
    } else { '' }
    if ([EqneditClipboardNative]::IsClipboardFormatAvailable($mathMlFormat) -or
        [EqneditClipboardNative]::IsClipboardFormatAvailable($mathMlPresentationFormat)) {
        throw 'Normal copy exposed registered MathML that PowerPoint centres.'
    }
    if ($officeHtml -notmatch '<math\b' -or
        $officeHtml -notmatch '</math>&#160;<!--EndFragment-->' -or
        $mathMl -notmatch 'display="inline"' -or
        $mathMl -notmatch 'mathsize="24pt"' -or
        $mathMl -notmatch '<mfrac>' -or $mathMl -notmatch '<msqrt>' -or
        $mathMl -notmatch '<munderover><mo[^>]*>&#x2211;</mo>' -or
        $mathMl -notmatch '<msubsup><mo[^>]*>&#x222B;</mo>') {
        throw 'Eqnedit64 did not publish inline 24 pt structural MathML in CF_HTML.'
    }
    $dibContract = Assert-DibV5OpaqueBlackOnWhite (Read-ClipboardBytes 17)

    $powerPoint = New-Object -ComObject PowerPoint.Application
    $powerPoint.DisplayAlerts = 0
    $presentation = $powerPoint.Presentations.Add(0)
    $slide = $presentation.Slides.Add(1, 12)
    $shapeRange = $slide.Shapes.Paste()
    if ($shapeRange.Count -ne 1) {
        throw "PowerPoint pasted an unexpected shape count: $($shapeRange.Count)"
    }
    $pastedShape = $shapeRange.Item(1)
    $powerPointFontSize = [double]$pastedShape.TextFrame2.TextRange.Font.Size
    $powerPointAlignment =
        [int]$pastedShape.TextFrame2.TextRange.ParagraphFormat.Alignment
    $powerPointLeft = [double]$pastedShape.Left
    if ([Math]::Abs($powerPointFontSize - 24.0) -gt 0.1) {
        throw "PowerPoint native equation is not 24 pt: $powerPointFontSize pt."
    }
    if ($powerPointAlignment -ne 1 -or $powerPointLeft -gt 1.0) {
        throw ("PowerPoint native equation is not left-aligned: " +
            "paragraphAlignment=$powerPointAlignment, left=$powerPointLeft pt.")
    }
    $presentation.SaveAs($pptxOutput, 24)
    $slideXml = Get-SlideXml $pptxOutput
    $hasInlineContainer = $slideXml -match '<a14:m(?:\s|>)'
    $hasInlineMath = $slideXml -match '<m:oMath(?:\s|>)'
    $hasDisplayMath = $slideXml -match '<m:oMathPara(?:\s|>)'
    $hasFraction = $slideXml -match '<m:f>'
    $hasRadical = $slideXml -match '<m:rad>'
    $naryCount = ([regex]::Matches($slideXml, '<m:nary>')).Count
    $barCount = ([regex]::Matches($slideXml, '<m:bar>')).Count
    $accentCount = ([regex]::Matches($slideXml, '<m:acc>')).Count
    if (-not $hasInlineContainer -or -not $hasInlineMath -or $hasDisplayMath -or
        -not $hasFraction -or -not $hasRadical -or $naryCount -lt 2 -or
        ($barCount + $accentCount) -lt 2) {
        throw (("PowerPoint paste contract failed: inlineContainer={0}, " +
            "inlineMath={1}, displayMath={2}, fraction={3}, radical={4}, nary={5}, " +
            "bars={6}, accents={7}, artifact={8}") -f $hasInlineContainer, $hasInlineMath,
            $hasDisplayMath, $hasFraction, $hasRadical, $naryCount, $barCount,
            $accentCount, $pptxOutput)
    }
    # XML and MathZones can both exist while PowerPoint draws no equation.
    # Shape.Export invokes PowerPoint's own renderer without a visible window;
    # the silhouette check below rejects blank output and replacement glyphs.
    $shapeRange.Item(1).Export($powerPointPngOutput, 2) # ppShapeFormatPNG
    if (-not (Test-Path -LiteralPath $powerPointPngOutput)) {
        throw 'PowerPoint did not export the pasted slide image.'
    }
    $powerPointImageSize = Assert-PowerPointEquationRendering $powerPointPngOutput

    $irfan = Start-Process -FilePath $irfanView `
        -ArgumentList @('/clippaste', "/convert=$pngOutput", '/silent') `
        -WindowStyle Hidden -Wait -PassThru
    if ($irfan.ExitCode -ne 0 -or -not (Test-Path -LiteralPath $pngOutput)) {
        throw "IrfanView clipboard paste failed: $($irfan.ExitCode)"
    }
    $imageSize = Assert-ImageHasInk $pngOutput 'IrfanView clipboard paste'

    $googlePublisher = Start-Process -FilePath $app `
        -ArgumentList @('--copy-google-slides-file', $cliTexInput) `
        -WorkingDirectory (Split-Path -Parent $app) -WindowStyle Hidden -Wait -PassThru
    if ($googlePublisher.ExitCode -ne 0) {
        throw "Eqnedit64 Google Slides clipboard publication failed: $($googlePublisher.ExitCode)"
    }
    $pngFormat = [EqneditClipboardNative]::RegisterClipboardFormat('PNG')
    $htmlFormat = [EqneditClipboardNative]::RegisterClipboardFormat('HTML Format')
    foreach ($requiredFormat in @($pngFormat, $htmlFormat)) {
        if (-not [EqneditClipboardNative]::IsClipboardFormatAvailable($requiredFormat)) {
            throw "Google Slides clipboard format is missing: $requiredFormat"
        }
    }
    $googlePng = Read-ClipboardBytes $pngFormat
    $pngContract = Get-PngContract $googlePng
    $googleHtml = Read-ClipboardUtf8 $htmlFormat
    $sizeMatch = [regex]::Match(
        $googleHtml, 'style="width:([0-9.]+)pt;height:([0-9.]+)pt"')
    $imageMatch = [regex]::Match(
        $googleHtml, 'src="data:image/png;base64,([^"]+)"')
    if (-not $sizeMatch.Success -or -not $imageMatch.Success -or
        $googleHtml -notmatch 'data-eqnedit-dpi="300"' -or
        $googleHtml -notmatch 'data-eqnedit-font-size="24pt"') {
        throw 'Google Slides HTML does not declare the 300 dpi / 24 pt image contract.'
    }
    $widthPoints = [double]::Parse(
        $sizeMatch.Groups[1].Value, [Globalization.CultureInfo]::InvariantCulture)
    $heightPoints = [double]::Parse(
        $sizeMatch.Groups[2].Value, [Globalization.CultureInfo]::InvariantCulture)
    if ([Math]::Abs($pngContract.Width * 72.0 / 300.0 - $widthPoints) -gt 0.001 -or
        [Math]::Abs($pngContract.Height * 72.0 / 300.0 - $heightPoints) -gt 0.001) {
        throw 'Google Slides PNG pixels and HTML point size describe different physical sizes.'
    }
    $embeddedPng = [Convert]::FromBase64String($imageMatch.Groups[1].Value)
    $clipboardHash = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($googlePng))
    $embeddedHash = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($embeddedPng))
    if ($clipboardHash -ne $embeddedHash) {
        throw 'Registered PNG and HTML-embedded PNG are not byte-identical.'
    }

    Set-ClipboardTexInput -UnicodeText $expectedOffice
    $texclip = Start-Process -FilePath $app -ArgumentList '--texclip' `
        -WorkingDirectory (Split-Path -Parent $app) -WindowStyle Hidden -Wait -PassThru
    if ($texclip.ExitCode -ne 0) {
        throw "Eqnedit64 --texclip Unicode conversion failed: $($texclip.ExitCode)"
    }
    if (-not [EqneditClipboardNative]::IsClipboardFormatAvailable($pngFormat) -or
        -not [EqneditClipboardNative]::IsClipboardFormatAvailable(17) -or
        [EqneditClipboardNative]::IsClipboardFormatAvailable($htmlFormat) -or
        [EqneditClipboardNative]::IsClipboardFormatAvailable(13)) {
        throw '--texclip did not replace text with PNG/DIBV5 image formats.'
    }
    $texclipUnicodePng = Read-ClipboardBytes $pngFormat
    $texclipContract = Get-PngContract $texclipUnicodePng
    $texclipUnicodeHash = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($texclipUnicodePng))
    if ($texclipUnicodeHash -ne $clipboardHash) {
        throw '--texclip Unicode wrapper normalization rendered a different equation.'
    }

    Set-ClipboardTexInput -UnicodeText '\[x^{1000}\]' -RawLatex $expectedRaw
    $texclip = Start-Process -FilePath $app `
        -ArgumentList '--clipboard-tex-to-png' `
        -WorkingDirectory (Split-Path -Parent $app) -WindowStyle Hidden -Wait -PassThru
    if ($texclip.ExitCode -ne 0) {
        throw "Eqnedit64 --texclip registered LaTeX conversion failed: $($texclip.ExitCode)"
    }
    $texclipRawPng = Read-ClipboardBytes $pngFormat
    $texclipRawHash = [Convert]::ToHexString(
        [Security.Cryptography.SHA256]::HashData($texclipRawPng))
    if ($texclipRawHash -ne $clipboardHash) {
        throw '--texclip did not prefer the registered raw LaTeX clipboard format.'
    }
    $texclipDib = Read-ClipboardBytes 17
    if ($texclipDib.Length -lt 124 -or
        [BitConverter]::ToInt32($texclipDib, 4) -ne $texclipContract.Width -or
        [BitConverter]::ToInt32($texclipDib, 8) -ne -$texclipContract.Height -or
        [Math]::Abs([BitConverter]::ToInt32($texclipDib, 24) - 11811) -gt 1 -or
        [Math]::Abs([BitConverter]::ToInt32($texclipDib, 28) - 11811) -gt 1) {
        throw '--texclip DIBV5 fallback does not match the PNG size and 300 dpi.'
    }
    $texclipDibContract = Assert-DibV5OpaqueBlackOnWhite $texclipDib
    $texclipIrfan = Start-Process -FilePath $irfanView `
        -ArgumentList @('/clippaste', "/convert=$texclipPngOutput", '/silent') `
        -WindowStyle Hidden -Wait -PassThru
    if ($texclipIrfan.ExitCode -ne 0 -or
        -not (Test-Path -LiteralPath $texclipPngOutput)) {
        throw "IrfanView could not paste --texclip PNG: $($texclipIrfan.ExitCode)"
    }
    $texclipImageSize = Assert-ImageHasInk $texclipPngOutput 'IrfanView texclip paste'

    Set-ClipboardTexInput -UnicodeText ''
    $emptyTexclip = Start-Process -FilePath $app -ArgumentList '--texclip' `
        -WorkingDirectory (Split-Path -Parent $app) -WindowStyle Hidden -Wait -PassThru
    if ($emptyTexclip.ExitCode -ne 83 -or
        -not [EqneditClipboardNative]::IsClipboardFormatAvailable(13) -or
        [EqneditClipboardNative]::IsClipboardFormatAvailable($pngFormat)) {
        throw '--texclip did not reject empty input without replacing the clipboard.'
    }

    if ($GooglePngArtifact) {
        [IO.File]::WriteAllBytes($GooglePngArtifact, $googlePng)
    }
    if ($PowerPointPngArtifact) {
        Write-ImageOnWhite $powerPointPngOutput $PowerPointPngArtifact
    }
    if ($PowerPointPptxArtifact) {
        Copy-Item -LiteralPath $pptxOutput -Destination $PowerPointPptxArtifact
    }

    Write-Host "PASS: normal DIBV5 is opaque black-on-white ($dibContract)"
    Write-Host 'PASS: clipboard contains raw LaTeX, inline MathML CF_HTML, Office TeX, EMF, and DIBV5 without centring registered MathML'
    Write-Host ("PASS: no-selection GUI copy -> visible editable Office Math " +
        "in PowerPoint ($powerPointFontSize pt, left=$powerPointLeft pt, " +
        "rendered $powerPointImageSize with ink)")
    Write-Host "PASS: IrfanView /clippaste produced a nonblank $imageSize PNG"
    Write-Host ("PASS: Google Slides clipboard contains a byte-identical " +
        "$($pngContract.Width)x$($pngContract.Height) 300 dpi PNG and " +
        "${widthPoints}x${heightPoints} pt HTML (24 pt base)")
    Write-Host ("PASS: --texclip replaced Unicode/LaTeX input with a nonblank " +
        "$texclipImageSize PNG/DIBV5 image at $($texclipContract.XPixelsPerMetre) px/m; $texclipDibContract")
    $completed = $true
} finally {
    Release-ComObject $shapeRange
    Release-ComObject $slide
    if ($presentation) { try { $presentation.Close() } catch {} }
    Release-ComObject $presentation
    if ($powerPoint) { try { $powerPoint.Quit() } catch {} }
    Release-ComObject $powerPoint
    try {
        $clipboardSnapshot.Restore()
        if ((Get-ClipboardFormatCount) -ne $originalFormatCount) {
            throw 'Clipboard restoration changed the number of materialised formats.'
        }
        if ([EqneditClipboardNative]::IsClipboardFormatAvailable(13) -ne $originalHadUnicodeText) {
            throw 'Clipboard restoration changed CF_UNICODETEXT availability.'
        }
        if ($originalHadUnicodeText -and (Read-ClipboardUnicodeText) -ne $originalUnicodeText) {
            throw 'Clipboard restoration changed the original Unicode text.'
        }
    } catch {
        $restoreFailure = $_.Exception
    }
    $clipboardSnapshot.Dispose()
    [EqneditClipboardNative]::OleUninitialize()
    [GC]::Collect()
    [GC]::WaitForPendingFinalizers()
    $cleanupPaths = @(
        $powerPointPngOutput, $pngOutput, $texclipPngOutput, $cliTexInput)
    if ($completed) { $cleanupPaths += $pptxOutput }
    foreach ($path in $cleanupPaths) {
        if (Test-Path -LiteralPath $path) { Remove-Item -LiteralPath $path -Force }
    }
    if ($restoreFailure) { throw $restoreFailure }
}
