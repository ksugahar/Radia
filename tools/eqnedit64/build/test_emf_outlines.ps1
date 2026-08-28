# The clipboard metafile must not depend on a font.
#
# The EMF is what PowerPoint receives when an equation is pasted as a picture.
# It used to record text plus a font name, and Latin Modern Math is loaded
# into Eqnedit64's own process alone -- so PowerPoint substituted and the
# equation stopped looking like TeX on the machine that made it.  Glyphs are
# recorded as outlines now.
#
# The check reads the metafile's records rather than looking at the result.
# Comparing rendered shapes was tried first and is not decisive: substitution
# changes the letterforms but barely moves the overall bounding box, so an
# aspect-ratio comparison passed happily with text records still in the file.
# "Contains no text record" is exact, and it is the property that matters.

$ErrorActionPreference = 'Stop'

$root = Split-Path -Parent $PSScriptRoot
$exe = Join-Path $root 'build\Eqnedit64.exe'
if (-not (Test-Path $exe)) { $exe = Join-Path $root 'dist\Eqnedit64.exe' }
if (-not (Test-Path $exe)) {
    Write-Output 'skip  Eqnedit64.exe has not been built'
    exit 0
}

# EMR record types that draw text, and therefore need the named font present.
$EMR_EXTTEXTOUTA = 83
$EMR_EXTTEXTOUTW = 84
$EMR_SMALLTEXTOUT = 108
$EMR_POLYPOLYGON = 8
$EMR_POLYPOLYGON16 = 91

$work = Join-Path ([System.IO.Path]::GetTempPath()) ('eqnemf-' + [guid]::NewGuid().ToString('N'))
New-Item -ItemType Directory -Path $work | Out-Null
$failures = @()

foreach ($tex in @('\int_a^b f(x)dx', '\frac{a}{b}', '\sqrt{x}', 'x^{2}+\alpha')) {
    $emf = Join-Path $work 'clip.emf'
    if (Test-Path $emf) { Remove-Item -LiteralPath $emf }
    & $exe --render-emf $tex $emf | Out-Null
    if (-not (Test-Path $emf)) { $failures += "$tex : no metafile was written"; continue }

    $bytes = [System.IO.File]::ReadAllBytes($emf)
    $at = 0
    $textRecords = 0
    $polygonRecords = 0
    while ($at + 8 -le $bytes.Length) {
        $type = [BitConverter]::ToUInt32($bytes, $at)
        $size = [BitConverter]::ToUInt32($bytes, $at + 4)
        if ($size -lt 8) { break }
        if ($type -eq $EMR_EXTTEXTOUTA -or $type -eq $EMR_EXTTEXTOUTW -or
            $type -eq $EMR_SMALLTEXTOUT) { $textRecords++ }
        if ($type -eq $EMR_POLYPOLYGON -or $type -eq $EMR_POLYPOLYGON16) {
            $polygonRecords++
        }
        $at += [int]$size
    }

    if ($textRecords -gt 0) {
        $failures += ("{0} : the metafile holds {1} text record(s) -- whatever opens it will substitute the font" -f $tex, $textRecords)
    }
    if ($polygonRecords -eq 0) {
        $failures += ("{0} : the metafile holds no filled outlines, so it draws nothing" -f $tex)
    }
}

[System.IO.Directory]::Delete($work, $true)

if ($failures.Count) {
    Write-Output ("FAIL  " + $failures.Count)
    $failures | ForEach-Object { Write-Output ("  " + $_) }
    exit 1
}
Write-Output 'ok    the clipboard metafile carries outlines and names no font'
exit 0
