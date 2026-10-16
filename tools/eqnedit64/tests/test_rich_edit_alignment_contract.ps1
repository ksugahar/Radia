# Exercise the actual RichEdit oracle without native, Office or clipboard use.
$ErrorActionPreference = 'Stop'
$source = Join-Path $PSScriptRoot '..\build\test_external_paste.ps1'
$tokens = $null; $errors = $null
$ast = [Management.Automation.Language.Parser]::ParseFile(
    (Resolve-Path $source).Path, [ref]$tokens, [ref]$errors)
if ($errors) { throw ($errors | Out-String) }
$definition = $ast.Find({ param($node)
    $node -is [Management.Automation.Language.FunctionDefinitionAst] -and
    $node.Name -eq 'Assert-RichEditUnanchoredAlignmentRows'
}, $true)
if (-not $definition) { throw 'The production alignment oracle is missing.' }
. ([scriptblock]::Create($definition.Extent.Text))
$marker = '<m:r><m:rPr><m:aln/></m:rPr><span>&amp;</span></m:r>'
$content = '<m:r><span>a</span></m:r>'
$row = '<m:e>' + $marker + $content + '</m:e>'
$valid = '<m:oMath xmlns:m="http://schemas.microsoft.com/office/2004/12/omml"><m:eqArr>' +
    ($row * 3) + '</m:eqArr></m:oMath>'
Assert-RichEditUnanchoredAlignmentRows ([xml]$valid) 3
$badCases = @{
    missing_marker = $valid.Replace('<m:aln/>', '')
    wrong_marker_text = $valid.Replace('&amp;', 'x')
    visible_ampersand = $valid.Replace($content, '<m:r><span>&amp;</span></m:r>')
    marker_not_first = $valid.Replace($marker + $content, $content + $marker)
    duplicate_marker = $valid.Replace($marker, $marker + $marker)
    missing_row = $valid.Replace($row * 3, $row * 2)
}
foreach ($name in $badCases.Keys) {
    $rejected = $false
    try { Assert-RichEditUnanchoredAlignmentRows ([xml]$badCases[$name]) 3 }
    catch { $rejected = $true }
    if (-not $rejected) { throw "Malformed RichEdit input was accepted: $name" }
}
Write-Output 'PASS: actual RichEdit span/alignment oracle accepts valid rows and rejects six malformed inputs.'
