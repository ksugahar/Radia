[CmdletBinding()]
param([Parameter(Mandatory)] [string]$DumpDirectory)
$ErrorActionPreference = 'Stop'
if ($env:EQNEDIT64_ISOLATED_TEST_SESSION -ne '1' -or $env:GITHUB_ACTIONS -ne 'true') {
    throw 'This diagnostic changes only disposable GitHub CI VMs, never LAB/100.'
}
$session = (Get-Process -Id $PID).SessionId
$hosts = @(Get-CimInstance Win32_Process -Filter "Name='fontdrvhost.exe'" | Where-Object SessionId -EQ $session)
if (-not $hosts.Count) { throw 'No font-host in the diagnostic session.' }
$sids = @(foreach ($hostProcess in $hosts) {
    $owner = Invoke-CimMethod -InputObject $hostProcess -MethodName GetOwnerSid
    if ($owner.ReturnValue -ne 0 -or $owner.Sid -notlike 'S-1-5-96-*') {
        throw 'Cannot establish the session UMFD account; do not broaden the ACL.'
    }
    $owner.Sid
}) | Select-Object -Unique
$directory = [IO.Path]::GetFullPath($DumpDirectory)
if (Test-Path -LiteralPath $directory) { throw 'Use a fresh restricted dump directory.' }
New-Item -ItemType Directory -Path $directory | Out-Null
$acl = Get-Acl -LiteralPath $directory
$acl.SetAccessRuleProtection($true, $false)
foreach ($sid in @('S-1-5-18', 'S-1-5-32-544') + $sids) {
    $identity = [Security.Principal.SecurityIdentifier]::new($sid)
    $rule = [Security.AccessControl.FileSystemAccessRule]::new(
        $identity, 'FullControl', 'ContainerInherit,ObjectInherit', 'None', 'Allow')
    $acl.AddAccessRule($rule)
}
Set-Acl -LiteralPath $directory -AclObject $acl
$key = 'HKLM:\SOFTWARE\Microsoft\Windows\Windows Error Reporting\LocalDumps\fontdrvhost.exe'
New-Item -Path $key -Force | Out-Null
New-ItemProperty -Path $key -Name DumpFolder -PropertyType ExpandString -Value $directory -Force | Out-Null
New-ItemProperty -Path $key -Name DumpType -PropertyType DWord -Value 2 -Force | Out-Null
New-ItemProperty -Path $key -Name DumpCount -PropertyType DWord -Value 1 -Force | Out-Null
[ordered]@{session=$session; umfd_sids=@($sids); dump_folder=$directory; dump_type=2; dump_count=1
    registry_key=$key; acl_sddl=(Get-Acl -LiteralPath $directory).Sddl
} | ConvertTo-Json -Depth 4 | Set-Content (Join-Path $directory 'configuration.json') -Encoding utf8
