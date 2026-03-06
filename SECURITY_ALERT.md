# Security Alert: PyPI API Token Exposure

**Date**: 2025-11-06
**Severity**: High
**Status**: Mitigated (Script Fixed)

## Issue

GitGuardian detected a Python Package Index (PyPI) API token exposed in the repository.

## Root Cause

The `Publish_to_PyPI.ps1` script contained a hardcoded PyPI API token on line 137:

```powershell
python -m twine upload dist/*.whl dist/*.tar.gz --username __token__ --password pypi-AgEI...
```

**Note**: This file is listed in `.gitignore` and was never committed to the GitHub repository. However, GitGuardian may have detected it through other means (e.g., local scanning, cloud sync, or temporary exposure).

## Immediate Actions Required

### 1. Revoke the Exposed Token

**CRITICAL**: The exposed token must be revoked immediately:

1. Visit https://pypi.org/manage/account/token/
2. Delete all existing API tokens
3. Generate a new token for the `radia-ngsolve` package

### 2. Verify Package Security

Check PyPI for unauthorized uploads:

1. Visit https://pypi.org/project/radia-ngsolve/
2. Review recent releases for any unauthorized versions
3. If suspicious versions exist, contact PyPI support immediately

## Remediation Completed

The `Publish_to_PyPI.ps1` script has been updated to use environment variables instead of hardcoded tokens:

### Before (Insecure)
```powershell
python -m twine upload dist/*.whl dist/*.tar.gz --username __token__ --password pypi-AgEI...
```

### After (Secure)
```powershell
if (-not $env:PYPI_TOKEN) {
	Write-Host "ERROR: PYPI_TOKEN environment variable not set" -ForegroundColor Red
	exit 1
}

python -m twine upload dist/*.whl dist/*.tar.gz --username __token__ --password $env:PYPI_TOKEN
```

### New Usage

```powershell
# Set the API token as an environment variable
$env:PYPI_TOKEN = "your-new-token-here"

# Run the build and upload script
.\Publish_to_PyPI.ps1
```

## Documentation Updates

The following files have been updated with secure usage instructions:

1. **PYPI_DISTRIBUTION.md**: Added environment variable setup instructions
2. **docs/CLAUDE.md**: Updated PyPI workflow with secure token handling
3. **Publish_to_PyPI.ps1**: Now requires environment variable (no hardcoded tokens)

## Prevention Measures

1. **Never hardcode credentials**: Always use environment variables or secure credential stores
2. **Verify .gitignore**: Ensure sensitive files are properly excluded from Git
3. **Use scoped tokens**: Create package-specific tokens with minimal permissions
4. **Regular audits**: Periodically review scripts for hardcoded secrets

## Files Modified

- `Publish_to_PyPI.ps1`: Updated to use environment variables
- `PYPI_DISTRIBUTION.md`: Added security documentation
- `docs/CLAUDE.md`: Updated usage instructions
- `SECURITY_ALERT.md`: This document (created)

## Timeline

- **2025-11-06**: GitGuardian alert received
- **2025-11-06**: Issue investigated and root cause identified
- **2025-11-06**: Script updated to use environment variables
- **2025-11-06**: Documentation updated
- **Pending**: User to revoke exposed token and generate new one

## References

- PyPI Token Management: https://pypi.org/manage/account/token/
- PyPI Security: https://pypi.org/security/
- Best Practices: https://packaging.python.org/guides/publishing-package-distribution-releases-using-github-actions-ci-cd-workflows/

---

**Action Required**: Please revoke the exposed token immediately and generate a new one.
