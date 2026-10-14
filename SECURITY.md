# Security Policy

## Supported Versions

| Version | Supported          |
|---------|--------------------|
| 4.x     | Yes                |
| 3.x and earlier | No        |

## Reporting a Vulnerability

If you discover a security vulnerability in Radia, please report it responsibly.

**Do NOT open a public issue.**

Instead, contact the maintainers privately:
- **Contact**: Open a [private security advisory](https://github.com/ksugahar/Radia/security/advisories/new) on GitHub

### What to Include

- Description of the vulnerability
- Steps to reproduce
- Potential impact
- Suggested fix (if any)

### Response Timeline

- **Acknowledgement**: Within 7 days
- **Assessment**: Within 14 days
- **Fix release**: Depends on severity

## Scope

This policy covers:
- The Radia Python package (`radia`)
- The `radia-mcp` and `cubit-mesh-export` packages in this monorepo
- C++ core library (`src/core/`)
- Build scripts and CI/CD workflows
- Example scripts (if they demonstrate insecure patterns)

## Known Considerations

- Radia uses Intel MKL shared libraries (`mkl_rt.dll`). MKL internally depends on `libiomp5md.dll` (Intel OpenMP). Ensure these are obtained from official Intel channels.
- The `_radia_pybind.pyd` binary is distributed via GitHub Releases. Verify checksums when downloading pre-built binaries.
