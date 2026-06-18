# Contributing to Radia

Thank you for your interest in contributing to Radia! This document provides guidelines and information for contributors.

## Getting Started

1. Fork the repository
2. Clone your fork and create a feature branch from `main`
3. Set up the development environment (see [BUILD.md](BUILD.md))

## Development Environment

**Requirements**:
- Python 3.12+
- MSVC (Visual Studio 2022) + Intel MKL
- CMake 3.20+

**Install dev dependencies**:
```bash
pip install -e ".[dev]"
```

## Coding Standards

### Python
- Formatter: **black** (line-length 100)
- Import sorting: **isort** (black profile)
- Target: Python 3.12

### C++
- Compiler: MSVC only (no GCC/Clang support currently)
- BLAS/LAPACK: Intel MKL only
- Matrix storage: **Row-major** `[target][source]` convention
- Green's function: **Laplace kernel only** (no Helmholtz)

### Commit Messages
- Use imperative mood: "Add feature" not "Added feature"
- Keep the first line under 72 characters
- Reference issues when applicable: "Fix #123"

## Pull Request Process

1. Ensure CI passes (build + basic tests)
2. Update `CHANGELOG.md` under `[Unreleased]` if your change is user-facing
3. Follow the PR template checklist
4. Request review from a maintainer

### PR Checklist
- [ ] Code follows the project coding standards
- [ ] No binary files (.pyd, .dll, .so) committed to git
- [ ] No `rad.FldUnits()` calls (removed; Radia always uses meters)
- [ ] No `CblasColMajor` in core source (except MMM/LAPACK)
- [ ] Tests pass locally: `pytest tests/ -m basic`

## Testing

```bash
# Run basic tests (fast, ~2 sec)
pytest tests/ -m basic -v

# Run comprehensive tests
pytest tests/ -m comprehensive -v

# Run all tests with coverage
pytest tests/ --cov=radia --cov-report=html
```

Test markers: `basic`, `comprehensive`, `advanced`, `performance`, `slow`, `benchmark`

## Project Policies

These are enforced by CI ([policy-lint.yml](.github/workflows/policy-lint.yml)):

1. **Unit System**: Radia always uses meters. No `FldUnits()` call needed
2. **No Binaries in Git**: .pyd/.dll/.so hosted on GitHub Releases
3. **Laplace Kernel Only**: No Helmholtz in C++ core
4. **Row-Major Matrices**: No `CblasColMajor` in core (except MMM/LAPACK)
5. **No Generated Files at Root**: .msh/.vtu/.vtk/.vol go next to source scripts
6. **No Legacy Paths**: Use `src/radia`, not `src/python`

### Example README Policy

Every example directory under `examples/` must have exactly **one README.md** that lists all Python files with descriptions.

- **Single-topic folder** (e.g., `vtk_export/`): Place `README.md` directly in the folder.
- **Multi-topic folder** (e.g., `kelvin_transformation/`): Place one top-level `README.md` that summarizes all sub-themes, plus optionally one `README.md` per subfolder.
- **README must match code**: Every `.py` file on disk must be listed in the README. Remove entries for deleted files. Update the README when adding or removing scripts.
- **Keep READMEs concise**: 1-line description per script, plus a brief overview of the folder's purpose.

## Reporting Issues

- Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.yml) for bugs
- Use the [feature request template](.github/ISSUE_TEMPLATE/feature_request.yml) for enhancements

## License

By contributing, you agree that your contributions will be licensed under the terms described in [LICENSE](LICENSE).
