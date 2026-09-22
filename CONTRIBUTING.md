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
- [ ] No `CblasColMajor` in core source (except LAPACK wrapper boundaries)
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

These are enforced by CI ([policy-lint.yml](.github/workflows/policy-lint.yml))
through `tools/policy_lint.py`, which is the single source of truth for the
list. A check that cannot run (git failing, a missing tool) fails the lint
rather than passing it. At the time of writing:

1. **Unit System**: Radia always uses meters. No `FldUnits()` call needed
2. **No Binaries in Git**: .pyd/.dll/.so hosted on GitHub Releases
3. **Laplace Kernel Only**: No Helmholtz in C++ core
4. **Row-Major Matrices**: No `CblasColMajor` in core (except LAPACK wrapper boundaries)
5. **No Generated Files at Root**: .msh/.vtu/.vtk/.vol go next to source scripts
6. **No Legacy Paths**: Use `src/radia`, not `src/python`
7. **Examples Tier Retired**: no tracked `examples/`
8. **Central HDiv Capabilities**: geometry/field order pairs come from one table
9. **No Product Names on Formulations**: cite the paper, not the vendor
10. **Commercial Interchange Stays Non-Public**: converters and their MCP
    servers live in the tool's non-public home, never here

### Examples Retired Policy

`examples/` is retired and must not be recreated. New exploratory work starts
outside the repository in `C:\temp`; tracked work enters only after promotion
to a durable lane:

- `tests/` for fast, deterministic bug and contract protection
- `validation_test/` for numerical, performance, native and multi-machine
  evidence with machine-readable result JSON
- `docs/<topic>/*.ipynb` for executed, result-bearing public demonstrations.
  Saved notebook output is the record; no adjacent JSON or runtime gate is
  required
- `src/` for reusable APIs and solver helpers

Standalone PySide/PyQt panels and notebook workbenches are retired and are not
a promotion target. The human production interface is the masked blocks in the
single Radia Simulink library; the AI interface is Python/MCP.

Historical references to `examples/` are migration blockers. Migrate the
owning code/docs to one of the lanes above, then delete the stale reference.

## Reporting Issues

- Use the [bug report template](.github/ISSUE_TEMPLATE/bug_report.yml) for bugs
- Use the [feature request template](.github/ISSUE_TEMPLATE/feature_request.yml) for enhancements

## License

By contributing, you agree that your contributions will be licensed under the terms described in [LICENSE](LICENSE).
