# Radia Documentation

This folder contains the official documentation for Radia.

## Documentation Organization

### User Documentation (Essential Reading)

#### API Documentation
- [API_REFERENCE.md](API_REFERENCE.md) - Complete Python API reference
- [API_EXTENSIONS.md](API_EXTENSIONS.md) - Extended features and new APIs

#### H-Matrix Acceleration
- [HMATRIX_USER_GUIDE.md](HMATRIX_USER_GUIDE.md) - User guide for H-matrix acceleration
- [HMATRIX_BENCHMARKS_RESULTS.md](HMATRIX_BENCHMARKS_RESULTS.md) - Performance benchmarks
- [HMATRIX_SERIALIZATION.md](HMATRIX_SERIALIZATION.md) - Disk cache feature

#### NGSolve Integration
- [NGSOLVE_USAGE_GUIDE.md](NGSOLVE_USAGE_GUIDE.md) - How to use Radia with NGSolve
- [NGSOLVE_INTEGRATION.md](NGSOLVE_INTEGRATION.md) - Integration overview

### Developer Documentation (Implementation Details)

#### H-Matrix Implementation
- [HMATRIX_IMPLEMENTATION_HISTORY.md](HMATRIX_IMPLEMENTATION_HISTORY.md) - Development history
- [HMATRIX_ENHANCEMENT_PROPOSAL_2025.md](HMATRIX_ENHANCEMENT_PROPOSAL_2025.md) - Future enhancements
- [hmatrix_field_design.md](hmatrix_field_design.md) - Field evaluation design
- [ML_PARAMETER_TUNING.md](ML_PARAMETER_TUNING.md) - Machine learning parameter optimization

#### Material and Field APIs
- [MATERIAL_API_IMPLEMENTATION.md](MATERIAL_API_IMPLEMENTATION.md) - Material property implementation
- [CF_BACKGROUND_FIELD_IMPLEMENTATION.md](CF_BACKGROUND_FIELD_IMPLEMENTATION.md) - Background field implementation
- [NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md](NGSOLVE_CF_BACKGROUND_FIELD_DESIGN.md) - NGSolve CoefficientFunction design

## Additional Resources

### Development Notes
Implementation notes, performance analysis, and troubleshooting are in [../dev/notes/](../dev/notes/):
- Implementation details (directory structure, conversion reports)
- Performance analysis (OpenMP, H-matrix optimizations)
- Issue resolutions (DLL issues, integration problems)

### Release History
Release notes and version history are in [../dev/notes/releases/](../dev/notes/releases/).

### Examples
Working code examples are in [../examples/](../examples/):
- [simple_problems/](../examples/simple_problems/) - Basic examples
- [solver_benchmarks/](../examples/solver_benchmarks/) - Performance benchmarks
- [ngsolve_integration/](../examples/ngsolve_integration/) - NGSolve integration examples

## Quick Start

1. **New users**: Start with [API_REFERENCE.md](API_REFERENCE.md) and [examples/simple_problems/](../examples/simple_problems/)
2. **Performance optimization**: Read [HMATRIX_USER_GUIDE.md](HMATRIX_USER_GUIDE.md)
3. **NGSolve users**: Read [NGSOLVE_USAGE_GUIDE.md](NGSOLVE_USAGE_GUIDE.md)
4. **Contributors**: Review [HMATRIX_IMPLEMENTATION_HISTORY.md](HMATRIX_IMPLEMENTATION_HISTORY.md)

## See Also

- [Main README](../README.md) - Project overview and installation
- [Build Instructions](../README_BUILD.md) - How to build from source
- [CHANGELOG](../CHANGELOG.md) - Version history
