# Radia Test Suite

Automated tests for the Radia Python 3.12 module.

## Directory Structure

```
tests/
├── __init__.py                      # Test package initialization
├── test_simple.py                   # Basic functionality tests (quick)
├── test_radia.py                    # Comprehensive test suite
├── test_advanced.py                 # Advanced features and edge cases
├── test_parallel_performance.py     # OpenMP parallelization tests
├── test_radia_ngsolve.py              # NGSolve integration tests
├── test_magpylib_comparison.py      # Comparison with magpylib library
├── conftest.py                      # pytest configuration
├── test_utils.py                    # Shared path utilities
├── benchmarks/                      # Performance benchmarks
│   ├── benchmark_openmp.py         # OpenMP scaling tests
│   ├── benchmark_correct.py        # Correctness vs performance
│   ├── benchmark_heavy.py          # Heavy computation tests
│   └── benchmark_threads.py        # Thread scaling tests
├── fixtures/                        # Test data and helper functions
└── README.md                        # This file
```

## Running Tests

### Prerequisites

```bash
# Ensure radia module is built and available
cd S:/Visual_Studio/02_Visual_Studio_2022_コマンドライン_コンパイル/04_Radia
powershell.exe -ExecutionPolicy Bypass -File Build.ps1

# Install test dependencies (optional)
pip install pytest pytest-cov
```

### Quick Test (Recommended for CI/CD)

```bash
# Run basic functionality test (fastest)
python tests/test_simple.py
```

**Expected output:**
```
============================================================
[OK] ALL TESTS PASSED!
============================================================
```

### Comprehensive Test Suite

```bash
# Run all functional tests
python tests/test_radia.py
```

**Expected output:**
```
*** ALL TESTS PASSED! ***
Total: 7/7 tests passed (100.0%)
```

### Advanced Tests

```bash
# Run advanced feature tests
python tests/test_advanced.py
```

### Performance Tests

```bash
# Test OpenMP parallelization performance
python tests/test_parallel_performance.py
```

### NGSolve Integration Test

```bash
# Test radia_ngsolve module (requires NGSolve)
python tests/test_radia_ngsolve.py
```

**Note:** This test requires NGSolve to be installed. If NGSolve is not available, the test will be skipped.

**Expected output:**
```
[OK] ALL TESTS PASSED!
```

### Comparison Test with magpylib

Compare Radia results with magpylib (another magnetic field library):

```bash
# Requires: pip install magpylib
python tests/test_magpylib_comparison.py
```

**Note**: This test validates Radia's magnetic field calculations by comparing with magpylib (an independent library). Key points:
- Both libraries use **SI units (Tesla)** for magnetization/polarization and field output
- For permanent magnets: M = Br (remanence)
- Agreement: **~0.5%** with 32 subdivisions (default setting)

**Subdivision convergence** (azimuthal segments):
- 16 segments: ~2% error
- 32 segments: ~0.5% error (current default)
- 64 segments: ~0.1% error
- 128 segments: ~0.03% error

This validates the accuracy of Radia's electromagnetic field calculations and demonstrates
second-order convergence with mesh refinement.

## Using pytest (Optional)

If you have pytest installed:

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest -v tests/

# Run specific test file
pytest tests/test_simple.py

# Run with coverage report
pytest --cov=radia tests/

# Run only tests matching pattern
pytest -k "material" tests/
```

## Benchmarks

Performance benchmarking scripts are in `tests/benchmarks/`:

### OpenMP Scaling Benchmark

```bash
python tests/benchmarks/benchmark_openmp.py
```

Tests field calculation performance with different thread counts (1, 2, 4, 8 cores).

**Expected speedup (8-core system):**
- Simple geometry: 1.8x
- Complex geometry: 2.7x

### Other Benchmarks

```bash
# Correctness validation
python tests/benchmarks/benchmark_correct.py

# Heavy computation benchmark
python tests/benchmarks/benchmark_heavy.py

# Thread scaling analysis
python tests/benchmarks/benchmark_threads.py
```

## Test Categories

### 1. Basic Functionality (`test_simple.py`)

Fast smoke tests covering:
- Module import
- Version checking
- Basic geometry creation
- Magnetization setting
- Field calculation
- Object deletion

**Duration**: ~1 second

### 2. Comprehensive Suite (`test_radia.py`)

Complete functional test suite:
- Module import and version
- Geometry creation (rectangles, polygons)
- Material definition (Steel37, NdFeB, etc.)
- Field calculations (B, H, A)
- Relaxation solver
- Geometric transformations
- Memory management

**Duration**: ~5-10 seconds

### 3. Advanced Features (`test_advanced.py`)

Edge cases and advanced functionality:
- Complex geometries
- Subdivision algorithms
- Non-linear materials
- Multiple object interactions
- VTK export functionality

**Duration**: ~30-60 seconds

### 4. Parallel Performance (`test_parallel_performance.py`)

OpenMP parallelization validation:
- Single-threaded baseline
- Multi-threaded scaling
- Performance regression detection
- Numerical consistency check

**Duration**: ~2-5 minutes

## Continuous Integration

### Minimal CI Configuration

For fast CI/CD pipelines, run only basic tests:

```yaml
# .github/workflows/test.yml (example)
- name: Run basic tests
  run: python tests/test_simple.py
```

### Full CI Configuration

For comprehensive validation:

```yaml
- name: Run all tests
  run: |
	python tests/test_simple.py
	python tests/test_radia.py
	python tests/test_parallel_performance.py
```

## Test Data

Test fixtures and helper data should be placed in `tests/fixtures/`:

```python
# Example usage
from tests.fixtures.geometries import create_test_magnet

magnet = create_test_magnet()
```

## Writing New Tests

### Test File Naming Convention

- `test_*.py`: Functional tests (auto-discovered by pytest)
- `benchmark_*.py`: Performance benchmarks (manual execution)

### Test Function Naming

```python
def test_feature_name():
	"""Test description"""
	# Arrange
	# Act
	# Assert
```

### Example Test

```python
def test_magnet_field_calculation():
	"""Test that rectangular magnet produces correct field"""
	import sys
	sys.path.insert(0, 'build/lib/Release')
	import radia as rad

	# Create magnet
	mag = rad.ObjRecMag([0,0,0], [10,10,10])
	rad.ObjSetM(mag, [0,0,1000])

	# Calculate field
	field = rad.Fld(mag, 'b', [0,0,20])

	# Verify
	assert field[2] > 0, "Field should be positive in Z direction"
	assert abs(field[0]) < 0.01, "Bx should be near zero"
	assert abs(field[1]) < 0.01, "By should be near zero"
```

## Troubleshooting

### Module Import Fails

```python
# Error: No module named 'radia'
```

**Solution**: Ensure module is built and in correct location:
```bash
ls build/lib/Release/radia.pyd  # Should exist
```

### Tests Fail After Code Changes

1. Rebuild module: `powershell.exe -ExecutionPolicy Bypass -File Build.ps1`
2. Clear Python cache: `find . -name "*.pyc" -delete`
3. Re-run tests

### Performance Tests Show Regression

Check:
- Build configuration (should be Release, not Debug)
- OpenMP enabled: Check CMakeLists.txt
- System load: Close other applications
- Thread count: Set `OMP_NUM_THREADS` environment variable

## Test Coverage Goals

- **Unit tests**: >80% code coverage
- **Integration tests**: All major API functions
- **Performance tests**: No regression >10%
- **Memory tests**: No leaks detected

## Reporting Issues

If tests fail, please report with:
1. Test output (full error message)
2. System info: OS, Python version, compiler
3. Build configuration: Debug/Release
4. Steps to reproduce

## References

- pytest documentation: https://docs.pytest.org/
- Radia documentation: See `docs/` directory
- Performance reports: `docs/OPENMP_PERFORMANCE_REPORT.md`
- Security fixes: `SECURITY_FIXES.md`

---

**Last Updated**: 2025-10-30
**Radia Version**: 4.32
**Python Version**: 3.12
