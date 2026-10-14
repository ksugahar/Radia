#!/usr/bin/env python
"""
Radia Test Runner - Convenience script for running pytest tests

Usage:
    python run_tests.py                    # Run basic tests (fast)
    python run_tests.py --all              # Run all tests
    python run_tests.py --comprehensive    # Run comprehensive tests
    python run_tests.py --benchmark        # Run benchmark tests (slow)
    python run_tests.py -v                 # Verbose output
    python run_tests.py -k "test_import"   # Run specific test by name

This script provides a simple interface to pytest without needing to
remember command-line options.
"""

import sys
import os
import argparse

# Add tests directory to path for conftest.py to be found
tests_dir = os.path.dirname(os.path.abspath(__file__))
project_root = os.path.dirname(tests_dir)

def main():
    parser = argparse.ArgumentParser(
        description='Radia Test Runner',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_tests.py                    Run basic tests (fast)
  python run_tests.py --all              Run all tests
  python run_tests.py --comprehensive    Run comprehensive tests
  python run_tests.py --benchmark        Run benchmarks (slow)
  python run_tests.py -v                 Verbose output
  python run_tests.py -k "import"        Run tests matching "import"
  python run_tests.py --ngsolve          Run NGSolve integration tests
"""
    )

    # Test selection
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--basic', action='store_true', default=True,
                       help='Run basic tests only (default)')
    group.add_argument('--all', action='store_true',
                       help='Run all tests (may be slow)')
    group.add_argument('--comprehensive', action='store_true',
                       help='Run comprehensive tests')
    group.add_argument('--benchmark', action='store_true',
                       help='Run benchmark tests (slow)')
    group.add_argument('--ngsolve', action='store_true',
                       help='Run NGSolve integration tests')

    # pytest options passthrough
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Verbose output')
    parser.add_argument('-k', '--keyword', type=str, default=None,
                        help='Only run tests matching the keyword expression')
    parser.add_argument('-x', '--exitfirst', action='store_true',
                        help='Exit on first failure')
    parser.add_argument('--tb', type=str, default='short',
                        choices=['short', 'long', 'auto', 'no', 'line', 'native'],
                        help='Traceback style (default: short)')
    parser.add_argument('--pdb', action='store_true',
                        help='Drop into debugger on failures')
    parser.add_argument('--collect-only', action='store_true',
                        help='Only collect tests, do not run')

    args = parser.parse_args()

    # Check if pytest is available
    try:
        import pytest
    except ImportError:
        print("ERROR: pytest is not installed.")
        print("Install with: pip install pytest")
        return 1

    # Build pytest arguments
    pytest_args = [tests_dir]

    # Test selection markers
    if args.all:
        pass  # No marker filter
    elif args.comprehensive:
        pytest_args.extend(['-m', 'comprehensive'])
    elif args.benchmark:
        pytest_args.extend(['-m', 'benchmark'])
    elif args.ngsolve:
        pytest_args.extend(['-m', 'ngsolve'])
    else:  # default: basic
        pytest_args.extend(['-m', 'basic'])

    # Options
    if args.verbose:
        pytest_args.append('-v')

    if args.keyword:
        pytest_args.extend(['-k', args.keyword])

    if args.exitfirst:
        pytest_args.append('-x')

    pytest_args.extend(['--tb', args.tb])

    if args.pdb:
        pytest_args.append('--pdb')

    if args.collect_only:
        pytest_args.append('--collect-only')

    # Print what we're doing
    print("=" * 60)
    print("Radia Test Runner")
    print("=" * 60)
    print(f"Project root: {project_root}")
    print(f"Tests directory: {tests_dir}")
    print(f"pytest args: {' '.join(pytest_args)}")
    print("=" * 60)
    print()

    # Run pytest
    return pytest.main(pytest_args)


if __name__ == '__main__':
    sys.exit(main())
