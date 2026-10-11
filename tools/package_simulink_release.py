"""Reject retired CPU MEX ZIP publication; use the independent GPU wheel."""
from pathlib import Path
import argparse


def build_package(mex_dir: Path, output_dir: Path, *, full_library: bool = False):
    raise RuntimeError(
        "CPU MEX ZIP publication is retired. Build the independent radia-gpu "
        "CUDA wheel through packages/radia-gpu; CPU execution uses Python/pybind11."
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mex-dir', required=True, type=Path)
    parser.add_argument('--output-dir', type=Path, default=Path('dist'))
    parser.add_argument('--full-library', action='store_true')
    args = parser.parse_args()
    build_package(args.mex_dir, args.output_dir, full_library=args.full_library)


if __name__ == '__main__':
    main()
