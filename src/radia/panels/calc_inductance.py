"""
Inductance extractor using ngbem LaplaceSL BEM.

Called as subprocess from Cubit panel:
    python calc_inductance.py --cub5 model.cub5 --source src --sink snk
           --sigma 5.8e7 --order 3 --freq 0

Workflow:
  1. Import NGSolve FIRST (before cubit)
  2. Open cub5, export STEP
  3. Use BEMExtractor to compute L (and R at given frequency)

IMPORTANT: NGSolve must be imported BEFORE cubit to avoid numpy DLL conflict.
Outputs JSON to stdout (last line).
"""

import argparse
import json
import os
import sys
import tempfile


def _setup_cubit():
    """Import cubit with path cleanup (NGSolve already imported)."""
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))
    from install_panels import find_cubit_bin

    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)

    # Block Cubit's bundled site-packages
    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    import cubit
    # Suppress cubit.init() banner on stderr
    # Use -noinit to prevent .cubit startup file from being played
    # (it tries to load register_toolbar.py which imports Qt)
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

    for p in list(sys.path):
        if "site-packages" in p and ("cubit" in p.lower() or "Cubit" in p):
            sys.path.remove(p)

    return cubit


def extract_inductance(cub5_file, source_block, sink_block, sigma, order, freq):
    """Extract inductance using BEMExtractor.

    Args:
        cub5_file: Cubit model file
        source_block: Block name for current source face
        sink_block: Block name for current sink face (empty string = closed loop)
        sigma: Conductivity [S/m]
        order: Curve order for high-order BEM
        freq: Frequency [Hz] (0 = DC)
    """
    # 1. Import NGSolve FIRST
    from ngsolve import Mesh, TaskManager

    # 2. Import cubit
    cubit = _setup_cubit()

    # 3. Load model
    cubit.cmd(f'open "{cub5_file}"')

    # 4. Export STEP
    tmpdir = tempfile.mkdtemp(prefix="radia_ind_")
    step_file = os.path.join(tmpdir, "geometry.step").replace("\\", "/")
    cubit.cmd(f'export step "{step_file}" overwrite')

    # 5. Get block names for reporting
    radia_src = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    if os.path.abspath(radia_src) not in sys.path:
        sys.path.insert(0, os.path.abspath(radia_src))

    from cubit_bem_extractor import BEMExtractor

    # 6. Create extractor and load mesh
    ext = BEMExtractor()
    ext.load_cubit_mesh(cubit, geometry_file=step_file, curve_order=order)

    # 7. Set conductor
    # Find the conductor block (first volume block)
    block_ids = cubit.get_block_count()
    conductor_block = None
    for bid in range(1, block_ids + 1):
        try:
            name = cubit.get_exodus_entity_name("block", bid)
            if name and name not in [source_block, sink_block]:
                conductor_block = name
                break
        except Exception:
            pass

    if conductor_block is None:
        # Use default conductor label
        conductor_block = "conductor"

    ext.set_conductor(conductor_block, sigma=sigma)

    # 8. Set port
    if sink_block:
        # Open conductor: source -> sink
        ext.set_port("port1", block=conductor_block,
                     source_block=source_block, sink_block=sink_block)
    else:
        # Closed loop
        ext.set_port("port1", block=conductor_block)

    # 9. Extract
    freqs = [freq] if freq > 0 else [0.0]
    with TaskManager():
        result = ext.extract(freqs=freqs, mode='mqs')

    # 10. Build output
    L = result.L[0, 0] if hasattr(result.L, 'shape') else result.L
    R = result.R_at(freqs[0]) if freq > 0 else 0.0

    return {
        "inductance_H": float(L),
        "resistance_ohm": float(R) if R else 0.0,
        "frequency_Hz": freq,
        "sigma": sigma,
        "order": order,
        "source_block": source_block,
        "sink_block": sink_block or "(closed loop)",
    }


def main():
    parser = argparse.ArgumentParser(description="ngbem inductance extractor")
    parser.add_argument("--cub5", required=True, help="Cubit .cub5 model file")
    parser.add_argument("--source", default="", help="Source block name")
    parser.add_argument("--sink", default="", help="Sink block name")
    parser.add_argument("--sigma", type=float, default=5.8e7, help="Conductivity [S/m]")
    parser.add_argument("--order", type=int, default=1, help="Curve order (1-5)")
    parser.add_argument("--freq", type=float, default=0.0, help="Frequency [Hz]")
    args = parser.parse_args()

    try:
        result = extract_inductance(
            args.cub5, args.source, args.sink,
            args.sigma, args.order, args.freq
        )
        print(json.dumps(result))
    except Exception as e:
        print(json.dumps({"error": str(e)}))
        sys.exit(1)


if __name__ == "__main__":
    main()
