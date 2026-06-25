"""Build the ELF CEFC-2020 C-type yoke as Cubit CAD volumes.

Each of the 13 hex elements in a private CEFC-2020 C-type .meg
reference becomes ONE Cubit CAD volume.  4 of the 13 (hex 6, 7, 8, 9)
have slanted faces at the pole-tip bevel; the rest are axis-aligned
bricks.  Coordinates are mm; we scale to meters when emitting `create
vertex` commands.

After the volumes are created they are added to block "yoke" and
each is meshed with curve interval = 1 to produce exactly 1 hex per
volume -- matching the ELF discretization.

Played from `em_elf_quarter_xz.jou` (or run standalone with
RADIA_PANELS_DIR set, but typically driven by the .jou).
"""
from __future__ import annotations
import os
import cubit


ELF_MEG = os.environ.get("ELF_MEG_PATH", "")

MM_TO_M = 0.001


def _read_meg(path):
    """Parse ELF .meg into (nodes, hex_elements).

    nodes: dict node_id -> (x, y, z) in mm
    hex_elements: list of (eid, [n1..n8])
    """
    nodes = {}
    hex_elements = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if line.startswith("MGR1"):
                parts = line.split()
                nid = int(parts[1])
                x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
                nodes[nid] = (x, y, z)
            elif line.startswith("MMB8T"):
                parts = line.split()
                eid = int(parts[1])
                nids = [int(parts[i]) for i in range(4, 12)]
                hex_elements.append((eid, nids))
    return nodes, hex_elements


def _build_yoke():
    if not ELF_MEG or not os.path.isfile(ELF_MEG):
        print(f"ERROR: ELF_MEG not found: {ELF_MEG}")
        print("Set ELF_MEG_PATH env var to the private .meg reference file.")
        return

    nodes, elements = _read_meg(ELF_MEG)
    print(f"[em_elf_yoke_builder] ELF mesh: {len(nodes)} nodes, "
          f"{len(elements)} hex elements")

    # `create volume surface ... heal` consumes the input vertices /
    # surfaces -- so vertices CAN'T be shared across hexes via a
    # global elf_to_v map.  Strategy: create FRESH vertices per hex
    # (8 vertices each, 13 * 8 = 104 vertices total, with duplicates
    # at shared positions), build 6 surfaces + 1 volume per hex,
    # then `imprint+merge volume all` at the end to fuse coincident
    # vertices/curves so adjacent volumes share faces.
    yoke_vols = []
    # Exodus hex8 face ordering (bottom CCW from outside, top
    # matching CCW so side quads stay consistent):
    face_pairs = [
        (0, 1, 2, 3),   # bottom
        (4, 5, 6, 7),   # top
        (0, 1, 5, 4),   # side 1
        (1, 2, 6, 5),   # side 2
        (2, 3, 7, 6),   # side 3
        (3, 0, 4, 7),   # side 4
    ]
    for eid, nids in elements:
        # Create 8 fresh vertices at the ELF node positions for THIS hex.
        cv = []
        for n in nids:
            x, y, z = nodes[n]
            cubit.cmd(f"create vertex {x * MM_TO_M:.10g} "
                      f"{y * MM_TO_M:.10g} {z * MM_TO_M:.10g}")
            cv.append(cubit.get_last_id("vertex"))
        face_sids = []
        for fp in face_pairs:
            v_str = " ".join(str(cv[i]) for i in fp)
            cubit.cmd(f"create surface vertex {v_str}")
            face_sids.append(cubit.get_last_id("surface"))
        s_str = " ".join(str(s) for s in face_sids)
        cubit.cmd(f"create volume surface {s_str} heal")
        vid = cubit.get_last_id("volume")
        yoke_vols.append(vid)

    # Unite the 13 yoke volumes into ONE CAD volume so that:
    #   - yoke-yoke interfaces are conformal by construction
    #   - air-yoke imprint+merge produces a single shared boundary
    #   - subsequent CAD subtraction (in the .jou) carves a single
    #     yoke-shaped hole from the air sphere
    # Cost: lose the EXACT 13-hex topology -- the unioned yoke
    # tetmeshes (or hex-meshes if topology allows) at the user's
    # chosen mesh size, NOT the per-region 1-hex-each ELF
    # discretization.  Trade-off: physics fidelity (correct B at
    # origin) vs mesh-count fidelity (exactly 13 hex).
    v_str = " ".join(str(v) for v in yoke_vols)
    cubit.cmd(f"unite volume {v_str}")
    yoke_united = cubit.get_last_id("volume")
    # Tetmesh with mesh size targeting the smallest ELF region
    # (z=5..13mm bevel layer = 8 mm thick).  Use 5 mm for adequate
    # discretization through that layer.
    cubit.cmd(f"volume {yoke_united} scheme tetmesh")
    # Empirical (2026-04-25): 5mm yoke gives Bz=-240mT (5.2% off ELF
    # -228.1mT).  Refining to 2mm gives Bz=-177mT (DIVERGES, -22%);
    # higher fes_order also diverges (Kelvin pair ID is linear-only).
    # 5mm + fes_order=1 is the empirical sweet spot for this geometry.
    cubit.cmd(f"volume {yoke_united} size 0.005")
    cubit.cmd(f"mesh volume {yoke_united}")
    yoke_vols = [yoke_united]
    v_str = str(yoke_united)

    # Block assignment.
    cubit.cmd(f"block 1 add volume {v_str}")
    cubit.cmd('block 1 name "yoke"')

    n_hex = len(list(cubit.parse_cubit_list("hex", "in volume " + v_str)))
    print(f"[em_elf_yoke_builder] yoke meshed: {len(yoke_vols)} volumes, "
          f"{n_hex} hex elements (target: 13)")


_build_yoke()
