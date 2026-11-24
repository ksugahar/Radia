"""Display of both magnetic-surface families on the two-hex analytic case.

Display only. Reuses two_families.py for the field, projection, cell-local
tracker, patch flux and surface flux, so the numbers in the JSON are the same
diagnostics as the convergence study, evaluated on the displayed surfaces.

The input is the divergence-free B_h = curl(A_h) copied into HDiv; no H is
involved. Each surface of family u (v) is the set of B_h lines seeded on the
inlet section z = 0 along x = label (y = label); its points come from the
tracker, which lands every crossing of the shared face x = 0 on that face and
continues from the same physical point in the neighbouring cell, so the
surfaces are connected across the face by construction.

Written: <out>.msh (MSH 4.2 ASCII, for Gmsh) and <out>.json.
  physical 2D "u_family_surfaces", "v_family_surfaces", "shared_face_x0"
  physical 1D "u_field_lines", "v_field_lines", "element_edges"
  NodeData   "distance_to_exact_surface" (oracle labels; position error of the
             displayed surface: projection + tracking + reconstruction)
  ElementData "family_label" and "local_leakage |Bh.n|/|Bh|" at quad centres
             (display only; the integrated leakage is in the JSON)

Usage: python display_surfaces.py --out display [--config 2x1x1:6]
"""
import argparse
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np
import ngsolve as ngs

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import two_families as tf  # noqa: E402


def write_msh42(path, groups, entities, node_blocks, elem_blocks, node_data, elem_data):
    """Minimal MSH 4.2 ASCII writer (same layout as radia.gmsh_post_export)."""
    with open(path, "w", newline="\n") as f:
        f.write("$MeshFormat\n4.2 0 8\n$EndMeshFormat\n")
        f.write("$PhysicalNames\n%d\n" % len(groups))
        for dim, tag, name in groups:
            f.write('%d %d "%s"\n' % (dim, tag, name))
        f.write("$EndPhysicalNames\n")
        cnt = [sum(1 for e in entities if e[0] == d) for d in range(4)]
        f.write("$Entities\n%d %d %d %d\n" % tuple(cnt))
        for d in range(4):
            for dim, tag, bbox, phys in entities:
                if dim != d:
                    continue
                ph = "1 %d" % phys if phys else "0"
                if dim == 0:  # point: tag x y z nphys [phys]
                    f.write("%d %s %s\n" % (tag, " ".join("%.15e" % v for v in bbox[:3]), ph))
                else:  # tag bbox nphys [phys] nbound [..] nembedded (MSH 4.2)
                    f.write("%d %s %s 0 0\n" % (tag, " ".join("%.15e" % v for v in bbox), ph))
        f.write("$EndEntities\n")
        nn = sum(len(c) for _, _, _, c in node_blocks)
        f.write("$Nodes\n%d %d 1 %d\n" % (len(node_blocks), nn, nn))
        for dim, tag, start, coords in node_blocks:
            f.write("%d %d 0 %d\n" % (dim, tag, len(coords)))
            f.write("".join("%d\n" % (start + i) for i in range(len(coords))))
            f.write("".join("%.15e %.15e %.15e\n" % tuple(c) for c in coords))
        f.write("$EndNodes\n")
        ne = sum(len(c) for _, _, _, c in elem_blocks)
        f.write("$Elements\n%d %d 1 %d\n" % (len(elem_blocks), ne, ne))
        eid, ids = 1, []
        for dim, tag, etype, conn in elem_blocks:
            f.write("%d %d %d %d\n" % (dim, tag, etype, len(conn)))
            ids.append(list(range(eid, eid + len(conn))))
            for c in conn:
                f.write("%d %s\n" % (eid, " ".join(str(int(n)) for n in c)))
                eid += 1
        f.write("$EndElements\n")
        for kind, name, items in [("$NodeData", n, it) for n, it in node_data] + [("$ElementData", n, it) for n, it in elem_data]:
            tags = [t for t, _ in items]
            f.write('%s\n1\n"%s"\n1\n0.0\n3\n0\n1\n%d\n' % (kind, name, len(tags)))
            f.write("1 %d %d %d\n%d\n" % (len(tags), min(tags), max(tags), len(tags)))
            f.write("".join("%d\n" % t for t in tags))
            f.write("".join("%.15e\n" % v for _, v in items))
            f.write(kind.replace("$", "$End") + "\n")
    return ids


def refine(pos, der, z_st, nsub):
    """Hermite along each tracked curve: (n_z_fine, Ns, 3) points."""
    out = [np.column_stack([pos[0], np.full(len(pos[0]), z_st[0])])]
    for k in range(len(z_st) - 1):
        dz = z_st[k + 1] - z_st[k]
        for i in range(1, nsub + 1):
            tau = i / nsub
            hv, _ = tf.hermite(pos[k], pos[k + 1], der[k], der[k + 1], dz, tau)
            out.append(np.column_stack([hv, np.full(len(hv), z_st[k] + tau * dz)]))
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output stem (writes .msh and .json)")
    ap.add_argument("--config", default="2x1x1:6")
    ap.add_argument("--labels", type=float, nargs="+", default=[-1.2, -0.6, 0.0, 0.6, 1.2])
    ap.add_argument("--seed-offset", type=float, default=0.0123)
    ap.add_argument("--ns", type=int, default=33)
    ap.add_argument("--mz", type=int, default=8, help="tracker output stations per z cell")
    ap.add_argument("--nsub", type=int, default=4, help="Hermite display points per station interval")
    ap.add_argument("--q", type=int, default=6)
    ap.add_argument("--loop-m", type=int, default=128)
    ap.add_argument("--hmax", type=float, default=0.05)
    ap.add_argument("--threads", type=int, default=1)
    ap.add_argument("--eps", type=float, default=1e-9)
    a = ap.parse_args()
    ngs.SetNumThreads(a.threads)
    t0 = time.perf_counter()
    dims, k = a.config.split(":")
    nx, ny, nz = [int(d) for d in dims.split("x")]
    k = int(k)
    with ngs.TaskManager():
        mesh = tf.build_mesh(nx, ny, nz)
        gfa, bh, pdiag = tf.project(mesh, k, a.eps)
    fx = np.linspace(-tf.HALF, tf.HALF, nx + 1)
    fy = np.linspace(-tf.HALF, tf.HALF, ny + 1)
    fz = np.linspace(0, tf.ZEND, nz + 1)
    field = tf.Field(mesh, bh, fx, fy, fz)
    afield = tf.Field(mesh, gfa, fx, fy, fz)
    bfun_h = lambda S: field(S[:, 0], S[:, 1], S[:, 2])
    bfun_x = lambda S: tf.exact_b(S[:, 0], S[:, 1], S[:, 2])

    # signed patch flux at the inlet, an interior cut and the outlet (Stokes, A_h)
    po = a.seed_offset
    patches = [(x0 + po, y0 + po) for x0 in np.arange(-1, 1, 0.5) for y0 in np.arange(-1, 1, 0.5)]
    edges = [tf.loop_points(x0, y0, x0 + .5, y0 + .5, a.loop_m) for x0, y0 in patches]
    zp = np.unique(np.concatenate([fz, [0.5 * tf.ZEND]]))
    imid = int(np.argmin(np.abs(zp - 0.5 * tf.ZEND)))
    with ngs.TaskManager():
        P, _, pcr = tf.track(field, np.vstack([e for es in edges for e in es]), zp, a.hmax, check_faces=False)
    prow, off = [], 0
    for (x0, y0), es in zip(patches, edges):
        sz = [len(e) for e in es]
        r = {"x0": x0, "y0": y0}
        for si, key in ((0, "inlet"), (imid, "mid"), (len(zp) - 1, "outlet")):
            ee, o = [], off
            for s_ in sz:
                ee.append(P[si][o:o + s_])
                o += s_
            r[key] = tf.loop_flux(afield, ee, zp[si])
        r["mid_rel"], r["outlet_rel"] = r["mid"] / r["inlet"] - 1, r["outlet"] / r["inlet"] - 1
        off += sum(sz)
        prow.append(r)
    patch = {"loop_points_per_edge": a.loop_m, "z_mid": float(zp[imid]), "area_exact": 0.25, "crossings": pcr["n"],
             "mid_rel_max": max(abs(r["mid_rel"]) for r in prow), "outlet_rel_max": max(abs(r["outlet_rel"]) for r in prow),
             "sum_inlet": sum(r["inlet"] for r in prow), "sum_mid": sum(r["mid"] for r in prow), "patches": prow}

    # both families
    z_st = np.linspace(0, tf.ZEND, nz * a.mz + 1)
    groups = [(2, 1, "u_family_surfaces"), (2, 2, "v_family_surfaces"), (2, 3, "shared_face_x0"),
              (1, 4, "u_field_lines"), (1, 5, "v_field_lines"), (1, 6, "element_edges")]
    entities, nblocks, eblocks = [], [], []
    node_dist, elem_label, elem_leak, surf_rec = [], [], [], []
    nid, etag = 1, 1
    quad_owner = []  # (block index, label, centres, local leakage)
    for fam in (0, 1):
        for lab in a.labels:
            label0 = lab + po
            tseed, X0 = tf.seeds_family(fam, label0, a.ns)
            tseed = tseed + po
            X0[:, 1 - fam] += po
            with ngs.TaskManager():
                pos, der, cr = tf.track(field, X0, z_st, a.hmax)
            uu, vv = (np.full(a.ns, label0), tseed) if fam == 0 else (tseed, np.full(a.ns, label0))
            pob = np.array([np.column_stack(tf.forward_map(uu, vv, zk)) for zk in z_st])
            dob = np.array([tf.exact_b(p_[:, 0], p_[:, 1], np.full(a.ns, zk))[:, :2] for p_, zk in zip(pob, z_st)])
            rec = {"family": "uv"[fam], "label": label0, "n_seed": a.ns, "stations": len(z_st), "q": a.q,
                   "tracked_Bh": tf.surface_flux(pos, der, z_st, tseed, bfun_h, fam, label0, a.q),
                   "oracle_exactB": tf.surface_flux(pob, dob, z_st, tseed, bfun_x, fam, label0, a.q),
                   "oracle_Bh": tf.surface_flux(pob, dob, z_st, tseed, bfun_h, fam, label0, a.q),
                   "shared_face_crossings": cr["n"], "normal_jump_max_rel": cr["normal_jump_max_rel"],
                   "tangent_jump_max_rel": cr["tangent_jump_max_rel"], "landing_residual_max": cr["landing_residual_max"]}
            surf_rec.append(rec)
            G = refine(pos, der, z_st, a.nsub)  # (nzf, ns, 3)
            nzf = G.shape[0]
            coords = G.reshape(-1, 3)
            ids = nid + np.arange(len(coords)).reshape(nzf, a.ns)
            labs = tf.labels(coords[:, 0], coords[:, 1], coords[:, 2])[fam]
            gr = np.linalg.norm(tf.label_grad(coords[:, 0], coords[:, 1], coords[:, 2], fam), axis=1)
            node_dist += list(zip(ids.ravel().tolist(), (np.abs(labs - label0) / gr).tolist()))
            bb = list(coords.min(0)) + list(coords.max(0))
            entities.append((2, etag, bb, 1 + fam))
            nblocks.append((2, etag, nid, coords))
            quads = np.stack([ids[:-1, :-1], ids[:-1, 1:], ids[1:, 1:], ids[1:, :-1]], -1).reshape(-1, 4)
            eblocks.append((2, etag, 3, quads))
            c = coords[quads - nid]  # (nq, 4, 3)
            cen = c.mean(1)
            nrm = np.cross(c[:, 2] - c[:, 0], c[:, 3] - c[:, 1])
            bc = bfun_h(cen)
            leak = np.abs(np.sum(bc * nrm, 1)) / (np.linalg.norm(bc, axis=1) * np.linalg.norm(nrm, axis=1))
            quad_owner.append((len(eblocks) - 1, label0, leak))
            nid += len(coords)
            etag += 1
            # field lines of this surface (one line per seed)
            lines = np.stack([ids[:-1, :], ids[1:, :]], -1).transpose(1, 0, 2).reshape(-1, 2)
            entities.append((1, etag, bb, 4 + fam))
            eblocks.append((1, etag, 1, lines))
            etag += 1
    # shared face and element edges
    corners = np.array([[x, y, z] for z in (0, tf.ZEND) for y in (-tf.HALF, tf.HALF) for x in fx])
    cid = nid + np.arange(len(corners))
    ix = lambda i, j, kz: cid[kz * 2 * (nx + 1) + j * (nx + 1) + i]
    shared = [[ix(i, 0, 0), ix(i, 1, 0), ix(i, 1, 1), ix(i, 0, 1)] for i in range(1, nx)]
    edg = []
    for kz in (0, 1):
        for j in (0, 1):
            edg += [[ix(i, j, kz), ix(i + 1, j, kz)] for i in range(nx)]
        edg += [[ix(i, 0, kz), ix(i, 1, kz)] for i in range(nx + 1)]
    edg += [[ix(i, j, 0), ix(i, j, 1)] for i in range(nx + 1) for j in (0, 1)]
    bbox = [-tf.HALF, -tf.HALF, 0.0, tf.HALF, tf.HALF, tf.ZEND]
    nblocks.append((0, etag, nid, corners))
    entities.append((0, etag, bbox, 0))
    etag += 1
    if shared:
        entities.append((2, etag, bbox, 3))
        eblocks.append((2, etag, 3, np.array(shared)))
        etag += 1
    entities.append((1, etag, bbox, 6))
    eblocks.append((1, etag, 1, np.array(edg)))

    msh = a.out + ".msh"
    ids = write_msh42_entities(msh, groups, entities, nblocks, eblocks, quad_owner, node_dist)
    res = {"script": os.path.basename(__file__),
           "script_sha256": hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest(),
           "two_families_sha256": hashlib.sha256(open(tf.__file__, "rb").read()).hexdigest(),
           "args": {**vars(a), "out": os.path.basename(a.out)},
           "config": {"n_x": nx, "n_y": ny, "n_z": nz, "order": k},
           "versions": {"ngsolve": ngs.__version__, "numpy": np.__version__, "python": platform.python_version()},
           "input_field": "B_h = curl(A_h), A_h in HCurl(order k), copied into HDiv(order k); no H, no constitutive law",
           "projection": pdiag, "patch_flux": patch, "surfaces": surf_rec,
           "display": {"msh": os.path.basename(msh), "msh_sha256": hashlib.sha256(open(msh, "rb").read()).hexdigest(),
                       "nodes": int(sum(len(b[3]) for b in nblocks)), "elements": int(sum(len(b[3]) for b in eblocks)),
                       "gmsh_role": "display only"},
           "t_s": time.perf_counter() - t0}
    with open(a.out + ".json", "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    for r in surf_rec:
        print("%s %+.4f abs/norm tracked %.2e exactB %.2e  signed %.1e  crossings %d nJ %.1e"
              % (r["family"], r["label"], r["tracked_Bh"]["abs_over_norm"], r["oracle_exactB"]["abs_over_norm"],
                 r["tracked_Bh"]["signed_over_norm"], r["shared_face_crossings"], r["normal_jump_max_rel"]))
    print("patch flux M=%d: mid %.2e outlet %.2e; wrote %s in %.1f s" % (a.loop_m, patch["mid_rel_max"], patch["outlet_rel_max"], msh, res["t_s"]))


def write_msh42_entities(msh, groups, entities, nblocks, eblocks, quad_owner, node_dist):
    """Assemble data fields after element ids are known, then write."""
    # element ids are sequential over blocks in order; compute them up front
    start, eid = [], 1
    for b in eblocks:
        start.append(eid)
        eid += len(b[3])
    lab_items, leak_items = [], []
    for bi, label0, leak in quad_owner:
        ids = start[bi] + np.arange(len(leak))
        lab_items += [(int(i), float(label0)) for i in ids]
        leak_items += list(zip(ids.tolist(), leak.tolist()))
    return write_msh42(msh, groups, entities, nblocks, eblocks,
                       [("distance_to_exact_surface", node_dist)],
                       [("family_label", lab_items), ("local_leakage |Bh.n|/|Bh|", leak_items)])


if __name__ == "__main__":
    main()
