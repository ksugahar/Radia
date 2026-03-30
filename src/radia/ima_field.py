"""
IMA (Image Method of Analysis) field evaluation helper.

After solving with IMA symmetry (e.g., image='+x-z'), the interaction matrix
correctly accounts for image elements during the solve. However, rad.Fld()
only evaluates field from the quarter model elements. This module provides
a function to create image elements with correct M for field evaluation.

Usage:
    import radia as rad
    from radia.ima_field import add_ima_images

    # Create quarter model
    objs = [rad.ObjTetrahedron(v, [0,0,0]) for v in quarter_verts]
    container = rad.ObjCnt(objs)
    model = rad.ObjCnt([container, coil])

    # Solve with IMA
    rad.Solve(model, 0.0001, 100, 0, image='+x-z')

    # Add image elements for field evaluation
    add_ima_images(objs, quarter_verts, image='+x-z')

    # Now rad.Fld includes image contributions
    B = rad.Fld(model, 'b', [0, 0, 0])
"""

import radia as rad


def parse_image_spec(image):
    """Parse image string (e.g., '+x-z') into axis flags and signs.

    Returns:
        list of (axis, sign) tuples. axis is 'x','y','z'. sign is +1 or -1.
    """
    mirrors = []
    i = 0
    while i < len(image):
        ch = image[i]
        if ch in '+-':
            sign = +1 if ch == '+' else -1
            i += 1
            if i < len(image) and image[i] in 'xyzXYZ':
                mirrors.append((image[i].lower(), sign))
                i += 1
        elif ch in 'xyzXYZ':
            mirrors.append((ch.lower(), +1))
            i += 1
        else:
            i += 1
    return mirrors


def _mirror_verts(verts, mx=False, my=False, mz=False):
    """Mirror vertex list about specified axes."""
    result = []
    for v in verts:
        x, y, z = v[0], v[1], v[2]
        if mx:
            x = -x
        if my:
            y = -y
        if mz:
            z = -z
        result.append([x, y, z])
    return result


def _image_magnetization(M, mx=False, my=False, mz=False):
    """Compute image magnetization from physical mirror symmetry.

    For any mirror about axis k:
      M_k flips (antisymmetric across mirror)
      M_other stays (symmetric across mirror)
    """
    Mx, My, Mz = M
    if mx:
        Mx = -Mx
    if my:
        My = -My
    if mz:
        Mz = -Mz
    return [Mx, My, Mz]


def add_ima_images(objs, vert_list, image, container=None):
    """Create IMA image elements with correct magnetization and add to container.

    After solving with IMA, call this function to create mirror copies of the
    solved elements. The images have the physically correct magnetization
    (component along mirror axis flips). They are added to the same container
    so that rad.Fld() includes their field contributions.

    Args:
        objs: list of Radia object handles (quarter model iron elements)
        vert_list: list of vertex lists, one per element (same order as objs).
            For tetrahedra: each entry is [[x,y,z], [x,y,z], [x,y,z], [x,y,z]].
            For hexahedra: each entry is [[x,y,z], ...] (8 vertices).
        image: IMA specification string, e.g. '+x-z'
        container: Radia container to add images to. If None, images are
            added to a new container and the handle is returned.

    Returns:
        list of image element handles (for reference/cleanup)
    """
    mirrors = parse_image_spec(image)
    if not mirrors:
        return []

    # Build list of mirror combinations (single, dual, triple)
    axes = [m[0] for m in mirrors]
    mirror_combos = []

    # Single axis mirrors
    for ax, _sign in mirrors:
        mirror_combos.append({ax})

    # Dual axis mirrors
    for i in range(len(mirrors)):
        for j in range(i + 1, len(mirrors)):
            mirror_combos.append({mirrors[i][0], mirrors[j][0]})

    # Triple axis mirror (if 3 axes)
    if len(mirrors) >= 3:
        mirror_combos.append({m[0] for m in mirrors})

    image_objs = []
    for i, obj in enumerate(objs):
        # Get solved magnetization
        info = rad.ObjM(obj)
        M = list(info['magnetization'])
        verts = vert_list[i]
        n_verts = len(verts)

        for combo in mirror_combos:
            mx = 'x' in combo
            my = 'y' in combo
            mz = 'z' in combo

            mir_verts = _mirror_verts(verts, mx=mx, my=my, mz=mz)
            mir_M = _image_magnetization(M, mx=mx, my=my, mz=mz)

            if n_verts == 4:
                img = rad.ObjTetrahedron(mir_verts, mir_M)
            elif n_verts == 8:
                img = rad.ObjHexahedron(mir_verts, mir_M)
            else:
                # Skip unsupported element types
                continue

            image_objs.append(img)

    # Add to container
    if container is not None and image_objs:
        for img in image_objs:
            rad.ObjAddToCnt(container, [img])

    return image_objs
