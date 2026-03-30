"""
IMA (Image Method of Analysis) field evaluation helper.

NOTE: As of 2026-03-31, rad.Fld() natively supports IMA mirror contributions
for ALL element types (hex, wedge, tet) via RadIMAFieldContext. This module
is no longer needed for basic IMA field evaluation.

This module remains available for cases where explicit image elements are
needed (e.g., adding images to a different container, or for visualization).

Usage (legacy, no longer required for field evaluation):
    import radia as rad
    from radia.ima_field import add_ima_images

    # Create quarter model
    objs = [rad.ObjTetrahedron(v, [0,0,0]) for v in quarter_verts]
    container = rad.ObjCnt(objs)
    model = rad.ObjCnt([container, coil])

    # Solve with IMA
    rad.Solve(model, 0.0001, 100, 0, image='+x-z')

    # rad.Fld() now automatically includes IMA mirror contributions
    B = rad.Fld(model, 'b', [0, 0, 0])

    # Explicit image elements (only needed for visualization or separate containers):
    # add_ima_images(objs, quarter_verts, image='+x-z', container=model)
"""

import radia as rad


def parse_image_spec(image):
    """Parse image string (e.g., '+x-z') into axis flags and signs.

    Returns:
        list of (axis, sign) tuples. axis is 'x','y','z'. sign is +1 or -1.
        sign is the IMA BC sign: +1 for symmetric (field parallel to mirror),
        -1 for antisymmetric (field perpendicular to mirror).
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


def _image_magnetization(M, mx=False, my=False, mz=False, bc_sign=1):
    """Compute image magnetization from physical mirror symmetry with BC sign.

    For any mirror about axis k:
      M_k flips (pseudo-vector reflection)
      M_other stays

    The BC sign (+1 or -1) is then applied to the entire mirrored M.
    This accounts for symmetric (+) vs antisymmetric (-) boundary conditions.

    Args:
        M: [Mx, My, Mz] magnetization
        mx, my, mz: which axes to mirror
        bc_sign: combined BC sign for this mirror combination
    """
    Mx, My, Mz = M
    if mx:
        Mx = -Mx
    if my:
        My = -My
    if mz:
        Mz = -Mz
    return [Mx * bc_sign, My * bc_sign, Mz * bc_sign]


def add_ima_images(objs, vert_list, image, container=None):
    """Create IMA image elements with correct magnetization and add to container.

    After solving with IMA, call this function to create mirror copies of the
    solved elements. The images have the physically correct magnetization
    (pseudo-vector flip + BC sign). They are added to the same container
    so that rad.Fld() includes their field contributions.

    NOTE: As of 2026-03-31, rad.Fld() natively handles IMA mirrors via
    RadIMAFieldContext. This function is only needed for explicit image
    element creation (e.g., visualization or separate containers).

    Args:
        objs: list of Radia object handles (quarter model iron elements)
        vert_list: list of vertex lists, one per element (same order as objs).
            For tetrahedra: each entry is [[x,y,z], [x,y,z], [x,y,z], [x,y,z]].
            For hexahedra: each entry is [[x,y,z], ...] (8 vertices).
            For wedges: each entry is [[x,y,z], ...] (6 vertices).
        image: IMA specification string, e.g. '+x-z'
        container: Radia container to add images to. If None, images are
            added to a new container and the handle is returned.

    Returns:
        list of image element handles (for reference/cleanup)
    """
    mirrors = parse_image_spec(image)
    if not mirrors:
        return []

    # Build sign lookup: axis -> BC sign
    sign_map = {ax: s for ax, s in mirrors}

    # Build list of mirror combinations with combined BC signs
    # Each combo is (set_of_axes, combined_bc_sign)
    mirror_combos = []

    # Single axis mirrors
    for ax, s in mirrors:
        mirror_combos.append(({ax}, s))

    # Dual axis mirrors
    for i in range(len(mirrors)):
        for j in range(i + 1, len(mirrors)):
            combo_axes = {mirrors[i][0], mirrors[j][0]}
            combo_sign = mirrors[i][1] * mirrors[j][1]
            mirror_combos.append((combo_axes, combo_sign))

    # Triple axis mirror (if 3 axes)
    if len(mirrors) >= 3:
        combo_axes = {m[0] for m in mirrors}
        combo_sign = 1
        for _, s in mirrors:
            combo_sign *= s
        mirror_combos.append((combo_axes, combo_sign))

    image_objs = []
    for i, obj in enumerate(objs):
        # Get solved magnetization
        info = rad.ObjM(obj)
        M = list(info['magnetization'])
        verts = vert_list[i]
        n_verts = len(verts)

        for combo_axes, bc_sign in mirror_combos:
            mx = 'x' in combo_axes
            my = 'y' in combo_axes
            mz = 'z' in combo_axes

            mir_verts = _mirror_verts(verts, mx=mx, my=my, mz=mz)
            mir_M = _image_magnetization(M, mx=mx, my=my, mz=mz, bc_sign=bc_sign)

            if n_verts == 4:
                img = rad.ObjTetrahedron(mir_verts, mir_M)
            elif n_verts == 6:
                img = rad.ObjWedge(mir_verts, mir_M)
            elif n_verts == 8:
                img = rad.ObjHexahedron(mir_verts, mir_M)
            else:
                continue

            image_objs.append(img)

    # Add to container
    if container is not None and image_objs:
        for img in image_objs:
            rad.ObjAddToCnt(container, [img])

    return image_objs
