"""Mirror-image utilities for HDiv-VIM symmetry reduction."""
import math


def parse_image_string(image):
    """Parse tokens such as ``+x-z`` into ``(axis, sign)`` pairs."""
    value = image.strip().lower().replace(" ", "")
    planes = []
    axis_of = {"x": 0, "y": 1, "z": 2}
    index = 0
    while index < len(value):
        if (value[index] not in "+-" or index + 1 >= len(value)
                or value[index + 1] not in axis_of):
            raise ValueError(
                "bad IMA image string %r (expected tokens like '+x','-z')" % image
            )
        planes.append((axis_of[value[index + 1]], 1 if value[index] == "+" else -1))
        index += 2
    axes = [axis for axis, _ in planes]
    if len(set(axes)) != len(axes):
        raise ValueError("IMA image string %r repeats an axis" % image)
    return planes


def image_group(planes):
    """Return every non-empty reflection subset and its product sign."""
    images = []
    for mask in range(1, 1 << len(planes)):
        axes = []
        sign = 1
        for index, (axis, plane_sign) in enumerate(planes):
            if mask & (1 << index):
                axes.append(axis)
                sign *= plane_sign
        images.append((tuple(sorted(axes)), sign))
    return images


def cyclic_image_terms(planes, n_fold, alternating=False):
    """Return the (masks, signs, angles) of an N-fold rotation group, with an optional ``-z`` midplane.

    Each term is T = M o R: rotate about +z, then mirror on the mask axes. The rotations
    come first, because the connected sector pairs its seam with the first angle. Only the
    midplane commutes with every rotation about +z, and only its antisymmetric ``-z`` form
    leaves the HDiv normal trace free, which is what a field crossing the midplane needs.
    """
    planes = list(planes)
    if planes and planes != [(2, -1)]:
        raise NotImplementedError(
            "image_cyclic combines only with the '-z' midplane mirror; got planes %r" % (planes,))
    masks, signs, angles = [], [], []
    for k in range(1, n_fold):
        masks.append(0)
        signs.append(-1.0 if alternating and k % 2 else 1.0)
        angles.append(2.0*math.pi*k/n_fold)
    if planes:
        masks.append(4)
        signs.append(-1.0)
        angles.append(0.0)
        for k in range(1, n_fold):
            masks.append(4)
            signs.append(1.0 if alternating and k % 2 else -1.0)
            angles.append(2.0*math.pi*k/n_fold)
    return masks, signs, angles
