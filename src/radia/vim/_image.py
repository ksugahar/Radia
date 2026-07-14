"""Mirror-image utilities for HDiv-VIM symmetry reduction."""


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
