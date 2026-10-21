"""Dependency-free distribution policy for native wheel payloads."""


class NativePayloadDistribution:
    """Tell the build backend to place MEX/Python binaries in platlib."""

    def has_ext_modules(self):
        return True
