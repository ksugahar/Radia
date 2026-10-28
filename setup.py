"""Declare the prebuilt CPU native payload without compiling extensions."""
from setuptools import Distribution, setup


class BinaryDistribution(Distribution):
    """The packaged pybind11 extension requires a platform/ABI wheel."""

    def has_ext_modules(self):
        return True


setup(distclass=BinaryDistribution)
