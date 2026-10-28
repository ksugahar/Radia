"""Declare the prebuilt CPU native payload without compiling extensions."""
from pathlib import Path
import shutil
from setuptools import Distribution, setup
from setuptools.command.build_py import build_py


class BinaryDistribution(Distribution):
    """The packaged pybind11 extension requires a platform/ABI wheel."""

    def has_ext_modules(self):
        return True


class BuildPython(build_py):
    """Ship the single shared GPU source closure with the Python distribution."""

    def run(self):
        super().run()
        source = Path(__file__).resolve().parent / "src/core/gpu"
        destination = Path(self.build_lib) / "radia/_gpu_sources"
        destination.mkdir(parents=True, exist_ok=True)
        for path in source.iterdir():
            if path.suffix in {".hpp", ".cuh", ".cpp", ".cu"}:
                shutil.copyfile(path, destination / path.name)


setup(distclass=BinaryDistribution, cmdclass={"build_py": BuildPython})
