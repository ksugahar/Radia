# Radia Python package
# This module re-exports all symbols from the C++ extension module (radia.pyd)
# so that 'import radia' works correctly when installed via pip

__version__ = "1.3.12"

# Import all symbols from the C++ extension module
try:
    from radia.radia import *
except ImportError:
    # Fallback for development: try importing from the same directory
    try:
        from .radia import *
    except ImportError as e:
        raise ImportError(
            "Failed to import radia C++ extension module (radia.pyd). "
            "Ensure the package was built correctly with Build.ps1 before installation."
        ) from e
