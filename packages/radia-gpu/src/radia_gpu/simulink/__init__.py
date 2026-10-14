"""Simulink boundary helpers and source adapters for Radia production models."""

from importlib import import_module

__all__ = [
    "IHOperatorAssemblyOptions",
    "assemble_ih_operators",
    "default_output_path",
]


def __getattr__(name):
    if name in __all__:
        module = import_module(".ih_operator_assembly", __name__)
        return getattr(module, name)
    raise AttributeError(name)
