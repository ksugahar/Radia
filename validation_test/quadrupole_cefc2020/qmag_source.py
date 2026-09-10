"""Check that a Q-mag run uses the requested checkout or installed wheel."""

from importlib import metadata
import json
from pathlib import Path


def _same_file(left: Path, right: Path) -> bool:
    """Do two names denote the same file?

    ``Path.resolve()`` alone cannot answer this on the validation hosts: mdx1
    and hibino stage ``<checkout>/src/radia`` as a DIRECTORY JUNCTION into the
    venv's site-packages, so resolving the checkout name yields the
    site-packages name and a plain equality test rejects the very tree the
    caller asked for.  ``samefile`` compares the underlying file and therefore
    passes through the junction in either direction.  It needs both paths to
    exist, so a non-existent pair (unit tests, a staged path that was never
    created) falls back to the resolved-name comparison.
    """
    try:
        return left.samefile(right)
    except (OSError, ValueError):
        return left.resolve() == right.resolve()


def require_radia_source(source: str, module_path, checkout_src) -> dict:
    """Validate import provenance before a run can claim wheel verification."""
    if source not in ("repo", "installed"):
        raise ValueError("radia source must be 'repo' or 'installed'")
    module = Path(module_path).resolve()
    checkout = Path(checkout_src).resolve()
    named = Path(module_path).absolute()
    # A junction makes the resolved name leave the checkout, so the raw name
    # the interpreter actually imported through counts as well: claiming a
    # wheel verification while Python loaded through the checkout is exactly
    # what this gate exists to catch.
    from_checkout = (module.is_relative_to(checkout)
                     or named.is_relative_to(Path(checkout_src).absolute()))
    result = {"requested": source, "resolved_module": str(module),
              "imported_as": str(named),
              "resolved_from_checkout": from_checkout}
    if source == "repo":
        if not _same_file(module, checkout / "radia" / "__init__.py"):
            raise RuntimeError(
                f"--radia-source repo resolved to {module}, not the requested "
                f"checkout at {checkout}")
        return result

    if from_checkout:
        raise RuntimeError(
            f"--radia-source installed was imported as {named} (resolving to "
            f"{module}), which is the checkout at {checkout}; use a dedicated "
            "environment with a non-editable Radia wheel")
    try:
        dist = metadata.distribution("radia")
    except metadata.PackageNotFoundError as exc:
        raise RuntimeError("No installed Radia distribution metadata was found") from exc
    direct_url_text = dist.read_text("direct_url.json")
    if direct_url_text:
        direct_url = json.loads(direct_url_text)
        if direct_url.get("dir_info", {}).get("editable", False):
            raise RuntimeError(
                "--radia-source installed requires a non-editable Radia wheel; "
                f"the installed distribution is editable ({module})")
    if not dist.read_text("WHEEL"):
        raise RuntimeError("The Radia distribution has no WHEEL metadata")

    # A different editable tree or PYTHONPATH entry can shadow a genuine wheel.
    # Match the imported module to the distribution's recorded file, not merely
    # to a directory outside this runner's checkout.
    recorded = [Path(dist.locate_file(item)).resolve()
                for item in (dist.files or ())
                if str(item).replace("\\", "/") == "radia/__init__.py"]
    if not recorded or not all(_same_file(item, module) for item in recorded):
        raise RuntimeError(
            f"Loaded Radia module {module} does not match the installed wheel's "
            f"recorded radia/__init__.py: {recorded}")
    result.update(distribution_version=dist.version,
                  distribution_module=str(recorded[0]),
                  editable=False)
    return result
