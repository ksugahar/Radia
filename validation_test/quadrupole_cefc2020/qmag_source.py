"""Check that a Q-mag run uses the requested checkout or installed wheel."""

from importlib import metadata
import json
from pathlib import Path


def require_radia_source(source: str, module_path, checkout_src) -> dict:
    """Validate import provenance before a run can claim wheel verification."""
    if source not in ("repo", "installed"):
        raise ValueError("radia source must be 'repo' or 'installed'")
    module = Path(module_path).resolve()
    checkout = Path(checkout_src).resolve()
    from_checkout = module.is_relative_to(checkout)
    result = {"requested": source, "resolved_module": str(module),
              "resolved_from_checkout": from_checkout}
    if source == "repo":
        if module != checkout / "radia" / "__init__.py":
            raise RuntimeError(
                f"--radia-source repo resolved to {module}, not the requested "
                f"checkout at {checkout}")
        return result

    if from_checkout:
        raise RuntimeError(
            f"--radia-source installed resolved to the checkout at {module}; "
            "use a dedicated environment with a non-editable Radia wheel")
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
    if recorded != [module]:
        raise RuntimeError(
            f"Loaded Radia module {module} does not match the installed wheel's "
            f"recorded radia/__init__.py: {recorded}")
    result.update(distribution_version=dist.version,
                  distribution_module=str(recorded[0]),
                  editable=False)
    return result
