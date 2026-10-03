"""Build isolated serial experiments; never modify/install production SparseSolv."""
import argparse
import difflib
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import pybind11


def replace(text, old, new, count):
    if text.count(old) != count:
        raise RuntimeError(f"prototype no longer matches source: {old!r}")
    return text.replace(old, new)


def hermitian(text):
    text = replace(text, "    static bool finite_scalar(const Scalar& v) {", """    Scalar adjoint(const Scalar& v) const {
        if constexpr (std::is_same_v<Scalar, double>) return v;
        else return config_.conjugate ? std::conj(v) : v;
    }

    static bool finite_scalar(const Scalar& v) {""", 1)
    text = replace(text, "L_.values[ii] * L_.values[jj] * inv_diag_[k]",
                   "L_.values[ii] * adjoint(L_.values[jj]) * inv_diag_[k]", 2)
    text = replace(text, "L_.values[kk] * L_.values[kk] * inv_diag_[k]",
                   "L_.values[kk] * adjoint(L_.values[kk]) * inv_diag_[k]", 2)
    text = replace(text, "        Lt_ = L_.transpose();",
                   "        Lt_ = L_.transpose();\n        for (auto& v : Lt_.values) v = adjoint(v);", 1)
    return text


def natural(text):
    # Re-evaluate the earlier natural-row prototype against the current contract.
    # Only the serial path changes; parallel scheduling remains untouched.
    for direction in ("forward", "backward"):
        signature = f"    void {direction}_substitution(const Scalar* x, Scalar* y) const {{"
        if direction == "forward":
            body = """
        if (get_num_threads() == 1 && size_ >= 49152) {
            for (index_t i = 0; i < size_; ++i) {
                Scalar s = x[i];
                const index_t end = L_.row_ptr[i + 1] - 1;
                for (index_t k = L_.row_ptr[i]; k < end; ++k)
                    s -= L_.values[k] * y[L_.col_idx[k]];
                y[i] = s / L_.values[end];
            }
            return;
        }"""
        else:
            body = """
        if (get_num_threads() == 1 && size_ >= 49152) {
            for (index_t i = size_; i-- > 0;) {
                Scalar s = Scalar(0);
                for (index_t k = Lt_.row_ptr[i] + 1; k < Lt_.row_ptr[i + 1]; ++k)
                    s -= Lt_.values[k] * y[Lt_.col_idx[k]];
                y[i] = s * inv_diag_[i] + x[i];
            }
            return;
        }"""
        text = replace(text, signature, signature + body, 1)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    repo = here.parents[2]
    source = repo / "src/ext/sparsesolv/include"
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    hashes = {}
    relative = Path("sparsesolv/preconditioners/ic_preconditioner.hpp")
    original = (source / relative).read_text(encoding="utf-8")
    for variant in ("baseline", "natural", "hermitian", "combined"):
        root = out / variant
        shutil.copytree(source, root)
        text = original
        if variant in ("hermitian", "combined"):
            text = hermitian(text)
        if variant in ("natural", "combined"):
            text = natural(text)
        (root / relative).write_text(text, encoding="utf-8", newline="\n")
        (out / f"{variant}.patch").write_text("".join(difflib.unified_diff(
            original.splitlines(True), text.splitlines(True),
            fromfile="a/" + relative.as_posix(), tofile="b/" + relative.as_posix()
        )), encoding="utf-8")
        hashes[variant] = hashlib.sha256(text.encode()).hexdigest()
    (out / "variants.json").write_text(json.dumps(hashes, indent=2), encoding="utf-8")
    pybind = pybind11.get_cmake_dir()
    subprocess.run(["cmake", "-S", str(here), "-B", str(out / "build"),
                    "-A", "x64", f"-DVARIANT_ROOT={out.as_posix()}",
                    f"-Dpybind11_DIR={pybind}", f"-DPython_EXECUTABLE={sys.executable}"], check=True)
    subprocess.run(["cmake", "--build", str(out / "build"), "--config", "Release", "--parallel", "2"], check=True)


if __name__ == "__main__":
    main()
