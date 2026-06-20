# Legacy FP64 hex VIM scripts (deprecated 2026-05-12)

These scripts are **FP64-precision** hex VIM K assembly and analysis tools, used to
generate the cuboid 5×2×1 / A1 / cylinder / sphere Cauer extraction results
documented in v2 progress reports through 2026-05-11.

**Status**: DEPRECATED. Superseded by **double-double (DD ~32 digit)** pipeline
in the parent directory (`dd_*.py` files).

## Policy (Sugahara 2026-05-12)

> radia-vim hex CUDA/GPU 実装は **DD 一本**で行く。FP64 versions は legacy
> validation reference 扱い、production research code は DD pipeline で生成。

See `memory/feedback_dd_canonical_policy.md` for the rationale.

## Migration map

| Legacy FP64 (this dir) | Canonical DD (parent dir) |
|---|---|
| `hex_vim_cupy.py` | `dd_basis.py` (HDiv basis eval) |
| `hex_vim_cupy_kassembly.py` | `dd_k_assembly_gpu.py` (Spherical Duffy K) |
| `hex_vim_gpu_order_sweep.py` | `dd_production_*.py` (production runs) |
| `hex_vim_kameari_from_foster.py` | `dd_full_pipeline.py` (Cauer extraction) |
| `hex_vim_cauer_interval.py` | Built into `dd_full_pipeline.py` verified-interval section |

## Algorithm equivalence experiments (FP64)

These scripts demonstrated the **4-method equivalence proof** (memory:
`project_cauer_stage_limit_4method_proof.md`) — that QD-Padé, Stieltjes,
Modified Chebyshev (Wheeler), and VIM-direct Lanczos+Arnoldi all give the
same Cauer rungs and are bounded by the same FP64 input-precision limit:

- `hex_vim_stieltjes_cauer.py` — Stieltjes 3-term recurrence
- `hex_vim_modified_chebyshev.py` — Wheeler's modified Chebyshev algorithm
- `hex_vim_lanczos_arnoldi.py` — VIM-direct Lanczos with Arnoldi reorthog
- `hex_vim_hiruma_arnoldi.py` — Hiruma 3-term initial attempt (unit mismatch noted)

These are **kept for research documentation** but not for production use.
DD pipeline subsumes all of them.

## When to use these legacy scripts

Only for:
1. Cross-validating new DD results against FP64 reference (FP64 precision floor)
2. Quick exploratory runs where stage 4-5 Cauer rungs are sufficient
3. Algorithm-equivalence reproductions

For production research output, **always use the DD pipeline** (`dd_*.py`).
