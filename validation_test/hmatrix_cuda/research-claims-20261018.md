# Preliminary research results suitable for a proposal

These are bounded research-candidate results, not released-product guarantees.
Different cases, CPU thread counts and timing boundaries must not be combined.

## Linear HDiv-MMM: complete warmed API

Source: [saved numerical and performance JSON](results/hdiv-sparse-spd-reuse-20261018.json),
introduced by commit `04277d411698571971b1b8016ca98b4d94fde451`.

| Item | Recorded result |
| --- | --- |
| Public case | P2 unit sphere, 10,344 HDiv DoFs, 4,208 charges |
| CPU baseline | Same compute host, NGSolve 16 threads, BLAS 1 thread |
| CPU warmed API median | 21.15150205 s |
| Single-GPU warmed API median | 2.06130690 s |
| CPU/GPU ratio | 10.2612 |
| Independent CPU coefficient relative difference | 4.5074e-10 |
| Original-equation true relative residual | 8.6866e-12 |
| Accounted native setup peak | 1,791,115,980 bytes; not global peak |

The two-cycle ABBA boundary includes FE setup, assembly, transfers, solve,
result recovery and destruction. Engine initialization, file exchange,
independent validation and controller work are excluded. The CPU and GPU are
different public API implementations, including their host-processing costs;
this is not an isolated kernel comparison. Calls use initialized modules,
not a complete cold process. Immutable CPU-wheel adoption remains provisional;
the measured native identities are recorded, but common-wheel release
acceptance is not established by this experiment.

Suggested Japanese wording:

> 公開球モデルを用いた10,344自由度の線形HDiv-MMMについて、単一GPUによる
> 構築・転送・求解・結果回収を含む解析API時間を約2.06秒とし、同一ホストの
> CPU 16スレッド実行に対して約10.3倍の短縮を研究用候補実装で予備検証した。
> 独立CPU解との差と元の方程式の真の残差を確認している。

Do not replace the CPU baseline with a one-thread result from another case.
Earlier 19.3x/1.68x observations may remain valid for their own declared case;
they are not the thread-scaling pair for this new 10,344-DoF result.

## Nonlinear BH: separate numerical evidence

[The 10,344-DoF nonlinear record](results/nonlinear-chunks-10344-20261018.json)
uses an independently authored synthetic monotone BH table, H=100,000 A/m,
explicitly matched all-high-order quadrature, and a fixed SPD polynomial
preconditioner. Its single complete GPU call is 12.5761 s, with two Newton
updates. Independent CPU coefficient difference is 1.4781e-9 and the
independently reconstructed nonlinear relative residual is 2.8633e-11.
This is a single-call diagnostic, not a CPU/GPU speed ratio or general
iteration-count promise. It does not certify measured steel, IH coupling,
topology optimization, arbitrary geometry or 100k-DoF scaling.

New prepared-context timing must be cited separately: it excludes initial
FE/material-map/H-matrix construction from each reused solve. Preparation and
release must be included when reporting a complete series or amortized cost.
Never relabel a prepared nonlinear solve as a two-second one-shot analysis.


## Subsequent nonlinear measurements (2026-10-19)

[chebyshev-preconditioner-20261019.json](results/chebyshev-preconditioner-20261019.json)
records the newer **10,344-DoF nonlinear complete API median of 5.12746 s**
on one A100, with two Newton updates, including FE/material preparation,
H-matrix construction, transfers, solve, recovery and destruction. The
retained-preparation nonlinear solve median is **2.08163 s**, excluding initial
preparation. Both are moderate synthetic-BH public-sphere results. The strong
saturation small case requires eleven Newton updates, so two updates is not a
general nonlinear convergence claim. This GPU/GPU comparison does not establish
a nonlinear CPU/GPU speed ratio. Native nonlinear residual is about 1e-14;
independent CPU reconstruction is below 2.864e-11 against its declared gate.
The controlled device-array upper bound remains below 2 GiB; total host,
CUDA library/driver and global peaks are not measured. These results do not
certify 100k DoF, RadFld-to-GridFunction projection at that size, measured steel,
IH coupling, topology optimization or released common-wheel distribution.
