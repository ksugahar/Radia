# PEEC: dense LU vs HACApK BiCGSTAB scaling

Helical solenoid, 20 turns x 25 seg/turn, nwinc^2 sub-filaments. Freq 10 kHz. BiCGSTAB rtol 1e-6 with Jacobi preconditioner.

| N_fil | Dense L | Dense build | Dense solve | Dense peak mem | H-mat | Compression | H-mat build | H-mat solve | Iters | Residual | H-mat peak mem |
|------:|--------:|------------:|------------:|---------------:|------:|------------:|------------:|------------:|------:|---------:|---------------:|
| 500 | 2 MB | 20.2 ms | 67.5 ms | 144 MB | 2 MB | 100.0% | 39.9 ms | 3.1 ms | 13 | 9.2e-07 | 110 MB |
| 2000 | 31 MB | 315.5 ms | 1.65 s | 395 MB | 27 MB | 90.0% | 575.8 ms | 235.4 ms | 36 | 9.9e-07 | 137 MB |
| 4500 | 154 MB | 1.60 s | 16.13 s | 1.49 GB | 77 MB | 50.0% | 1.61 s | 842.8 ms | 43 | 5.7e-07 | 188 MB |
| 8000 | 488 MB | 5.12 s | 1.40 min | 4.42 GB | 64 MB | 13.1% | 1.27 s | 1.09 s | 67 | 8.6e-07 | 177 MB |
| 12500 | 1.16 GB | 14.89 s | 5.81 min | 10.60 GB | 80 MB | 6.7% | 1.46 s | 1.37 s | 61 | 8.8e-07 | 195 MB |
| 18000 | -- | -- | -- | -- | 90 MB | 3.6% | 1.61 s | 1.77 s | 71 | 9.9e-07 | 209 MB |
| 24500 | -- | -- | -- | -- | 157 MB | 3.4% | 2.68 s | 3.05 s | 70 | 9.3e-07 | 281 MB |

## Speedup

| N_fil | Build speedup | Solve speedup | H-mat / dense L memory |
|------:|--------------:|--------------:|-----------------------:|
| 500 | 0.5x | 21.9x | 100.0% |
| 2000 | 0.5x | 7.0x | 90.0% |
| 4500 | 1.0x | 19.1x | 50.0% |
| 8000 | 4.0x | 77.0x | 13.1% |
| 12500 | 10.2x | 254.7x | 6.7% |
| 18000 | -- | -- | -- |
| 24500 | -- | -- | -- |
