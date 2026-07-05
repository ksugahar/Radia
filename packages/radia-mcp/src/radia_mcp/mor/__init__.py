"""radia_mcp.mor — Model Order Reduction for eddy-current FEM.

Central focus: **Cauer Ladder Network (CLN)** — a LAB SPECIALTY where
菅原研 (Sugahara Lab) is a primary contributor. Multiple foundational
papers feature K. Sugahara as co-author:
  - Kameari-Ebrahimi-Sugahara-Shindo-Matsuo 2018 (CLN foundation)
  - Sugahara-Kameari-Ebrahimi-Shindo-Matsuo 2018 (unbounded CLN)
  - Ebrahimi-Sugahara-Matsuo-Kaimori-Kameari 2020 (3D quasi-static)
  - Kuriyama-Kameari-Ebrahimi-Fujiwara-Sugahara-Shindo-Matsuo 2019
    (multiple expansion points)
  - Matsuo-Kameari-Sugahara-Shindo 2018 (matrix formulation)

Also covers general MOR for industrial inductors, motor + hybrid twins,
PRIMA, Arnoldi-Krylov, POD methods.

Distilled from:
  - All of public-safe curated corpus (~32 PDFs)
  - Industrial inductor application (Köster-König-Birò Graz 2021)
"""
