# Standalone Arc Kernel Check

From the repository root, with a C++17 compiler and CMake:

```powershell
cmake -S tests/cpp/arc_section -B build/arc-section
cmake --build build/arc-section --config Release
ctest --test-dir build/arc-section -C Release --output-on-failure
```

This project compiles the committed arc header directly, without Python,
NGSolve, or the rest of the native extension. It checks an axial closed form,
four independent volume-integration references with a convergence check,
five boundary partition identities, and ten two-sided boundary limits.
The partition and continuity checks are consistency checks, not independent
field references. Nonfinite errors fail explicitly, including Release builds.

This target does not certify the Python extension build, CoilBuilder placement,
or vector-potential accuracy. Those require the corresponding integration tests.
