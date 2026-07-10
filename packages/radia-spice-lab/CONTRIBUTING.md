# Contributing to ltspice-converter

Thanks for your interest. This project is intentionally narrow in
scope (file-format conversion only), so contributions are most welcome
in four shapes:

1. **Failing-circuit reports**
2. **Pull requests that improve pass rate on a documented test corpus**
3. **New supported element classes** with both a round-trip test and a
   benchmark delta
4. **Documentation** (clarity, examples, edge cases)

## 1. Reporting a failing circuit

The single most useful contribution is a minimal `.asc` or `.cir` file
the converter mishandles, posted as a [GitHub issue](
https://github.com/ksugahar/ltspice-converter/issues/new/choose).

Please redact anything you cannot share publicly (vendor IP, NDA
material). The smaller the failing example, the faster the fix.

## 2. Dev setup

```bash
git clone https://github.com/ksugahar/ltspice-converter
cd ltspice-converter
python -m venv .venv
. .venv/Scripts/activate          # or `source .venv/bin/activate` on Linux/Mac
pip install -e ".[mcp,test]"
```

This project is distributed via GitHub only (no PyPI release). End
users install with:

```bash
pip install git+https://github.com/ksugahar/ltspice-converter
```

Verify the install:

```bash
python -m pytest tests/ -q
# Should report "32 passed" (31 round-trip + 1 vendor-symbol regression)
```

## 3. Running the public test suite

```bash
python -m pytest tests/ -v
```

All tests use author-authored minimal circuits committed under
`tests/fixtures/`. No external corpora are required for the public
suite.

## 4. Running the private benchmark

`bench/baseline.py` (gitignored) runs the converter over a
LAB-private corpus of 5,000+ real-world `.asc` files. The directory
structure expected by the script is in `bench/baseline.py:CORPUS_ROOTS`.
You can adapt the script to point at your own corpus:

```bash
python bench/baseline.py --per-source 100 --out bench/my_baseline.json
```

The corpus itself MUST NOT be added to the repository; the converter
is the deliverable, not the training data. See [BENCHMARKS.md](
docs/BENCHMARKS.md) for methodology.

## 5. Pull request checklist

Before opening a PR:

- [ ] `python -m pytest tests/ -q` passes locally
- [ ] If the change targets a class of circuits, a new regression test
      is added to `tests/test_conversions.py` (author-authored minimal
      example, no third-party content)
- [ ] If pass-rate is impacted, run `bench/baseline.py` before and
      after and quote the delta in the PR description
- [ ] No vendor-bundled / textbook / GitHub-aggregate content is
      added to the repo (`.gitignore` should catch this; double-check
      with `git diff --cached --name-only`)
- [ ] Commit message explains the **WHY**, not just the WHAT

## 6. Style

- Python 3.10+
- Stdlib + `schemdraw` + `numpy`. No new heavy dependencies without
  prior discussion.
- ASCII in source code (no Unicode `→ ≤ ²`, write `->`, `<=`, `^2`).
  This codebase is exercised on Japanese-locale Windows consoles
  (cp932) where stray Unicode in `print()` causes immediate failures.
- Type hints encouraged.

## 7. Scope reminders

What this project does:
- `.asc` <-> `.cir` <-> schemdraw Python script round-trip

What this project does NOT do (and won't accept PRs for):
- Circuit simulation (use LTspice, ngspice, PySpice instead)
- Circuit analysis / symbolic solving (use sympy / lcapy)
- Schematic editing GUI
- LTspice symbol library redistribution
- **Textbook-derived test fixtures** -- see below

A converter PR that adds a circuit-analysis feature will be redirected
to a separate project.

### Why textbook-derived fixtures are out of scope

The converter is developed against a **private** training corpus that
includes circuits transcribed from copyrighted textbooks.  Those
circuits live ONLY on the author's LAB-private storage and are NOT
mirrored here.  PRs that add textbook-sourced `.cir` / `.asc` files
to `tests/fixtures/` (or anywhere else under this repo) will be
declined regardless of how cleanly they round-trip.  The policy is
simple: **trained converter is public, training data is private.**

For new test fixtures, author your own minimal circuit from scratch
(an RC/RLC/filter under a dozen components is plenty for exercising
any single feature) and submit that.

## 8. License

By contributing, you agree that your contributions will be licensed
under [MIT](LICENSE).
