# Development Scripts

This directory contains utility scripts used during development and maintenance of the Radia codebase.

## Scripts

### apply_security_fixes.py

**Purpose**: Apply security vulnerability fixes to the codebase

**Usage**:
```bash
python docs/scripts/apply_security_fixes.py
```

**What it does**:
1. Fixes buffer overflow in `src/python/radpy.cpp` (CombErStr function)
2. Fixes array bounds check in `src/python/pyparse.h` (CopyPyStringToC function)
3. Removes unnecessary `Py_XINCREF` calls that cause memory leaks (43 locations)

**Applied**: 2025-10-30 (commit f3e776a)

**Status**: ✅ Already applied - kept for reference only

---

### fix_test_imports.py

**Purpose**: Update import paths in test files after reorganization to tests/ directory

**Usage**:
```bash
python docs/scripts/fix_test_imports.py
```

**What it does**:
- Updates all test files in `tests/` directory
- Updates all benchmark files in `tests/benchmarks/`
- Changes import paths from `build/lib/Release` to `../build/lib/Release`

**Applied**: 2025-10-30 (commit 5670f4b)

**Status**: ✅ Already applied - kept for reference only

---

### security_fixes.patch

**Purpose**: Git patch file showing the security fixes in unified diff format

**Usage**:
```bash
# View the patch
cat docs/scripts/security_fixes.patch

# Apply the patch (if needed)
git apply docs/scripts/security_fixes.patch
```

**Contents**:
- Buffer overflow fix in radpy.cpp
- Array bounds fix in pyparse.h

**Status**: ✅ Already applied - kept for reference only

---

## Notes

These scripts were created during the security audit and code reorganization efforts in October 2025. They are preserved for:

1. **Documentation**: Understanding what changes were made
2. **Reference**: Future similar fixes
3. **Reproducibility**: Being able to re-apply fixes if needed
4. **History**: Tracking development process

**⚠️ Warning**: Do not run these scripts again on the current codebase as the fixes have already been applied. Running them again may cause errors or duplicate changes.

---

## Related Documentation

- [SECURITY_FIXES.md](../../SECURITY_FIXES.md) - Detailed security vulnerability documentation
- [tests/README.md](../../tests/README.md) - Testing documentation
- Git history: `git log --oneline -10`

---

**Last Updated**: 2025-10-30
