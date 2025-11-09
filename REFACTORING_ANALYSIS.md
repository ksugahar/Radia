# Legacy File Refactoring Analysis

## Files to Refactor (27 files)

### Category 1: Header-only files (Type definitions, constants)
These are likely simple to refactor - mostly renaming:

1. **radappl.h** - Application definitions
2. **radauxst.h** - Auxiliary structures
3. **radcnvrg.h** - Convergence definitions
4. **radexcep.h** - Exception definitions
5. **radg.h** - Geometry base definitions
6. **radg3da1.h** - 3D geometry auxiliary definitions
7. **radhandl.h** - Handle definitions
8. **radmtra1.h** - Material auxiliary definitions
9. **radopnam.h** - Operation names
10. **radstlon.h** - String/long definitions
11. **radtrans.h** - Transform definitions
12. **radyield.h** - Yield definitions

**Refactoring approach**: Simple rename to `rad_xxx.h` format

---

### Category 2: Core implementation files (需要慎重)
These contain critical logic and need careful refactoring:

13. **radcast.cpp / radcast.h** - Type casting utilities
14. **radg3d.cpp / radg3d.h** - 3D geometry base class
15. **radgroup.cpp / radgroup.h** - Object grouping
16. **radinter.cpp** - Interaction calculations
17. **radmater.cpp / radmater.h** - Material definitions
18. **radmaterial.cpp** - Material implementation
19. **radsend.cpp / radsend.h** - Data serialization/send
20. **radtransform.cpp** - Transformation operations

**Refactoring approach**:
- Rename with analysis
- Add unit tests first
- Ensure no breaking changes

---

### Category 3: Recently added H-matrix files (新しいファイル)
These were recently added for H-matrix support:

21. **radhmat.cpp / radhmat.h** - H-matrix field source (NEW - 2025)

**Note**: These may already follow better practices, check if renaming is needed

---

## Proposed New Names

| Old Name | New Name | Category | Priority |
|----------|----------|----------|----------|
| radappl.h | rad_application.h | Header | Low |
| radauxst.h | rad_auxiliary_structures.h | Header | Low |
| radcast.cpp/h | rad_type_cast.cpp/h | Core | High |
| radcnvrg.h | rad_convergence.h | Header | Low |
| radexcep.h | rad_exception.h | Header | Low |
| radg.h | rad_geometry_base.h | Header | Medium |
| radg3d.cpp/h | rad_geometry_3d.cpp/h | Core | High |
| radg3da1.h | rad_geometry_3d_aux.h | Header | Low |
| radgroup.cpp/h | rad_group.cpp/h | Core | High |
| radhandl.h | rad_handle.h | Header | Low |
| radhmat.cpp/h | rad_hmatrix.cpp/h | Core | Medium |
| radinter.cpp | rad_interaction.cpp | Core | High |
| radmater.cpp/h | rad_material_def.cpp/h | Core | High |
| radmaterial.cpp | rad_material_impl.cpp | Core | High |
| radmtra1.h | rad_material_aux.h | Header | Low |
| radopnam.h | rad_operation_names.h | Header | Low |
| radsend.cpp/h | rad_serialization.cpp/h | Core | Medium |
| radstlon.h | rad_string_long.h | Header | Low |
| radtrans.h | rad_transform_def.h | Header | Low |
| radtransform.cpp | rad_transform_impl.cpp | Core | High |
| radyield.h | rad_yield.h | Header | Low |

---

## Unit Test Coverage Analysis

### Existing Tests
Check `tests/` directory for current coverage.

### Files Lacking Unit Tests (Priority for test addition)

1. **radcast.cpp** - Type casting utilities
   - Need: Test all casting functions

2. **radgroup.cpp** - Object grouping
   - Need: Test group creation, manipulation

3. **radinter.cpp** - Interaction calculations
   - Need: Test interaction matrix computation

4. **radmaterial.cpp** - Material implementation
   - Need: Test material application, retrieval

5. **radtransform.cpp** - Transformations
   - Need: Test all transform types (translation, rotation, etc.)

6. **radsend.cpp** - Serialization
   - Need: Test serialization/deserialization

---

## Refactoring Strategy

### Phase 1: Analysis (Current)
- [x] List all legacy files
- [ ] Analyze dependencies
- [ ] Categorize by complexity

### Phase 2: Add Unit Tests FIRST
- [ ] Add tests for critical functions
- [ ] Ensure 80%+ code coverage for core files
- [ ] Validate current behavior

### Phase 3: Refactor High-Priority Files
- [ ] radcast.cpp/h
- [ ] radg3d.cpp/h
- [ ] radgroup.cpp/h
- [ ] radinter.cpp
- [ ] radmater.cpp/h
- [ ] radtransform.cpp

### Phase 4: Refactor Medium-Priority Files
- [ ] radhmat.cpp/h
- [ ] radsend.cpp/h

### Phase 5: Refactor Low-Priority Files (Headers)
- [ ] All header-only files

---

## Risk Assessment

**High Risk** (需要慎重テスト):
- radg3d.cpp/h (core geometry class)
- radinter.cpp (interaction calculations)
- radmaterial.cpp (material system)

**Medium Risk**:
- radcast.cpp/h (type casting)
- radgroup.cpp/h (grouping)
- radtransform.cpp (transformations)

**Low Risk**:
- Header-only files (mostly definitions)

---

**Date**: 2025-11-09
**Status**: Analysis Phase
