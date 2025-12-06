# Hybrid Variable DOF Solver Design

## Overview

This document describes the design for extending `radTInteraction` and the LU/BiCGSTAB solvers to support **variable degrees of freedom (DOF)** per element, enabling **hybrid analysis** with:

- **MSC hexahedra (6 DOF)**: Surface charge density on each face (sigma_1...sigma_6)
- **Standard polyhedra (3 DOF)**: Magnetization components (Mx, My, Mz)

**Reference**: H. Yano, K. Sugahara, "Magnetic Moment Method with the Idea of Magnetic Surface Charge Method", J. Magn. Soc. Jpn., 2023

---

## Current Architecture (3 DOF Fixed)

### rad_interaction.h

```cpp
// Current: Fixed 3x3 interaction matrices
std::vector<std::vector<TMatrix3df>> vInteractMatrix;  // N x N array of 3x3 matrices
TMatrix3df** InteractMatrix;

// Current: Fixed TVector3d arrays (3 components)
std::vector<TVector3d> vExternFieldArray;
std::vector<TVector3d> vNewMagnArray;
std::vector<TVector3d> vNewFieldArray;
TVector3d* ExternFieldArray;
TVector3d* NewMagnArray;
TVector3d* NewFieldArray;
```

### rad_relaxation_methods.cpp

```cpp
// Current: Assumes 3 DOF per element
int ndof = 3 * AmOfMainElem;

// Current: Builds dense matrix from 3x3 blocks
for(int i = 0; i < AmOfMainElem; i++) {
    for(int comp_i = 0; comp_i < 3; comp_i++) {
        int row = 3 * i + comp_i;
        // ...
    }
}
```

---

## Proposed Architecture (Variable DOF)

### Key Design Decisions

1. **DOF Per Element**: Each element reports its DOF via `NumberOfDegOfFreedom()` virtual method
   - Standard (radTg3dRelax): returns 3
   - MSC hexahedra (radTExtrPolygonMSC): returns 6

2. **Flattened Storage**: Replace fixed-size structures with variable-size arrays
   - Total DOF = sum of all element DOFs
   - Each element has a starting offset in the flattened array

3. **Block-Based Interaction Matrix**: Store as dense row-major array with variable block sizes
   - Block (i,j) size = DOF[i] x DOF[j]

### New Data Structures (rad_interaction.h)

```cpp
class radTInteraction : public radTg {
private:
    // Variable DOF support
    int m_totalDOF;                           // Total degrees of freedom (sum of all element DOFs)
    std::vector<int> m_elemDOF;               // DOF for each element (3 or 6)
    std::vector<int> m_elemDOFOffset;         // Starting offset in flattened arrays

    // Variable-size interaction matrix (row-major, dense)
    // Size: m_totalDOF x m_totalDOF
    std::vector<double> m_flatInteractMatrix;

    // Variable-size field arrays
    std::vector<double> m_flatExternFieldArray;  // Size: m_totalDOF
    std::vector<double> m_flatMagnArray;         // Size: m_totalDOF (M or sigma values)
    std::vector<double> m_flatFieldArray;        // Size: m_totalDOF

    // Legacy compatibility pointers (may keep for backward compatibility)
    TMatrix3df** InteractMatrix;  // Points to 3x3 view for standard elements only
    TVector3d* ExternFieldArray;  // Points to TVector3d view for standard elements only

public:
    // DOF management
    int GetTotalDOF() const { return m_totalDOF; }
    int GetElementDOF(int elemIdx) const { return m_elemDOF[elemIdx]; }
    int GetElementDOFOffset(int elemIdx) const { return m_elemDOFOffset[elemIdx]; }

    // Matrix access (block view)
    double* GetInteractBlock(int row_elem, int col_elem);
    const double* GetInteractBlock(int row_elem, int col_elem) const;

    // Array access
    double* GetMagnArray() { return m_flatMagnArray.data(); }
    double* GetFieldArray() { return m_flatFieldArray.data(); }
    double* GetExternFieldArray() { return m_flatExternFieldArray.data(); }

    // Setup methods
    void ComputeDOFOffsets();
    int SetupInteractMatrix_VariableDOF();
};
```

### DOF Offset Computation

```cpp
void radTInteraction::ComputeDOFOffsets()
{
    m_elemDOF.resize(AmOfMainElem);
    m_elemDOFOffset.resize(AmOfMainElem + 1);

    m_totalDOF = 0;
    for(int i = 0; i < AmOfMainElem; i++) {
        m_elemDOFOffset[i] = m_totalDOF;
        m_elemDOF[i] = g3dRelaxPtrVect[i]->NumberOfDegOfFreedom();
        m_totalDOF += m_elemDOF[i];
    }
    m_elemDOFOffset[AmOfMainElem] = m_totalDOF;  // End sentinel
}
```

### Interaction Matrix Setup (Variable DOF)

The interaction matrix element N[i][j] computes the field at element i due to element j.

For **hybrid analysis**, we need to handle four cases:

| Element i | Element j | Block Size | Field Computation |
|-----------|-----------|------------|-------------------|
| Standard (3 DOF) | Standard (3 DOF) | 3x3 | B_comp (existing) |
| Standard (3 DOF) | MSC (6 DOF) | 3x6 | Field from sigma faces |
| MSC (6 DOF) | Standard (3 DOF) | 6x3 | B at 6 eval points |
| MSC (6 DOF) | MSC (6 DOF) | 6x6 | B at 6 eval points from sigma |

```cpp
int radTInteraction::SetupInteractMatrix_VariableDOF()
{
    ComputeDOFOffsets();

    // Allocate flattened interaction matrix
    m_flatInteractMatrix.resize(m_totalDOF * m_totalDOF, 0.0);

    for(int col = 0; col < AmOfMainElem; col++) {
        radTg3dRelax* elem_col = g3dRelaxPtrVect[col];
        int dof_col = m_elemDOF[col];
        int offset_col = m_elemDOFOffset[col];

        for(int row = 0; row < AmOfMainElem; row++) {
            radTg3dRelax* elem_row = g3dRelaxPtrVect[row];
            int dof_row = m_elemDOF[row];
            int offset_row = m_elemDOFOffset[row];

            // Compute interaction block
            // Block is stored at: m_flatInteractMatrix[offset_row * m_totalDOF + offset_col]
            // Block size: dof_row x dof_col

            ComputeInteractBlock(elem_row, elem_col, row, col,
                                &m_flatInteractMatrix[offset_row * m_totalDOF + offset_col],
                                m_totalDOF);  // Leading dimension
        }
    }

    return 1;
}
```

### Interaction Block Computation

```cpp
void radTInteraction::ComputeInteractBlock(
    radTg3dRelax* elem_row, radTg3dRelax* elem_col,
    int row_idx, int col_idx,
    double* block, int ldim)
{
    int dof_row = m_elemDOF[row_idx];
    int dof_col = m_elemDOF[col_idx];

    bool row_is_msc = (dof_row == 6);
    bool col_is_msc = (dof_col == 6);

    if(!row_is_msc && !col_is_msc) {
        // Standard 3x3: existing B_comp method
        ComputeInteractBlock_3x3(elem_row, elem_col, block, ldim);
    }
    else if(!row_is_msc && col_is_msc) {
        // 3x6: Field at row center from col's 6 sigma faces
        ComputeInteractBlock_3x6(elem_row, elem_col, block, ldim);
    }
    else if(row_is_msc && !col_is_msc) {
        // 6x3: Field at row's 6 eval points from col's magnetization
        ComputeInteractBlock_6x3(elem_row, elem_col, block, ldim);
    }
    else {
        // 6x6: Field at row's 6 eval points from col's 6 sigma faces
        ComputeInteractBlock_6x6(elem_row, elem_col, block, ldim);
    }
}
```

---

## LU Solver Update

### Current LU System (3 DOF)

```
(1/chi - N) * M = H_ext
```

Where:
- M = [Mx_1, My_1, Mz_1, Mx_2, My_2, Mz_2, ..., Mx_N, My_N, Mz_N]^T (3N components)
- N = Interaction matrix (3N x 3N)
- chi = Susceptibility (diagonal blocks)

### Updated LU System (Variable DOF)

For hybrid MSC + standard:

```
(D - N) * X = F
```

Where:
- X = [M_1, M_2, ..., sigma_j1, sigma_j2, ..., sigma_j6, ...]^T (total DOF components)
- N = Interaction matrix (total DOF x total DOF)
- D = Diagonal blocks with 1/chi (or material-dependent terms for MSC)
- F = External field contribution

### Material Handling for MSC Elements

For MSC elements, the constitutive relation differs from standard MMM:

**Standard (3 DOF)**:
```
M = chi(H) * H
```

**MSC (6 DOF)**:
```
sigma_n = M dot n_n = chi(H_n) * H_n dot n_n
```

Where H_n is the field at the evaluation point for face n.

---

## BiCGSTAB Solver Update

The BiCGSTAB matvec operation needs updating:

```cpp
void radTRelaxationMethNo_1::MatVec_VariableDOF(
    const std::vector<double>& x,
    std::vector<double>& y,
    const std::vector<double>& inv_chi,
    int totalDOF)
{
    // y = (1/chi - N) * x
    // Using variable DOF blocks

    std::fill(y.begin(), y.end(), 0.0);

    #pragma omp parallel for if(IntrctPtr->AmOfMainElem > 50)
    for(int row_elem = 0; row_elem < IntrctPtr->AmOfMainElem; row_elem++) {
        int dof_row = IntrctPtr->GetElementDOF(row_elem);
        int offset_row = IntrctPtr->GetElementDOFOffset(row_elem);

        // Diagonal contribution: (1/chi) * x_row
        for(int k = 0; k < dof_row; k++) {
            y[offset_row + k] = inv_chi[offset_row + k] * x[offset_row + k];
        }

        // Off-diagonal: -N * x
        for(int col_elem = 0; col_elem < IntrctPtr->AmOfMainElem; col_elem++) {
            int dof_col = IntrctPtr->GetElementDOF(col_elem);
            int offset_col = IntrctPtr->GetElementDOFOffset(col_elem);

            const double* block = IntrctPtr->GetInteractBlock(row_elem, col_elem);

            // y_row -= N_block * x_col
            for(int i = 0; i < dof_row; i++) {
                for(int j = 0; j < dof_col; j++) {
                    y[offset_row + i] -= block[i * totalDOF + j] * x[offset_col + j];
                }
            }
        }
    }
}
```

---

## HACApK Compatibility (Future)

The design is compatible with future HACApK acceleration:

1. **H-matrix Clustering**: Works with variable DOF since clusters are based on geometry, not DOF
2. **ACA+ Approximation**: Low-rank blocks are independent of DOF structure
3. **rad.Fld() Acceleration**: ACA+ can accelerate field evaluation independently of solver

---

## Implementation Plan

### Phase 1: Infrastructure
1. Add DOF management to `radTInteraction`
2. Add flattened storage arrays
3. Add DOF offset computation
4. Update memory allocation

### Phase 2: Interaction Matrix
1. Implement `SetupInteractMatrix_VariableDOF()`
2. Implement block computation methods (3x3, 3x6, 6x3, 6x6)
3. Keep backward compatibility for pure 3-DOF systems

### Phase 3: LU Solver
1. Update `radTRelaxationMethNo_0::AutoRelax()` for variable DOF
2. Update system matrix construction
3. Update solution extraction

### Phase 4: BiCGSTAB Solver
1. Update `radTRelaxationMethNo_1::AutoRelax()` for variable DOF
2. Update `MatVec()` for variable DOF
3. Update Jacobi preconditioner for variable DOF

### Phase 5: Testing
1. Unit tests for pure MSC (6 DOF only)
2. Unit tests for hybrid (3 DOF + 6 DOF)
3. Validation against ELF_MAGIC reference

---

## API Changes

No new Python API required. Hybrid analysis is automatic based on element types:

```python
import radia as rad

# Create hexahedra (automatically uses MSC 6 DOF)
hex1 = rad.ObjThckPgn(...)  # Will use radTExtrPolygonMSC internally

# Create polyhedra (uses standard 3 DOF)
poly1 = rad.ObjPolyhdr(...)  # Uses radTPolyhedron

# Hybrid solve works automatically
grp = rad.ObjCnt([hex1, poly1])
rad.Solve(grp, 0.0001, 1000)  # Hybrid solver detects mixed DOF
```

---

**Last Updated**: 2025-12-06
**Author**: Claude Code
**Project**: Radia MSC Hybrid Solver
