#!/usr/bin/env python
"""Restore original 2025-12-21 version (17 iterations)."""

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Current broken state
old_code = '''					// Update chi from M using simple approach
					// Note: Using DefineInstantKsiTensor for now
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

						// Get M magnitude for H estimation
						double M_mag = std::sqrt(g3dRelaxPtr->Magn.x*g3dRelaxPtr->Magn.x +
						                         g3dRelaxPtr->Magn.y*g3dRelaxPtr->Magn.y +
						                         g3dRelaxPtr->Magn.z*g3dRelaxPtr->Magn.z);

						// Estimate H from M using current chi
						double chi_current = poly->CurrentChi;
						if(chi_current < 1.0e-6) chi_current = 1.0e-6;

						TVector3d H_est;
						H_est.x = g3dRelaxPtr->Magn.x / chi_current;
						H_est.y = g3dRelaxPtr->Magn.y / chi_current;
						H_est.z = g3dRelaxPtr->Magn.z / chi_current;

						// Get new chi from B-H curve
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
						double chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;

						// Store chi for next iteration's matrix build
						poly->CurrentChi = chi_new;

						// H = M / chi (constitutive relation)
						IntrctPtr->NewFieldArray[elem].x = g3dRelaxPtr->Magn.x / chi_new;
						IntrctPtr->NewFieldArray[elem].y = g3dRelaxPtr->Magn.y / chi_new;
						IntrctPtr->NewFieldArray[elem].z = g3dRelaxPtr->Magn.z / chi_new;
					}'''

# Restore original
new_code = '''					// FIX (2025-12-21): Update NewFieldArray with H = M / chi
					// This matches ELF's approach: H_int = M / chi (constitutive relation)
					// The updated H is used in the next iteration for chi(H) computation.
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
						TVector3d H_est = IntrctPtr->NewFieldArray[elem];
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
						double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi < 1.0e-6) chi = 1.0e-6;

						// H = M / chi (constitutive relation)
						IntrctPtr->NewFieldArray[elem].x = g3dRelaxPtr->Magn.x / chi;
						IntrctPtr->NewFieldArray[elem].y = g3dRelaxPtr->Magn.y / chi;
						IntrctPtr->NewFieldArray[elem].z = g3dRelaxPtr->Magn.z / chi;
					}'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Restored original 2025-12-21 version')
else:
    print('Old code not found - checking current state')
