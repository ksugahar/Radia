#!/usr/bin/env python
"""Fix HACApK chi update to use ComputeChiDualMethod."""

import re

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

old_code = '''				// FIX (2025-12-21): Update NewFieldArray with H = M / chi
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

new_code = '''				// FIX (2025-12-22): Use ELF-style dual-method chi update
					// Same algorithm as LU solver (ComputeChiDualMethod)
					// This is essential for convergence matching ELF
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
						TVector3d H_est = IntrctPtr->NewFieldArray[elem];
						double H_mag = std::sqrt(H_est.x*H_est.x + H_est.y*H_est.y + H_est.z*H_est.z);

						// Get current chi for dual-method comparison
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);
						double chi_old = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi_old < 1.0e-6) chi_old = 1.0e-6;
						double mu_old = chi_old + 1.0;

						// Try ELF-style dual-method chi update (same as LU solver)
						double chi_new;
						radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
						if(NonlinMater != nullptr)
						{
							chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
						}
						else
						{
							chi_new = chi_old;  // Fallback for non-BH curve materials
						}
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;

						// Store chi for next iteration's matrix build
						poly->CurrentChi = chi_new;

						// H = M / chi (constitutive relation)
						IntrctPtr->NewFieldArray[elem].x = g3dRelaxPtr->Magn.x / chi_new;
						IntrctPtr->NewFieldArray[elem].y = g3dRelaxPtr->Magn.y / chi_new;
						IntrctPtr->NewFieldArray[elem].z = g3dRelaxPtr->Magn.z / chi_new;
					}'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Successfully applied fix!')
else:
    print('Old code not found - checking if already fixed...')
    if 'ComputeChiDualMethod(H_mag, mu_old)' in content and 'FIX (2025-12-22)' in content:
        print('Already fixed!')
    else:
        print('ERROR: Could not find code to replace')
