#!/usr/bin/env python
"""Fix HACApK chi update to exactly match LU solver."""

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Current simple version
old_code = '''					// FIX (2025-12-21): Update NewFieldArray with H = M / chi
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

# Match LU solver exactly (lines 1693-1739)
new_code = '''					// FIX (2025-12-22): Match LU solver chi update exactly
					// See LU solver lines 1693-1739 for reference
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

						// Step 1: Compute H_new = M / chi_current (same as LU line 1695-1700)
						// chi_current is the chi used in this iteration's matrix
						double chi_current = poly->CurrentChi;
						if(chi_current < 1.0e-6) chi_current = 1.0e-6;
						double mu_old = chi_current + 1.0;

						TVector3d H_new;
						H_new.x = g3dRelaxPtr->Magn.x / chi_current;
						H_new.y = g3dRelaxPtr->Magn.y / chi_current;
						H_new.z = g3dRelaxPtr->Magn.z / chi_current;
						double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

						// Step 2: Update chi using dual-method (same as LU lines 1708-1736)
						radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
						double chi_new;
						if(NonlinMater != nullptr)
						{
							chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
						}
						else
						{
							// Fallback for linear materials
							TMatrix3d KsiTensor;
							TVector3d MrVect;
							MaterPtr->DefineInstantKsiTensor(H_new, KsiTensor, MrVect);
							chi_new = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						}
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;

						// Step 3: Store chi for next iteration's matrix (same as LU line 1739)
						poly->CurrentChi = chi_new;

						// Step 4: Update NewFieldArray for convergence check
						// Use H_new (computed with chi_current), not H with chi_new
						IntrctPtr->NewFieldArray[elem] = H_new;
					}'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Applied LU-matching chi update fix')
else:
    print('Old code not found - checking current state')
    if 'Match LU solver chi update exactly' in content:
        print('Already fixed!')
    else:
        print('Unknown state')
