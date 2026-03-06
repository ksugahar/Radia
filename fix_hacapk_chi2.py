#!/usr/bin/env python
"""Fix HACApK chi update to use poly->CurrentChi like LU solver."""

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

old_code = '''					// FIX (2025-12-22): Use ELF-style dual-method chi update
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

new_code = '''					// FIX (2025-12-22): Use ELF-style dual-method chi update
					// Same algorithm as LU solver (ComputeChiDualMethod)
					// This is essential for convergence matching ELF
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

						// Compute H_new from solved M using H = M / chi_current
						// This matches LU solver behavior
						double chi_current = poly->CurrentChi;
						if(chi_current < 1.0e-6) chi_current = 1.0e-6;

						TVector3d H_new;
						H_new.x = g3dRelaxPtr->Magn.x / chi_current;
						H_new.y = g3dRelaxPtr->Magn.y / chi_current;
						H_new.z = g3dRelaxPtr->Magn.z / chi_current;
						double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

						// Get mu_old from CurrentChi (same as LU solver line 1699-1700)
						double mu_old = chi_current + 1.0;

						// Try ELF-style dual-method chi update (same as LU solver)
						double chi_new;
						radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
						if(NonlinMater != nullptr)
						{
							chi_new = NonlinMater->ComputeChiDualMethod(H_mag, mu_old);
						}
						else
						{
							chi_new = chi_current;  // Fallback for non-BH curve materials
						}
						if(chi_new < 1.0e-6) chi_new = 1.0e-6;

						// Store chi for next iteration's matrix build
						poly->CurrentChi = chi_new;

						// Update NewFieldArray with H_new for convergence check
						IntrctPtr->NewFieldArray[elem] = H_new;
					}'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Successfully applied fix 2!')
else:
    print('Old code not found')
    # Check current state
    if 'chi_current = poly->CurrentChi' in content:
        print('Already fixed with poly->CurrentChi!')
    else:
        print('ERROR: Code state unknown')
