#!/usr/bin/env python
"""Add ELF-style chi initialization to HACApK solver."""

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

old_code = '''	// Cache polyhedron pointers for fast access (same as LU solver)
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		polyCache[elem] = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
	}

	// B-field convergence tracking (same as LU solver, ELF mucal2)
	std::vector<double> OldBnorm(AmOfMainElem, 0.0);
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;'''

new_code = '''	// Cache polyhedron pointers for fast access (same as LU solver)
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		polyCache[elem] = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
	}

	// Initialize CurrentChi with ELF-style initial value (same as LU solver lines 1494-1522)
	// This uses BH curve's 2nd point: chi = B2/(mu0*H2) - 1
	// Without this, CurrentChi starts at 1.0 causing slow convergence
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		int dof = IntrctPtr->GetElementDOF(elem);
		if(dof == 6)
		{
			radTPolyhedron* poly = polyCache[elem];
			if(poly && poly->Use6DOF_MSC)
			{
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);
				radTNonlinearIsotropMaterial* NonlinMater = dynamic_cast<radTNonlinearIsotropMaterial*>(MaterPtr);
				if(NonlinMater != nullptr)
				{
					double chi_init = NonlinMater->GetInitialChi_ELF_Style();
					if(chi_init > 0)
					{
						poly->CurrentChi = chi_init;
					}
				}
			}
		}
	}

	// B-field convergence tracking (same as LU solver, ELF mucal2)
	std::vector<double> OldBnorm(AmOfMainElem, 0.0);
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Added ELF-style chi initialization to HACApK solver!')
else:
    print('Old code not found - checking if already fixed')
    if 'GetInitialChi_ELF_Style' in content and 'Initialize CurrentChi with ELF-style' in content:
        # Check HACApK section specifically
        hacapk_start = content.find('AutoRelax_VariableDOF(double PrecOnMagnetiz, int MaxIterNumber, char MagnResetIsNotNeeded)')
        if hacapk_start > 0:
            hacapk_section = content[hacapk_start:hacapk_start+3000]
            if 'GetInitialChi_ELF_Style' in hacapk_section:
                print('Already has chi initialization in HACApK section')
            else:
                print('Chi init exists but not in HACApK - needs manual fix')
        else:
            print('HACApK section not found')
    else:
        print('Unknown state')
