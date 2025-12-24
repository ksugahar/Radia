#!/usr/bin/env python
"""
Complete rewrite of HACApK solver to match LU solver structure.

Changes:
1. Add OldBnorm array for B-field convergence tracking
2. Use ComputeChiDualMethod for chi update (same as LU)
3. Use B-field convergence criterion (same as LU/ELF)
4. Proper H_new computation from solved M
"""

filepath = r's:\Radia\01_GitHub\src\core\rad_relaxation_methods.cpp'

with open(filepath, 'r', encoding='utf-8') as f:
    content = f.read()

# Find and replace the entire HACApK nonlinear loop section
# From "// Outer nonlinear iteration" to "IntrctPtr->RelaxStatusParam.MisfitM"

old_code = '''	// Outer nonlinear iteration
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}

		// Update H from M for nonlinear materials
		// For 6DOF MSC hexahedra: compute effective M from sigma, then H = M/chi
		if(outerIter > 0)
		{
			for(int elem = 0; elem < AmOfMainElem; elem++)
			{
				int dof = IntrctPtr->GetElementDOF(elem);
				int offset = IntrctPtr->GetElementDOFOffset(elem);
				radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
				radTMaterial* MaterPtr = (radTMaterial*)(g3dRelaxPtr->MaterHandle.rep);

				if(dof == 3)
				{
					// 3DOF case: H = M / chi
					TVector3d H_est = IntrctPtr->NewFieldArray[elem];
					TMatrix3d KsiTensor;
					TVector3d MrVect;
					MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

					double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
					if(chi < 1.0e-6) chi = 1.0e-6;

					IntrctPtr->NewFieldArray[elem].x = FlatMagn[offset + 0] / chi;
					IntrctPtr->NewFieldArray[elem].y = FlatMagn[offset + 1] / chi;
					IntrctPtr->NewFieldArray[elem].z = FlatMagn[offset + 2] / chi;
				}
				else if(dof == 6)
				{
					// 6DOF MSC hexahedra: compute effective M from sigma, then H = M/chi
					radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
					if(poly && poly->Use6DOF_MSC)
					{
						// Compute effective magnetization from sigma using weighted average
						double Mx = 0.0, My = 0.0, Mz = 0.0;
						double wx = 0.0, wy = 0.0, wz = 0.0;
						for(int face = 0; face < 6; face++)
						{
							double sigma = FlatMagn[offset + face];
							TVector3d& n = poly->FaceNormal[face];
							double nx2 = n.x * n.x;
							double ny2 = n.y * n.y;
							double nz2 = n.z * n.z;
							Mx += sigma * n.x;
							My += sigma * n.y;
							Mz += sigma * n.z;
							wx += nx2;
							wy += ny2;
							wz += nz2;
						}
						double M_eff_x = (wx > 1.0e-10) ? Mx / wx : 0.0;
						double M_eff_y = (wy > 1.0e-10) ? My / wy : 0.0;
						double M_eff_z = (wz > 1.0e-10) ? Mz / wz : 0.0;

						// Get chi from current H estimate
						TVector3d H_est = IntrctPtr->NewFieldArray[elem];
						TMatrix3d KsiTensor;
						TVector3d MrVect;
						MaterPtr->DefineInstantKsiTensor(H_est, KsiTensor, MrVect);

						double chi = (KsiTensor.Str0.x + KsiTensor.Str1.y + KsiTensor.Str2.z) / 3.0;
						if(chi < 1.0e-6) chi = 1.0e-6;

						// Update H = M_eff / chi
						IntrctPtr->NewFieldArray[elem].x = M_eff_x / chi;
						IntrctPtr->NewFieldArray[elem].y = M_eff_y / chi;
						IntrctPtr->NewFieldArray[elem].z = M_eff_z / chi;
					}
				}
			}
		}

		// Solve with BiCGSTAB using H-matrix
		double residual = 0.0;
		const double bicg_tol = 1.0e-6;
		int n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(totalDOF, bicg_tol, MaxIterNumber - totalIterCount, residual);
		totalIterCount += n_iter;

		// Update element magnetization
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double diff = FlatMagn[i] - OldMagn[i];
			M_diff_sq += diff * diff;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

		// Sync magnetization to element objects
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];

			if(dof == 3)
			{
				g3dRelaxPtr->Magn.x = FlatMagn[offset];
				g3dRelaxPtr->Magn.y = FlatMagn[offset + 1];
				g3dRelaxPtr->Magn.z = FlatMagn[offset + 2];
			}
			else if(dof == 6)
			{
				radTPolyhedron* poly = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
				if(poly && poly->Use6DOF_MSC)
				{
					// Store sigma values
					for(int k = 0; k < 6; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma using weighted average
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < 6; face++)
					{
						double sigma = poly->Sigma[face];
						TVector3d& n = poly->FaceNormal[face];
						double nx2 = n.x * n.x;
						double ny2 = n.y * n.y;
						double nz2 = n.z * n.z;
						Mx += sigma * n.x;
						My += sigma * n.y;
						Mz += sigma * n.z;
						wx += nx2;
						wy += ny2;
						wz += nz2;
					}
					if(wx > 1.0e-10) g3dRelaxPtr->Magn.x = Mx / wx;
					if(wy > 1.0e-10) g3dRelaxPtr->Magn.y = My / wy;
					if(wz > 1.0e-10) g3dRelaxPtr->Magn.z = Mz / wz;

					// FIX (2025-12-21): Update NewFieldArray with H = M / chi
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
					}
				}
			}
		}

		double rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		MisfitE2 = rel_change * rel_change;

		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		if(radYield.Check() == 0) return outerIter;
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);'''

new_code = '''	// Cache polyhedron pointers for fast access (same as LU solver)
	std::vector<radTPolyhedron*> polyCache(AmOfMainElem, nullptr);
	for(int elem = 0; elem < AmOfMainElem; elem++)
	{
		radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];
		polyCache[elem] = dynamic_cast<radTPolyhedron*>(g3dRelaxPtr);
	}

	// B-field convergence tracking (same as LU solver, ELF mucal2)
	std::vector<double> OldBnorm(AmOfMainElem, 0.0);
	const double MU_0 = 4.0 * 3.14159265358979323846 * 1.0e-7;

	// Outer nonlinear iteration (rewritten to match LU solver structure)
	for(outerIter = 0; outerIter < MaxIterNumber; outerIter++)
	{
		// Store old values
		for(int i = 0; i < totalDOF; i++)
		{
			OldMagn[i] = FlatMagn[i];
		}

		// Store old B norm for convergence check (same as LU solver)
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC)
				{
					double chi = poly->CurrentChi;
					if(chi < 1.0e-6) chi = 1.0e-6;
					TVector3d& M = poly->Magn;
					TVector3d H(M.x / chi, M.y / chi, M.z / chi);
					TVector3d B(MU_0 * (H.x + M.x), MU_0 * (H.y + M.y), MU_0 * (H.z + M.z));
					OldBnorm[elem] = std::sqrt(B.x*B.x + B.y*B.y + B.z*B.z);
				}
			}
		}

		// Solve with BiCGSTAB using H-matrix
		double residual = 0.0;
		const double bicg_tol = 1.0e-6;
		int n_iter = SolveBiCGSTAB_HMatrix_VariableDOF(totalDOF, bicg_tol, MaxIterNumber - totalIterCount, residual);
		totalIterCount += n_iter;

		// Update element magnetization from flat array
		double M_diff_sq = 0.0;
		double M_norm_sq = 0.0;
		for(int i = 0; i < totalDOF; i++)
		{
			double diff = FlatMagn[i] - OldMagn[i];
			M_diff_sq += diff * diff;
			M_norm_sq += FlatMagn[i] * FlatMagn[i];
		}

		// Sync magnetization to element objects and compute H_new
		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			int offset = IntrctPtr->GetElementDOFOffset(elem);
			radTg3dRelax* g3dRelaxPtr = IntrctPtr->g3dRelaxPtrVect[elem];

			if(dof == 3)
			{
				g3dRelaxPtr->Magn.x = FlatMagn[offset];
				g3dRelaxPtr->Magn.y = FlatMagn[offset + 1];
				g3dRelaxPtr->Magn.z = FlatMagn[offset + 2];
			}
			else if(dof == 6)
			{
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC)
				{
					// Store sigma values
					for(int k = 0; k < 6; k++)
					{
						poly->Sigma[k] = FlatMagn[offset + k];
					}

					// Compute effective magnetization from sigma (same as LU solver)
					double Mx = 0.0, My = 0.0, Mz = 0.0;
					double wx = 0.0, wy = 0.0, wz = 0.0;
					for(int face = 0; face < 6; face++)
					{
						double sigma = poly->Sigma[face];
						TVector3d& n = poly->FaceNormal[face];
						double nx2 = n.x * n.x;
						double ny2 = n.y * n.y;
						double nz2 = n.z * n.z;
						Mx += sigma * n.x;
						My += sigma * n.y;
						Mz += sigma * n.z;
						wx += nx2;
						wy += ny2;
						wz += nz2;
					}
					if(wx > 1.0e-10) poly->Magn.x = Mx / wx;
					if(wy > 1.0e-10) poly->Magn.y = My / wy;
					if(wz > 1.0e-10) poly->Magn.z = Mz / wz;

					// Compute H_new = M / chi_current (same as LU solver lines 1657-1668)
					if(IntrctPtr->NewFieldArray != nullptr)
					{
						double chi_used = poly->CurrentChi;
						if(chi_used < 1.0e-6) chi_used = 1.0e-6;
						IntrctPtr->NewFieldArray[elem].x = poly->Magn.x / chi_used;
						IntrctPtr->NewFieldArray[elem].y = poly->Magn.y / chi_used;
						IntrctPtr->NewFieldArray[elem].z = poly->Magn.z / chi_used;
					}
				}
			}
		}

		// Compute convergence and update chi (same structure as LU solver lines 1686-1793)
		double max_B_rel_change = 0.0;
		bool has_6dof_elements = false;

		for(int elem = 0; elem < AmOfMainElem; elem++)
		{
			int dof = IntrctPtr->GetElementDOF(elem);
			if(dof == 6)
			{
				has_6dof_elements = true;
				radTPolyhedron* poly = polyCache[elem];
				if(poly && poly->Use6DOF_MSC && IntrctPtr->NewFieldArray != nullptr)
				{
					TVector3d H_new = IntrctPtr->NewFieldArray[elem];
					radTMaterial* MaterPtr = (radTMaterial*)(IntrctPtr->g3dRelaxPtrVect[elem]->MaterHandle.rep);

					// Get chi used for this iteration's matrix (same as LU line 1699)
					double chi_matrix = poly->CurrentChi;
					double mu_old = chi_matrix + 1.0;

					// Compute H magnitude for chi update
					double H_mag = std::sqrt(H_new.x*H_new.x + H_new.y*H_new.y + H_new.z*H_new.z);

					// Use ELF-style dual-method chi update (same as LU lines 1708-1736)
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

					// Update chi for next iteration (same as LU line 1739)
					poly->CurrentChi = chi_new;

					// B-field convergence (same as LU lines 1741-1765)
					TVector3d& M_new = poly->Magn;
					double chi_for_B = chi_new;
					if(chi_for_B < 1.0e-6) chi_for_B = 1.0e-6;
					TVector3d H_for_B(M_new.x / chi_for_B, M_new.y / chi_for_B, M_new.z / chi_for_B);
					TVector3d B_new_vec(MU_0 * (H_for_B.x + M_new.x),
					                    MU_0 * (H_for_B.y + M_new.y),
					                    MU_0 * (H_for_B.z + M_new.z));
					double B_new_norm = std::sqrt(B_new_vec.x*B_new_vec.x + B_new_vec.y*B_new_vec.y + B_new_vec.z*B_new_vec.z);

					// Get B_sat from BH curve (same as LU lines 1753-1759)
					double B_sat = 1.0;
					if(NonlinMater != nullptr)
					{
						B_sat = NonlinMater->GetBsaturation();
						if(B_sat < 1.0e-10) B_sat = 1.0;
					}

					// B-field convergence: |B_new - B_old| / B_sat (same as LU lines 1761-1765)
					double B_old_norm = OldBnorm[elem];
					double B_rel_change = std::fabs(B_new_norm - B_old_norm) / B_sat;
					if(B_rel_change > max_B_rel_change)
						max_B_rel_change = B_rel_change;
				}
			}
		}

		// Convergence criterion (same as LU lines 1779-1794)
		double rel_change;
		if(has_6dof_elements)
		{
			// For 6DOF MSC: use ELF-style B-field change (mucal2)
			rel_change = max_B_rel_change;
		}
		else
		{
			// For 3DOF MMM: use M change
			rel_change = (M_norm_sq > 1.0e-30) ? std::sqrt(M_diff_sq / M_norm_sq) : std::sqrt(M_diff_sq);
		}
		MisfitE2 = rel_change * rel_change;

		if(rel_change <= PrecOnMagnetiz)
		{
			outerIter++;
			break;
		}

		if(radYield.Check() == 0) return outerIter;
	}

	IntrctPtr->RelaxStatusParam.MisfitM = std::sqrt(MisfitE2);'''

if old_code in content:
    content = content.replace(old_code, new_code)
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    print('Successfully rewrote HACApK solver to match LU structure!')
else:
    print('Old code not found - checking current state')
    # Try to find a substring
    if '// Outer nonlinear iteration' in content and 'SolveBiCGSTAB_HMatrix_VariableDOF' in content:
        print('HACApK solver exists but code structure different')
        print('Please check manually')
    else:
        print('HACApK solver not found')
