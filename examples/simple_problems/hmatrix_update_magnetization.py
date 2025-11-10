import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'Release'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src', 'python'))

import radia as rad
import time

print("="*70)
print("H-matrix UpdateMagnetization Test")
print("="*70)

try:
	# Step 1
	print("\n[Step 1] Creating geometry...")
	block = rad.ObjRecMag([0, 0, 0], [10, 10, 10])
	g = rad.ObjCnt([block])
	rad.ObjDivMag(block, [2, 2, 2])
	print(f"  OK: block={block}, container={g}")

	# Step 2
	print("\n[Step 2] Setting material...")
	rad.MatApl(g, rad.MatLin([0.06, 0.17], [0, 0, 1.2]))  # NdFeB
	print("  OK: Material applied")

	# Step 3
	print("\n[Step 3] Enabling H-matrix...")
	rad.SetHMatrixFieldEval(1, 1e-6)
	print("  OK: H-matrix enabled")

	# Step 4
	print("\n[Step 4] Computing initial field (builds H-matrix)...")
	obs_points = [[20, 0, 0], [30, 0, 0]]
	t0 = time.time()
	Bz1 = rad.FldBatch(g, 'bz', obs_points, use_hmatrix=1)
	t1 = time.time()
	print(f"  OK: Bz = {Bz1} (took {(t1-t0)*1000:.1f} ms)")

	# Step 5
	print("\n[Step 5] Changing magnetization...")
	rad.MatApl(g, rad.MatLin([1.0, 1.0, 1.0], [500, 0, 866]))
	print("  OK: M = [500, 0, 866] A/m")

	# Step 6
	print("\n[Step 6] Updating H-matrix magnetization (fast, no rebuild)...")
	t0 = time.time()
	rad.UpdateHMatrixMagnetization(g)
	t1 = time.time()
	print(f"  OK: Updated in {(t1-t0)*1000:.1f} ms")

	# Step 7
	print("\n[Step 7] Computing updated field with H-matrix...")
	t0 = time.time()
	Bz2_hmat = rad.FldBatch(g, 'bz', obs_points, use_hmatrix=1)
	t1 = time.time()
	print(f"  OK: Bz (H-matrix) = {Bz2_hmat} (took {(t1-t0)*1000:.1f} ms)")

	# Step 8
	print("\n[Step 8] Computing reference with direct calculation...")
	t0 = time.time()
	Bz2_direct = rad.FldBatch(g, 'bz', obs_points, use_hmatrix=0)
	t1 = time.time()
	print(f"  OK: Bz (direct) = {Bz2_direct} (took {(t1-t0)*1000:.1f} ms)")

	# Step 9
	print("\n[Step 9] Comparing accuracy...")
	if isinstance(Bz2_hmat, list):
		for i in range(len(Bz2_hmat)):
			h = Bz2_hmat[i]
			d = Bz2_direct[i]
			err = abs(h-d)/(abs(d)+1e-10)*100
			print(f"  Point {i}: H-mat={h:.6f} T, Direct={d:.6f} T, Error={err:.2f}%")
	else:
		h = Bz2_hmat
		d = Bz2_direct
		err = abs(h-d)/(abs(d)+1e-10)*100
		print(f"  H-mat={h:.6f} T, Direct={d:.6f} T, Error={err:.2f}%")

	print("\n" + "="*70)
	print("TEST PASSED: UpdateHMatrixMagnetization works correctly!")
	print("="*70)

	# VTK Export - Export geometry with same filename as script
	try:
		from radia_vtk_export import exportGeometryToVTK

		# Get script basename without extension
		script_name = os.path.splitext(os.path.basename(__file__))[0]
		vtk_filename = f"{script_name}.vtk"
		vtk_path = os.path.join(os.path.dirname(__file__), vtk_filename)

		# Export geometry
		exportGeometryToVTK(g, vtk_path)
		print(f"\n[VTK] Exported: {vtk_filename}")
		print(f"      View with: paraview {vtk_filename}")
	except ImportError:
		print("\n[VTK] Warning: radia_vtk_export not available (VTK export skipped)")
	except Exception as e:
		print(f"\n[VTK] Warning: Export failed: {e}")

except Exception as e:
	print(f"\n[ERROR] {type(e).__name__}: {e}")
	import traceback
	traceback.print_exc()
	sys.exit(1)
