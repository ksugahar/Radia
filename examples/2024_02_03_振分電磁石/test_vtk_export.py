#!/usr/bin/env python
"""
Test if exportGeometryToVTK and rad.ObjDrwVTK work
"""

import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'dist'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'build', 'lib', 'Release'))

import radia as rad

print("="*60)
print("Test: VTK Export Functionality")
print("="*60)

# Check if ObjDrwVTK exists
if hasattr(rad, 'ObjDrwVTK'):
	print('\nSUCCESS: rad.ObjDrwVTK is implemented')

	# Test with a simple object
	g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
	print(f'Created test object: ID {g}')

	try:
		vtk_data = rad.ObjDrwVTK(g, 'Axes->False')
		print(f'SUCCESS: ObjDrwVTK works!')
		print(f'  Data keys: {list(vtk_data.keys())}')

		if 'polygons' in vtk_data:
			poly = vtk_data['polygons']
			print(f'  Polygons keys: {list(poly.keys())}')
			print(f'  Number of vertices: {len(poly["vertices"])}')
			print(f'  Number of polygons: {len(poly["lengths"])}')

			# Test exportGeometryToVTK
			print('\nTesting exportGeometryToVTK...')
			from my_module import exportGeometryToVTK
			exportGeometryToVTK(g, 'test_export')

			if os.path.exists('test_export.vtk'):
				size = os.path.getsize('test_export.vtk')
				print(f'SUCCESS: VTK file created: test_export.vtk ({size} bytes)')

				# Show first few lines
				with open('test_export.vtk', 'r') as f:
					lines = f.readlines()[:10]
				print(f'  First 10 lines:')
				for line in lines:
					print(f'    {line.rstrip()}')
			else:
				print('ERROR: VTK file was not created')

	except Exception as e:
		print(f'ERROR: ObjDrwVTK failed: {e}')
		import traceback
		traceback.print_exc()
else:
	print('ERROR: rad.ObjDrwVTK is NOT implemented')

print()
print("="*60)
