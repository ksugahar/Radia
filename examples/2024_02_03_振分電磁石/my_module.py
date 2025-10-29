import radia as rad
import csv
from itertools import accumulate

def chunks(lst, n):
	"""Yield successive n-sized chunks from a list called 'lst'."""
	for i in range(0, len(lst), n):
		yield lst[i:i + n]

def exportGeometryToVTK(obj, fileName='radia_Geometry'):
	'''
	Writes the geometry of RADIA object "obj" to file fileName.vtk for use in Paraview. The format is VTK legacy, because it's simple. The file consists of polygons only (no cells).	
	'''
	
	vtkData = rad.ObjDrwVTK(obj, 'Axes->False')

	lengths = vtkData['polygons']['lengths']
	nPoly = len(lengths)
	offsets = list(accumulate(lengths))
	offsets.insert(0, 0) # prepend list with a zero
	points = vtkData['polygons']['vertices']
	nPnts = int(len(points)/3)

	# format the points array to be floats rather than double
	points = [round(num, 8) for num in points]		
	# define the connectivity list
	conn = list(range(nPnts)) 
	# define colours array
	colors = vtkData['polygons']['colors']

	# pre-process the output lists to have chunkLength items per line
	chunkLength = 9 # this writes 9 numbers per line (9 is the number used in Paraview if data is saved as the VTK Legacy format)
	offsets = list(chunks(offsets, chunkLength))
	points = list(chunks(points, chunkLength))
	conn = list(chunks(conn, chunkLength))
	colors = list(chunks(colors, chunkLength))

	# write the data to file
	with open(fileName + ".vtk", "w", newline="") as f:
		f.write('# vtk DataFile Version 5.1\n')
		f.write('vtk output\nASCII\nDATASET POLYDATA\n')
		f.write('POINTS ' + str(nPnts) + ' float\n')
		writer = csv.writer(f, delimiter=" ")
		writer.writerows(points)
		f.write('\n')
		f.write('POLYGONS ' + str(nPoly+1) + ' ' + str(nPnts) + '\n')
		f.write('OFFSETS vtktypeint64\n')
		writer.writerows(offsets)
		f.write('CONNECTIVITY vtktypeint64\n')  
		writer.writerows(conn)
		f.write('\n')
		f.write('CELL_DATA ' + str(nPoly) + '\n')
		f.write('COLOR_SCALARS Radia_colours 3\n')
		writer.writerows(colors)
