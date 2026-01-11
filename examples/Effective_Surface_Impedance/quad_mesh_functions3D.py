from netgen.meshing import *
from netgen.csg import *

import math as math

def MacroElement3D(delta,numElements=1,manual=False,manualVec=[0.5,1]):
	return {'delta':delta,'numElements':numElements,'manual':manual,'manualVec':manualVec}
	
def CuboidDomain(startX,endX,startY,endY,startZ,endZ,domainNumber=1):
	return {'startX':startX,'endX':endX,'startY':startY,'endY':endY,'startZ':startZ,'endZ':endZ,'domainNumber':domainNumber}
	
def MetaInfo3D(originX=0,originY=0,originZ=0,centered=False,scale=1,allBoundariesNumbered=False,refinements=1,maxh="Inf",maxRatio="Inf"):
	return {'originX':originX,'originY':originY,'originZ':originZ,'centered':centered,'scale':scale,'allBoundariesNumbered':allBoundariesNumbered,'refinements':refinements,'maxh':maxh,'maxRatio':maxRatio}
	
def BoundaryCuboid(startX,endX,startY,endY,startZ,endZ,left=-1,right=-1,top=-1,bottom=-1,front=-1,back=-1):
	return {'startX':startX,'endX':endX,'startY':startY,'endY':endY,'startZ':startZ,'endZ':endZ,'left':left,'right':right,'top':top,'bottom':bottom,'front':front,'back':back}
	
def OuterBoundary3D(left=1,right=1,top=1,bottom=1,front=1,back=1):
	return {'left':left,'right':right,'top':top,'bottom':bottom,'front':front,'back':back}
	
def BoundaryInfo3D(outer=OuterBoundary3D(),rectangles=[]):
	return {'outer':outer,'rectangles':rectangles}
	
def makeCuboidMesh(name="cuboidmesh",metaInfo=MetaInfo3D(),macroElementsX=[],macroElementsY=[],macroElementsZ=[],cuboidDomains=[],boundaryInfo=BoundaryInfo3D(),materialNames=[],boundaryNames=[]):
	mesh = Mesh()
	mesh.dim = 3
	
	originX = metaInfo['originX']
	originY = metaInfo['originY']
	originZ = metaInfo['originZ']
	centered = metaInfo['centered']
	scale = metaInfo['scale']
	allBoundariesNumbered = metaInfo['allBoundariesNumbered']
	refinements = metaInfo['refinements']
	maxh = metaInfo['maxh']
	maxRatio = metaInfo['maxRatio']	
	
	
	
	#apply scale
	
	for i in range(len(macroElementsX)):
		macroElementsX[i]['delta'] *= scale
		
	for i in range(len(macroElementsY)):
		macroElementsY[i]['delta'] *= scale
		
	for i in range(len(macroElementsZ)):
		macroElementsZ[i]['delta'] *= scale
	
	
	originX *= scale
	originY *= scale
	originZ *= scale


		
		
	#ensure maxh
	
	if maxh != "Inf":
		for i in range(len(macroElementsX)):
			ME=macroElementsX[i]
			if ME['delta']/ME['numElements'] > maxh:
				ME['numElements'] = math.ceil(ME['delta']/maxh)
					
		for i in range(len(macroElementsY)):
			ME=macroElementsY[i]
			if ME['delta']/ME['numElements'] > maxh:
				ME['numElements'] = math.ceil(ME['delta']/maxh)
				
		for i in range(len(macroElementsZ)):
			ME=macroElementsZ[i]
			if ME['delta']/ME['numElements'] > maxh:
				ME['numElements'] = math.ceil(ME['delta']/maxh)
					
					
	#ensure ratio
	
	if maxRatio != "Inf":
		def getLength(ME):
			return ME['delta']/ME['numElements']
		
		def doubleElements(ME):
			ME['numElements'] *= 2
	
	
		allRatiosOK = False
		loopCounter = 0
		while not allRatiosOK:
			loopCounter += 1
			if loopCounter in [10**i for i in range(2,6)]:
				print()
				print()
				print(loopCounter,"loops in ensure ratio routine, likely caused by the loop being infinite. If no result is reached, try again with a higher ratio tolerance")
				print()
				print()
			allRatiosOK = True
			for i in range(len(macroElementsX)):
				for j in range(len(macroElementsY)):
					for k in range(len(macroElementsZ)):
						if getLength(macroElementsX[i]) / getLength(macroElementsY[j]) < 1 / maxRatio:
							macroElementsY[j]['numElements'] *= 2
							allRatiosOK = False
						if getLength(macroElementsX[i]) / getLength(macroElementsY[j]) > maxRatio:
							macroElementsX[j]['numElements'] *= 2
							allRatiosOK = False
						if getLength(macroElementsX[i]) / getLength(macroElementsZ[k]) < 1 / maxRatio:
							macroElementsZ[j]['numElements'] *= 2
							allRatiosOK = False
						if getLength(macroElementsX[i]) / getLength(macroElementsZ[k]) > maxRatio:
							macroElementsX[j]['numElements'] *= 2
							allRatiosOK = False
						if getLength(macroElementsY[j]) / getLength(macroElementsZ[k]) < 1 / maxRatio:
							macroElementsZ[k]['numElements'] *= 2
							allRatiosOK = False
						if getLength(macroElementsY[j]) / getLength(macroElementsZ[k]) > maxRatio:
							macroElementsY[j]['numElements'] *= 2
							allRatiosOK = False
	
	
	#apply refinements
	
	for i in range(len(macroElementsX)):
		macroElementsX[i]['numElements'] *= refinements
			
	for i in range(len(macroElementsY)):
		macroElementsY[i]['numElements'] *= refinements
	
	for i in range(len(macroElementsZ)):
		macroElementsZ[i]['numElements'] *= refinements
		
		
		
	
	#calculate x coordinates of grid points; shift center if necessary
	
	totalX = 0
	for i in range(len(macroElementsX)):
		totalX += macroElementsX[i]['delta']
	
	currentX = 0
	if centered:
		currentX = -totalX/2
	currentX += originX
	
	x_coords=[currentX]
	
	for i in range(len(macroElementsX)):
		ME=macroElementsX[i]
		delta=ME['delta']
		numElements = ME['numElements']
		if ME['manual']:
			manualVec=ME['manualVec']
			for j in range(len(ME['manualVec'])):
				x_coords.append(currentX+delta*manualVec[j])
		else:
			for j in range(numElements):
				x_coords.append(currentX+(j+1)*delta/(numElements))
		currentX += delta
				
	
	
	
	#calculate y coordinates of grid points; shift center if necessary
	
	totalY= 0
	for i in range(len(macroElementsY)):
		totalY += macroElementsY[i]['delta']
	
	currentY = 0
	if centered:
		currentY = -totalY/2
	currentY += originY
	
	y_coords=[currentY]
	
	for i in range(len(macroElementsY)):
		ME=macroElementsY[i]
		delta=ME['delta']
		numElements = ME['numElements']
		if ME['manual']:
			manualVec=ME['manualVec']
			for j in range(len(ME['manualVec'])):
				y_coords.append(currentY+delta*manualVec[j])
		else:
			for j in range(numElements):
				y_coords.append(currentY+(j+1)*delta/(numElements))
		currentY += delta
		
		
		
	#calculate z coordinates of grid points; shift center if necessary
	
	totalZ= 0
	for i in range(len(macroElementsZ)):
		totalZ += macroElementsZ[i]['delta']
	
	currentZ = 0
	if centered:
		currentZ = -totalZ/2
	currentZ += originZ
	
	z_coords=[currentZ]
	
	for i in range(len(macroElementsZ)):
		ME=macroElementsZ[i]
		delta=ME['delta']
		numElements = ME['numElements']
		if ME['manual']:
			manualVec=ME['manualVec']
			for j in range(len(ME['manualVec'])):
				z_coords.append(currentZ+delta*manualVec[j])
		else:
			for j in range(numElements):
				z_coords.append(currentZ+(j+1)*delta/(numElements))
		currentZ += delta
			
			
	#add mesh points to the mesh
			
	NEX=len(x_coords)
	NEY=len(y_coords)
	NEZ=len(z_coords)
	
	pids = []
	for k in range(NEZ):
		for j in range(NEY):
			for i in range(NEX):
				pids.append (mesh.Add (MeshPoint(Pnt(x_coords[i], y_coords[j], z_coords[k]))))
	
	
	
	
	
	
	#make a table of domain numbers in each crossed MacroElement
	
	domainTable = [[[1 for l in range(len(macroElementsZ))] for j in range(len(macroElementsY))] for i in range(len(macroElementsX))]
	
	for c in range(len(cuboidDomains)):
		RD=cuboidDomains[c]
		for i in range(RD['endX']-RD['startX']+1):
			for j in range(RD['endY']-RD['startY']+1):
				for k in range(RD['endZ']-RD['startZ']+1):
					domainTable[i+RD['startX']-1][j+RD['startY']-1][k+RD['startZ']-1]= RD['domainNumber']
					
					
	#caculate the bottom left point in each crossed MacroElement
					
	bottomleftPoints = []
	
	currentPoint = 0
	startThisColumn = 0
	for i in range(len(macroElementsX)):
		bottomleftPoints.append([])
		currentPoint = startThisColumn
		startThisHeight = currentPoint
		
		for j in range(len(macroElementsY)):
			bottomleftPoints[i].append([])
			currentPoint = startThisHeight
		
			for k in range(len(macroElementsZ)):
				bottomleftPoints[i][j].append(currentPoint)
				currentPoint += macroElementsZ[k]['numElements']*NEX*NEY
				
			startThisHeight += macroElementsY[j]['numElements'] * NEX
			
		startThisColumn += macroElementsX[i]['numElements']
		
	
	#add the quad elements to the mesh
	
	for i in range(len(macroElementsX)):
		for j in range(len(macroElementsY)):
			for k in range(len(macroElementsZ)):
				for ix in range(macroElementsX[i]['numElements']):
					for iy in range(macroElementsY[j]['numElements']):
						for iz in range(macroElementsZ[k]['numElements']):
							base = bottomleftPoints[i][j][k] + ix + iy*NEX + iz*NEX*NEY
							pnum = [base,base+1,base+NEX+1,base+NEX,base+NEX*NEY,base+1+NEX*NEY,base+NEX+1+NEX*NEY,base+NEX+NEX*NEY]
							elpids = [pids[p] for p in pnum]
							mesh.Add (Element3D(domainTable[i][j][k],elpids))
					
			
	
	#make a table of the boundary conditions on each crossed MacroElement
	
	boundaryTable = [[[{"left":-1,"right":-1,"top":-1,"bottom":-1,"front":-1,"back":-1} for k in range(len(macroElementsZ))] for j in range(len(macroElementsY))] for i in range(len(macroElementsX))]
	
	for i in range(len(macroElementsX)):
		for j in range(len(macroElementsY)):
			boundaryTable[i][j][0]["bottom"] = boundaryInfo['outer']['bottom']
			boundaryTable[i][j][len(macroElementsZ)-1]["top"] = boundaryInfo['outer']['top']
		
	for i in range(len(macroElementsX)):
		for k in range(len(macroElementsZ)):
			boundaryTable[i][0][k]["front"] = boundaryInfo['outer']['front']
			boundaryTable[i][len(macroElementsY)-1][k]["back"] = boundaryInfo['outer']['back']
			
	for j in range(len(macroElementsY)):
		for k in range(len(macroElementsZ)):
			boundaryTable[0][j][k]["left"] = boundaryInfo['outer']['left']
			boundaryTable[len(macroElementsX)-1][j][k]["right"] = boundaryInfo['outer']['right']
	
	boundRects=boundaryInfo['rectangles']
	
	for bnd in range(len(boundRects)):
		BR=boundRects[bnd]
		startX = BR['startX']
		endX = BR['endX']
		startY = BR['startY']
		endY = BR['endY']
		startZ = BR['startZ']
		endZ = BR['endZ']
		
		for i in range(endX-startX+1):
			for j in range(endY-startY+1):
				if BR['bottom'] > 0:
					boundaryTable[startX+i-1][startY+j-1][startZ-1]["bottom"] = BR['bottom']
					if startZ > 1:
						boundaryTable[startX+i-1][startY+j-1][startZ-2]["top"] = BR['bottom']
					
				if BR["top"] > 0:
					boundaryTable[startX+i-1][startY+j-1][endZ-1]["top"] = BR["top"]
					if endZ < len(macroElementsZ):
						boundaryTable[startX+i-1][startY+j-1][endZ]["bottom"] = BR["top"]
				
		
		for i in range(endX-startX+1):
			for k in range(endZ-startZ+1):
				if BR['front'] > 0:
					boundaryTable[startX+i-1][startY-1][startZ+k-1]["front"] = BR['front']
					if startY > 1:
						boundaryTable[startX+i-1][startY-2][startZ+k-1]["back"] = BR['front']
					
				if BR["back"] > 0:
					boundaryTable[startX+i-1][endY-1][startZ+k-1]["back"] = BR["back"]
					if endY < len(macroElementsY):
						boundaryTable[startX+i-1][endY][startZ+k-1]["front"] = BR["back"]
						
						
		for j in range(endY-startY+1):
			for k in range(endZ-startZ+1):
				if BR['left'] > 0:
					boundaryTable[startX-1][startY+j-1][startZ+k-1]["left"] = BR['left']
					if startX > 1:
						boundaryTable[startX-2][startY+j-1][startZ+k-1]["right"] = BR['left']
					
				if BR["right"] > 0:
					boundaryTable[endX-1][startY+j-1][startZ+k-1]["right"] = BR["right"]
					if endX < len(macroElementsX):
						boundaryTable[endX][startY+j-1][startZ+k-1]["left"] = BR["right"]
					
	
	
	
	maxNumber = 0
	maxNumber = max(maxNumber,boundaryInfo['outer']['right'],boundaryInfo['outer']['left'],boundaryInfo['outer']['top'],boundaryInfo['outer']['bottom'],boundaryInfo['outer']['front'],boundaryInfo['outer']['back'])
	for i in range(len(boundRects)):
		BR=boundRects[i]
		maxNumber = max(maxNumber,BR['left'],BR['right'],BR['top'],BR['bottom'],BR['front'],BR['back'])
	
	if allBoundariesNumbered:
		maxNumber += 1
		
		for i in range(len(macroElementsX)):
			for j in range(len(macroElementsY)):
				for k in range(len(macroElementsZ)):
					for key in boundaryTable[i][j][k]:
						if(boundaryTable[i][j][k][key] < 1):
							boundaryTable[i][j][k][key] = maxNumber
							

	#add necessary face descriptors
	
	for i in range(maxNumber):
		mesh.Add (FaceDescriptor(surfnr=i,domin=1,domout=0,bc=i+1))
							
	#add surface elements at boundaries with boundary condition
	
	for i in range(len(macroElementsX)):
		for j in range(len(macroElementsY)):
			for k in range(len(macroElementsZ)):
				if boundaryTable[i][j][k]["bottom"]>0:
					for ix in range(macroElementsX[i]['numElements']):
						for iy in range(macroElementsY[j]['numElements']):
							base = bottomleftPoints[i][j][k] + ix + iy*NEX
							pnum = [base,base+1,base+NEX+1,base+NEX]
							elpids = [pids[p] for p in pnum]
							mesh.Add (Element2D(boundaryTable[i][j][k]["bottom"],elpids))
							
				if boundaryTable[i][j][k]["front"]>0:
					for ix in range(macroElementsX[i]['numElements']):
						for iz in range(macroElementsZ[k]['numElements']):
							base = bottomleftPoints[i][j][k] + ix + iz*NEX*NEY
							pnum = [base,base+1,base+NEX*NEY+1,base+NEX*NEY]
							elpids = [pids[p] for p in pnum]
							mesh.Add (Element2D(boundaryTable[i][j][k]["front"],elpids))
							
				if boundaryTable[i][j][k]["left"]>0:
					for iy in range(macroElementsY[j]['numElements']):
						for iz in range(macroElementsZ[k]['numElements']):
							base = bottomleftPoints[i][j][k] + iy*NEX + iz*NEX*NEY
							pnum = [base,base+NEX,base+NEX*NEY+NEX,base+NEX*NEY]
							elpids = [pids[p] for p in pnum]
							mesh.Add (Element2D(boundaryTable[i][j][k]["left"],elpids))
							
	for i in range(len(macroElementsX)):
		for j in range(len(macroElementsY)):		
			for ix in range(macroElementsX[i]['numElements']):
				for iy in range(macroElementsY[j]['numElements']):
					base = bottomleftPoints[i][j][0] + NEX*NEY*(NEZ-1) + ix + iy*NEX
					pnum = [base,base+1,base+NEX+1,base+NEX]
					elpids = [pids[p] for p in pnum]
					mesh.Add (Element2D(boundaryTable[i][j][len(macroElementsZ)-1]["top"],elpids))
				
	for i in range(len(macroElementsX)):
		for k in range(len(macroElementsZ)):
			for ix in range(macroElementsX[i]['numElements']):
				for iz in range(macroElementsZ[k]['numElements']):
					base = bottomleftPoints[i][0][k] + NEX*(NEY-1) + ix + iz*NEX*NEY
					pnum = [base,base+1,base+NEX*NEY+1,base+NEX*NEY]
					elpids = [pids[p] for p in pnum]
					mesh.Add (Element2D(boundaryTable[i][len(macroElementsY)-1][k]["back"],elpids))
	
	for j in range(len(macroElementsY)):
		for k in range(len(macroElementsZ)):
			for iy in range(macroElementsY[j]['numElements']):
				for iz in range(macroElementsZ[k]['numElements']):
					base = bottomleftPoints[0][j][k] + NEX-1 + iy*NEX + iz*NEX*NEY
					pnum = [base,base+NEX,base+NEX*NEY+NEX,base+NEX*NEY]
					elpids = [pids[p] for p in pnum]
					mesh.Add (Element2D(boundaryTable[len(macroElementsX)-1][j][k]["right"],elpids))
	
	
	#material and boundary names
	
	
	for i in range(len(materialNames)):
		mesh.SetMaterial(i+1,materialNames[i])
		
	for i in range(len(boundaryNames)):
		mesh.SetBCName(i,boundaryNames[i])
	
	mesh.Save(name+".vol")