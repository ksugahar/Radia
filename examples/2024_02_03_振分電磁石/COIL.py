from numpy import *
from scipy.spatial.transform import Rotation

class cCOIL:
	def __init__(self,type,I,x0,V,W,H,R,phi,tilt,L):
		self.type = type
		self.I = I
		M = array([[cos(tilt/180*pi),0,-sin(tilt/180*pi)],[0,1.0,0],[sin(tilt/180*pi),0,cos(tilt/180*pi)]])
		self.V = M@V
		rot = Rotation.from_matrix(self.V)
		self.rot = rot.as_euler('ZXZ')*(-180/pi)
		self.W = abs( cos(tilt/180*pi)*W + sin(tilt/180*pi)*H)
		self.H = abs(-sin(tilt/180*pi)*W + cos(tilt/180*pi)*H)
		if (type=='ARC'):
			self.R = R
			self.x0 = x0
			self.c = x0 - R*self.V[0,0:3]
			self.phi = phi
			M = array([[cos(phi/180*pi),sin(phi/180*pi),0],[-sin(phi/180*pi),cos(phi/180*pi),0],[0,0,1]])
			self.V1 = M@self.V
			self.x1 = self.c + R*self.V1[0,0:3]
		else:
			self.L = L
			self.x0 = x0
			self.V1 = self.V
			self.xa = self.W/2*self.V[0,0:3]
			self.x1 = x0 + L*self.V[1,0:3]
			self.xb = self.H/2*self.V[2,0:3]
			self.c  = (self.x0+self.x1)/2
