import sys
import glob
import os
import win32com.client
import numpy as np
import array

PPT = win32com.client.Dispatch("PowerPoint.Application")
PPT.Visible = 1
PPT.Presentations.Add();
PPT.ActivePresentation.PageSetup.SlideSize = 7;
PPT.ActivePresentation.PageSetup.SlideWidth = 700;
PPT.ActivePresentation.PageSetup.SlideHeight = 400;

PPT.ActivePresentation.PageSetup.FirstSlideNumber = 1;
layout = PPT.ActivePresentation.SlideMaster.CustomLayouts.Item(7);
newSlide = PPT.ActivePresentation.Slides.AddSlide(1,layout);

for n in range(24):
	x0 = 50+25*n
	y0 = 200-150
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0-10, y0-10,20,20)
	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
	a.Line.ForeColor.RGB = 0
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0- 4, y0- 4,8,8)
	a.Fill.ForeColor.RGB = 0

for n in range(15):
	x0 = 50+40*n
	y0 = 200-50
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0-10,y0-10,20,20)
	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
	a.Line.ForeColor.RGB = 0
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddLine(x0-7,y0-7,x0+7,y0+7)
	a.Line.ForeColor.RGB = 0
	a.Line.Weight = 2
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddLine(x0-7,y0+7,x0+7,y0-7)
	a.Line.ForeColor.RGB = 0
	a.Line.Weight = 2

for n in range(15):
	x0 = 50+40*n
	y0 = 200+50
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0-10, y0-10,20,20)
	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
	a.Line.ForeColor.RGB = 0
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0- 4, y0- 4,8,8)
	a.Fill.ForeColor.RGB = 0

for n in range(24):
	x0 = 50+25*n
	y0 = 200+150
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddShape(9,x0-10,y0-10,20,20)
	a.Fill.ForeColor.RGB = 255*255*255+3*255*255+3*255
	a.Line.ForeColor.RGB = 0
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddLine(x0-7,y0-7,x0+7,y0+7)
	a.Line.ForeColor.RGB = 0
	a.Line.Weight = 2
	a = PPT.ActiveWindow.Selection.SlideRange.Shapes.AddLine(x0-7,y0+7,x0+7,y0-7)
	a.Line.ForeColor.RGB = 0
	a.Line.Weight = 2

ncoil = len(PPT.ActiveWindow.Selection.SlideRange.Shapes)

for n in range(24):
	x0 = 50+25*n
	y0 = 200
	x =  45*np.cos(np.linspace(-np.pi/2,np.pi/2-0.2,11))+x0
	y = 150*np.sin(np.linspace(-np.pi/2,np.pi/2-0.2,11))+y0
	myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(1,x[0],y[0]);
	for n in range(len(x)-1):
		myFreeForm.AddNodes(1,0,x[n+1],y[n+1])
	a = myFreeForm.ConvertToShape()

for n in range(15):
	x0 = 50+40*n
	y0 = 200
	x =  30*np.cos(np.linspace(np.pi/2,-np.pi/2+0.4,11))+x0
	y =  50*np.sin(np.linspace(np.pi/2,-np.pi/2+0.4,11))+y0
	myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(1,x[0],y[0]);
	for n in range(len(x)-1):
		myFreeForm.AddNodes(1,0,x[n+1],y[n+1])
	a = myFreeForm.ConvertToShape()

for n in range(ncoil,len(PPT.ActiveWindow.Selection.SlideRange.Shapes)):
	PPT.ActiveWindow.Selection.SlideRange.Shapes(n+1).Line.ForeColor.RGB = 0
	PPT.ActiveWindow.Selection.SlideRange.Shapes(n+1).Line.EndArrowheadStyle  = 2
	PPT.ActiveWindow.Selection.SlideRange.Shapes(n+1).Line.EndArrowheadLength = 3
	PPT.ActiveWindow.Selection.SlideRange.Shapes(n+1).Line.EndArrowheadWidth  = 2

ncoil = len(PPT.ActiveWindow.Selection.SlideRange.Shapes)

for n in range(1):
	x0 = 50+25*n
	y0 = 200
	x = -45*np.cos(np.linspace(-np.pi/2,np.pi/2,11))+x0
	y = 150*np.sin(np.linspace(-np.pi/2,np.pi/2,11))+y0
	myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(1,x[0],y[0]);
	for n in range(len(x)-1):
		myFreeForm.AddNodes(1,0,x[n+1],y[n+1])
	a = myFreeForm.ConvertToShape()

for n in range(1):
	x0 = 50+40*n
	y0 = 200
	x = -30*np.cos(np.linspace(np.pi/2,-np.pi/2,11))+x0
	y =  50*np.sin(np.linspace(np.pi/2,-np.pi/2,11))+y0
	myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(1,x[0],y[0]);
	for n in range(len(x)-1):
		myFreeForm.AddNodes(1,0,x[n+1],y[n+1])
	a = myFreeForm.ConvertToShape()

for n in range(ncoil,len(PPT.ActiveWindow.Selection.SlideRange.Shapes)):
	PPT.ActiveWindow.Selection.SlideRange.Shapes(n+1).Line.ForeColor.RGB = 0

#myFreeForm = PPT.ActiveWindow.Selection.SlideRange.Shapes.BuildFreeform(0,50,200);
#x = linspace(50,450,401)
#y = 50*sin(4*pi*(x-x[0])/100)+200
#for n in range(len(x)):
#	myFreeForm.AddNodes(0,0,x[n],y[n])
#myFreeForm.ConvertToShape()

PPT.ActivePresentation.SaveAs(os.getcwd() + u"test.pptx")

